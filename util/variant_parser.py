"""
fetch_variants(): resolve a mixed list of variant identifiers to their
missense consequences using the Ensembl REST API (VEP).

Handles, in addition to all VEP-standard inputs:
  - dbSNP rsIDs                  e.g.  rs699, rs6265
  - UniProt acc + protein sub    e.g.  P00533 R132C
                                       P00533:R132C
                                       P00533/R132C
  - genomic SPDI-ish / region    e.g.  chr14 89993420 A/G

Design notes
------------
* Inputs are bucketed by type, each bucket hits the most appropriate VEP
  endpoint, all via POST batch calls (<=200 ids/call) so a few-thousand-row
  input stays to a few dozen HTTP requests.
* UniProt protein substitutions are not directly VEP-able. They are mapped
  UniProt-accession -> Ensembl translation (ENSP) via the xrefs endpoint
  (batched, cached), then submitted as protein-level HGVS (ENSP:p.Xaa#Yaa)
  to the VEP /hgvs endpoint, which back-maps to genomic coordinates.
* Only consequences containing "missense_variant" are returned.
* 429 rate-limit responses are honoured via Retry-After with backoff.

Requires: requests
"""

from __future__ import annotations

import re
import time
import itertools
from typing import Iterable, Iterator

import requests

SERVER = "https://rest.ensembl.org"
SPECIES = "homo_sapiens"
BATCH = 200                 # Ensembl POST limit
MAX_RETRIES = 5
SESSION = requests.Session()
SESSION.headers.update({"Content-Type": "application/json",
                        "Accept": "application/json"})

# Three-letter <-> one-letter amino acids for HGVS protein notation
_AA1to3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln",
    "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
    "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val", "U": "Sec", "*": "Ter", "X": "Xaa",
}

# UniProt accession: O/P/Q + digits/letters, or the newer 10-char form
_UNIPROT = r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}"
_PROTSUB = re.compile(
    rf"^\s*(?P<acc>{_UNIPROT})\s*[\s:/]\s*"
    r"(?P<ref>[A-Z*])(?P<pos>\d+)(?P<alt>[A-Z*])\s*$"
)
_RSID = re.compile(r"^\s*rs\d+\s*$", re.IGNORECASE)
# "chr14 89993420 A/G"  or "14 89993420 A/G" -> VEP default region format
_GENOMIC = re.compile(
    r"^\s*(?:chr)?(?P<chrom>[\dXYMT]+)\s+(?P<pos>\d+)\s+"
    r"(?P<ref>[ACGT-]+)\s*/\s*(?P<alt>[ACGT-]+)\s*$",
    re.IGNORECASE,
)


# --------------------------------------------------------------------------- #
# low-level HTTP with rate-limit handling
# --------------------------------------------------------------------------- #
def _post(endpoint: str, payload: dict, params: dict | None = None) -> list | dict:
    url = f"{SERVER}{endpoint}"
    for attempt in range(MAX_RETRIES):
        r = SESSION.post(url, json=payload, params=params or {})
        if r.status_code == 429:
            wait = float(r.headers.get("Retry-After", 1.0))
            time.sleep(wait + 0.1)
            continue
        if r.status_code in (500, 502, 503, 504):
            time.sleep(2 ** attempt)
            continue
        r.raise_for_status()
        return r.json()
    r.raise_for_status()
    return []


def _chunks(seq: list, n: int) -> Iterator[list]:
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


# --------------------------------------------------------------------------- #
# UniProt accession -> Ensembl ENSP translation id
# --------------------------------------------------------------------------- #
_ensp_cache: dict[str, str | None] = {}


def _map_uniprot_to_ensp(accs: Iterable[str]) -> dict[str, str | None]:
    """Resolve UniProt accessions to an Ensembl translation (ENSP) id.

    Uses GET /xrefs/symbol/<species>/<acc>?object_type=translation. There is no
    POST batch for xrefs, so this is per-accession but cached.
    """
    out: dict[str, str | None] = {}
    for acc in {a for a in accs if a not in _ensp_cache}:
        url = f"{SERVER}/xrefs/symbol/{SPECIES}/{acc}"
        for attempt in range(MAX_RETRIES):
            r = SESSION.get(url, params={"object_type": "translation"})
            if r.status_code == 429:
                time.sleep(float(r.headers.get("Retry-After", 1.0)) + 0.1)
                continue
            if r.status_code == 400:           # nothing found
                _ensp_cache[acc] = None
                break
            if r.status_code in (500, 502, 503, 504):
                time.sleep(2 ** attempt)
                continue
            r.raise_for_status()
            hits = r.json()
            _ensp_cache[acc] = hits[0]["id"] if hits else None
            break
        else:
            _ensp_cache[acc] = None
    for acc in accs:
        out[acc] = _ensp_cache.get(acc)
    return out


# --------------------------------------------------------------------------- #
# parsing / bucketing
# --------------------------------------------------------------------------- #
def _classify(raw: str):
    """Return (kind, normalized) where kind in
    {'rsid','protsub','genomic','vep'}."""
    s = raw.strip()
    if _RSID.match(s):
        return "rsid", s.lower()
    m = _PROTSUB.match(s)
    if m:
        return "protsub", m
    m = _GENOMIC.match(s)
    if m:
        gd = m.groupdict()
        ref, alt = gd["ref"].upper(), gd["alt"].upper()
        # VEP default region format: "chr start end allele strand"
        pos = int(gd["pos"])
        end = pos + len(ref) - 1
        return "genomic", f'{gd["chrom"]} {pos} {end} {ref}/{alt} 1'
    return "vep", s            # assume it is VEP-standard (HGVS, region, etc.)


def _protsub_to_hgvs(m: re.Match, ensp: str | None) -> str | None:
    if ensp is None:
        return None
    ref3 = _AA1to3.get(m.group("ref").upper())
    alt3 = _AA1to3.get(m.group("alt").upper())
    if not ref3 or not alt3:
        return None
    return f'{ensp}:p.{ref3}{m.group("pos")}{alt3}'


# --------------------------------------------------------------------------- #
# missense extraction
# --------------------------------------------------------------------------- #
def _extract_missense(vep_record: dict, original: str,
                      uniprot_variant_id: str | None = None) -> list[dict]:
    rows = []
    for tc in vep_record.get("transcript_consequences", []):
        cons = tc.get("consequence_terms", [])
        if "missense_variant" not in cons:
            continue
        rows.append({
            "input": original,
            "variant_id": vep_record.get("id"),
            "uniprot_variant_id": uniprot_variant_id,
            "location": vep_record.get("seq_region_name") and
                        f'{vep_record["seq_region_name"]}:'
                        f'{vep_record.get("start")}-{vep_record.get("end")}',
            "allele": vep_record.get("allele_string"),
            "gene_symbol": tc.get("gene_symbol"),
            "gene_id": tc.get("gene_id"),
            "transcript_id": tc.get("transcript_id"),
            "protein_id": tc.get("protein_id"),
            "amino_acids": tc.get("amino_acids"),
            "protein_position": tc.get("protein_start"),
            "hgvsp": tc.get("hgvsp"),
            "sift_prediction": tc.get("sift_prediction"),
            "polyphen_prediction": tc.get("polyphen_prediction"),
            "consequence_terms": cons,
        })
    return rows


# --------------------------------------------------------------------------- #
# public API
# --------------------------------------------------------------------------- #
def fetch_variants(identifiers: list[str], extra_vep_params: dict | None = None
                   ) -> dict:
    """
    Fetch missense variant consequences for a mixed list of identifiers.

    Returns
    -------
    {
      "missense":   [ {row}, ... ],   # one row per missense transcript_consequence
      "unmatched":  [ original_str, ... ],   # inputs VEP returned no missense for
      "unresolved": [ original_str, ... ],   # inputs that could not be parsed/mapped
    }
    """
    params = {"hgvs": 1, "sift": 1, "polyphen": 1, "canonical": 1}
    if extra_vep_params:
        params.update(extra_vep_params)

    # bucket inputs --------------------------------------------------------- #
    rsids:    list[tuple[str, str]] = []   # (id_for_api, original)
    genomic:  list[tuple[str, str]] = []
    vep_std:  list[tuple[str, str]] = []
    hgvs:     list[tuple[str, str]] = []   # protein-sub derived HGVS
    unresolved: list[str] = []

    protsub_matches: list[tuple[re.Match, str]] = []
    # original input -> normalized "ACC/refPOSalt" id (only for protein subs)
    uniprot_vid: dict[str, str] = {}
    for raw in identifiers:
        kind, val = _classify(raw)
        if kind == "rsid":
            rsids.append((val, raw))
        elif kind == "genomic":
            genomic.append((val, raw))
        elif kind == "protsub":
            protsub_matches.append((val, raw))
            uniprot_vid[raw] = (f'{val.group("acc")}/{val.group("ref").upper()}'
                                f'{val.group("pos")}{val.group("alt").upper()}')
        else:
            vep_std.append((val, raw))

    # resolve protein subs to ENSP HGVS ------------------------------------- #
    if protsub_matches:
        accs = [m.group("acc") for m, _ in protsub_matches]
        mapping = _map_uniprot_to_ensp(accs)
        for m, raw in protsub_matches:
            h = _protsub_to_hgvs(m, mapping.get(m.group("acc")))
            if h:
                hgvs.append((h, raw))
            else:
                unresolved.append(raw)

    missense: list[dict] = []
    matched_inputs: set[str] = set()

    def _run(endpoint, key, items, extra=None):
        """items: list of (api_value, original). POST in batches."""
        # keep original mapping; api values may collide so track per-batch order
        for chunk in _chunks(items, BATCH):
            payload = {key: [v for v, _ in chunk]}
            if extra:
                payload.update(extra)
            data = _post(endpoint, payload, params)
            if isinstance(data, dict):
                data = [data]
            # VEP echoes "input"; match it back to originals
            by_input: dict[str, str] = {}
            for v, orig in chunk:
                by_input.setdefault(v, orig)
            for rec in data:
                orig = by_input.get(rec.get("input"), rec.get("input"))
                rows = _extract_missense(rec, orig, uniprot_vid.get(orig))
                if rows:
                    missense.extend(rows)
                    matched_inputs.add(orig)

    if rsids:
        _run(f"/vep/{SPECIES}/id", "ids", rsids)
    if genomic:
        _run(f"/vep/{SPECIES}/region", "variants", genomic)
    if vep_std:
        # /region accepts the VEP default + many HGVS/SPDI strings; if some are
        # HGVS-only they'll error per-id, so route anything containing ':' that
        # looks like HGVS to /hgvs, the rest to /region.
        hgvs_like = [(v, o) for v, o in vep_std if re.search(r":[cgmnpr]\.", v)]
        region_like = [(v, o) for v, o in vep_std if (v, o) not in hgvs_like]
        if region_like:
            _run(f"/vep/{SPECIES}/region", "variants", region_like)
        if hgvs_like:
            hgvs.extend(hgvs_like)
    if hgvs:
        _run(f"/vep/{SPECIES}/hgvs", "hgvs_notations", hgvs)

    all_inputs = {raw for raw in identifiers}
    unmatched = sorted(all_inputs - matched_inputs - set(unresolved))

    return {
        "missense": missense,
        "unmatched": unmatched,
        "unresolved": unresolved,
    }


if __name__ == "__main__":
    import json
    test = [
        "rs699", "rs6265",
        "P00533 R132C", "P09874 S568F", "P00451 G41C",
        "P00533:R132C", "P09874:S568F", "P00451:G41C",
        "P00533/R132C", "P09874/S568F", "P00451/G41C",
        "chr14 89993420 A/G",
    ]
    result = fetch_variants(test)
    print(f"missense rows : {len(result['missense'])}")
    print(f"unmatched     : {result['unmatched']}")
    print(f"unresolved    : {result['unresolved']}")
    print(json.dumps(result["missense"][:3], indent=2))
