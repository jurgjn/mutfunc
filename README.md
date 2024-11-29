Download code & dependencies:
```
git clone git@github.com:evocellnet/mutfunc.git
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Download table of missense variants (~25G):
```
scp euler.ethz.ch:/cluster/work/beltrao/jjaenes/23.11.01_human_protein_map/missense_23.11.1.sqlite resources/
```

Run:
```
streamlit run Mutfunc_Home.py
```
