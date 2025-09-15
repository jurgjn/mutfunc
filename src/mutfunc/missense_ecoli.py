
import importlib.resources, os.path, subprocess, sqlite3
from pathlib import Path
import pandas as pd

def uf(x):
    return '{:,}'.format(x)

path_data = path_gz = importlib.resources.files('mutfunc') / 'data/'
path_gz = importlib.resources.files('mutfunc') / 'data/missense_ecoli.sqlite.gz'
path = importlib.resources.files('mutfunc') / 'data/missense_ecoli.sqlite'
if not os.path.isfile(path):
    print('Downloading missense_ecoli.sqlite')
    path_data.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        args=(
            'curl', '--progress-bar', *('--output', path_gz),
            'https://polybox.ethz.ch/index.php/s/zj45yowi9D2ybQA/download/missense_ecoli.sqlite.gz',
            *('--stderr', '/dev/stdout'),
        ),
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    subprocess.run(args=('gunzip', path_gz))

with sqlite3.connect(path) as db:
    df_ = pd.read_sql_query(sql='SELECT COUNT() from missense', con=db)
    print(uf(df_.squeeze()), 'missense variants in', path)

def query_missense(variants):
    with sqlite3.connect(path) as db:
        df_ = pd.read_sql_query(sql='SELECT * FROM missense WHERE variant_id in ' + str(tuple(variants)), con=db)
    return df_
