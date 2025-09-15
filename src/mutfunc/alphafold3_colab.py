
def install():
    get_ipython().system('''
%%shell
# Adapted from the official dockerfile: https://github.com/google-deepmind/alphafold3/blob/main/docker/Dockerfile
rm -rf /content/sample_data

# /app/alphafold - alphafold3 git repo checkout
mkdir /app
cd /app
git clone https://github.com/google-deepmind/alphafold3.git alphafold
apt-get install --yes --quiet python3.12 python3-pip python3.12-venv python3.12-dev git wget gcc g++ make zlib1g-dev zstd
''')