# Based on: https://docs.streamlit.io/deploy/tutorials/docker
# Updated to python 3.13 & removed software-properties-common as deprecated

FROM python:3.13-slim

# python3-pymol via apt-get as pip package missing for linux/arm64
# libboost-all-dev required by CombFold
RUN DEBIAN_FRONTEND=noninteractive \
apt-get update --quiet \
&& apt-get install --yes --quiet build-essential curl git python3-pymol libboost-all-dev \
&& rm -rf /var/lib/apt/lists/*

COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

# Suppress `warning: Failed to hardlink files; falling back to full copy. This may lead to degraded performance.`
ENV UV_LINK_MODE=copy

# https://docs.astral.sh/uv/guides/integration/docker/#installing-a-package
ENV UV_SYSTEM_PYTHON=1

# https://docs.astral.sh/uv/guides/integration/docker/#compiling-bytecode
ENV UV_COMPILE_BYTECODE=1

WORKDIR /app/mutfunc
COPY data/* data/

COPY pyproject.toml .

RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync

COPY *.py .
COPY assets/* assets/
COPY pages/* pages/
COPY util/* util/

EXPOSE 8050

ENTRYPOINT ["uv", "run", "gunicorn", "app:server", "--bind", "0.0.0.0:8050"]
