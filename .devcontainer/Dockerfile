FROM node:20

RUN apt-get update && apt-get install -y --no-install-recommends \
  less \
  git \
  procps \
  sudo \
  fzf \
  zsh \
  man-db \
  unzip \
  gnupg2 \
  gh \
  iptables \
  ipset \
  iproute2 \
  dnsutils \
  aggregate \
  jq \
  nano \
  vim \
  && apt-get clean && rm -rf /var/lib/apt/lists/*

# Set up uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/
ENV UV_LINK_MODE=copy
ENV UV_PROJECT_ENVIRONMENT=/home/node/venv

# Download foldcomp platform-specific binary
ARG TARGETPLATFORM
RUN if [ "$TARGETPLATFORM" = "linux/amd64" ]; then \
        curl https://mmseqs.com/foldcomp/foldcomp-linux-x86_64.tar.gz | tar xvfz - -C /usr/bin/; \
    elif [ "$TARGETPLATFORM" = "linux/arm64" ]; then \
        curl https://mmseqs.com/foldcomp/foldcomp-linux-arm64.tar.gz | tar xvfz - -C /usr/bin/; \
    else \
        echo Unsupported platform; \
        exit 1; \
    fi

# Set `DEVCONTAINER` environment variable to help with orientation
ENV DEVCONTAINER=true

# Create workspace and config directories and set permissions
RUN mkdir -p /workspace /home/node/.claude /home/node/venv && \
  chown -R node:node /workspace /home/node/.claude /home/node/venv

WORKDIR /workspace

# Install global packages
ENV NPM_CONFIG_PREFIX=/usr/local/share/npm-global
ENV PATH=$PATH:/usr/local/share/npm-global/bin

# Set the default shell to zsh rather than sh
ENV SHELL=/bin/zsh

# Set the default editor and visual
ENV EDITOR=nano
ENV VISUAL=nano

# Install Claude
RUN npm install -g @anthropic-ai/claude-code

# Set up non-root user
USER node
