FROM ubuntu:22.04

RUN apt-get -qq update \
  && DEBIAN_FRONTEND=noninteractive apt-get install -yq \
  curl \
  python3-dev \
  python3-pip \
  wget \
  && \
  rm -rf /var/lib/apt/lists/*

RUN mkdir -p /opt/truvari-source/truvari/
COPY setup.py README.md pyproject.toml /opt/truvari-source/
COPY truvari/ /opt/truvari-source/truvari/
WORKDIR /opt/truvari-source

RUN wget https://mafft.cbrc.jp/alignment/software/mafft_7.505-1_amd64.deb \
    && dpkg -i mafft_7.505-1_amd64.deb && rm mafft_7.505-1_amd64.deb

RUN python3 -m pip install --upgrade pip && \
    python3 -m pip install setproctitle pylint anybadge coverage && \
    python3 -m pip install --upgrade setuptools && \
    python3 -m pip install ./

WORKDIR /data

ENTRYPOINT ["truvari"]
