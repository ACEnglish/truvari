FROM ubuntu:20.04

RUN apt-get -qq update \
  && DEBIAN_FRONTEND=noninteractive apt-get install -yq \
  bcftools \
  curl \
  python3-dev \
  python3-pip \
  samtools \
  tabix \
  vcftools \
  wget \
  && \
  rm -rf /var/lib/apt/lists/*

ADD . /opt/truvari-source
WORKDIR /opt/truvari-source

RUN wget https://mafft.cbrc.jp/alignment/software/mafft_7.505-1_amd64.deb \
    && dpkg -i mafft_7.505-1_amd64.deb && rm mafft_7.505-1_amd64.deb

RUN python3 -m pip install --upgrade pip && \
    python3 -m pip install setproctitle pylint anybadge coverage && \
    python3 -m pip install --upgrade setuptools && \
    python3 -m pip install ./

WORKDIR /data

ENTRYPOINT ["truvari"]
