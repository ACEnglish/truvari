FROM ubuntu:20.04

RUN apt-get -qq update && apt-get install -yq \
  python3-pip \
  python3-dev \
  curl \
  && \
  rm -rf /var/lib/apt/lists/*

ADD . /opt/truvari-source
WORKDIR /opt/truvari-source

RUN python3 -m pip install --upgrade pip && \
    python3 -m pip install setproctitle pylint anybadge coverage && \
    python3 -m pip install --upgrade setuptools && \
    python3 -m pip install ./

WORKDIR /data

ENTRYPOINT ["truvari"]
