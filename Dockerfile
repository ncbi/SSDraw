FROM python:3.11-slim

ARG VERSION=0.1.6

WORKDIR /app

RUN apt-get update && apt-get install -y \
    dssp \
    make \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip
RUN pip install ssdraw==${VERSION}

ENTRYPOINT ["ssdraw"]