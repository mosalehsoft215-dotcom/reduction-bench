# Workstation for reducing HP/IP turbine test campaigns: Python with the IF97 steam tables.
FROM python:3.12-slim-bookworm
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps curl tmux \
 && rm -rf /var/lib/apt/lists/*
ENV OMP_NUM_THREADS=2 \
    OPENBLAS_NUM_THREADS=2 \
    MKL_NUM_THREADS=2 \
    PIP_NO_CACHE_DIR=1
RUN python3 -m pip install iapws==1.5.4 numpy==2.2.6 scipy==1.15.3
COPY data/ /app/data/
WORKDIR /app
RUN mkdir -p /app/solution /app/output
