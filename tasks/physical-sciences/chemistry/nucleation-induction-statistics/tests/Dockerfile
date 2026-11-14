FROM python:3.12-slim-bookworm
# verifier image for the held-back crystalliser runs: the workspace's NumPy and SciPy, procps for the process sweep,
# pytest with its CTRF reporter
ENV PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
RUN apt-get update -qq && apt-get install -y -qq --no-install-recommends procps && rm -rf /var/lib/apt/lists/*
RUN mkdir -p /app
RUN python3 -m pip install --no-cache-dir numpy==2.2.6 scipy==1.15.3 pytest==8.4.1 pytest-json-ctrf==0.3.5
RUN mkdir -p /app/solution /srv/nr-work /logs/verifier
COPY . /tests/
