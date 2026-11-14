FROM python:3.12-slim-bookworm
# test-stand workspace: two example firings and the notes on their records in /app/data
RUN apt-get update && apt-get install -y --no-install-recommends procps curl tmux && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir --disable-pip-version-check numpy==2.2.6 scipy==1.15.3
ENV OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMBA_NUM_THREADS=2
WORKDIR /app
COPY data/ /app/data/
RUN mkdir -p /app/output /app/solution
