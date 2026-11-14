# Counting-room workstation: the bench notes and the two released assays go to /app/data, with Python, NumPy and SciPy.
FROM python:3.12-slim-bookworm
RUN apt-get update && apt-get install -y --no-install-recommends procps tmux && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir numpy==2.2.6 scipy==1.15.3
ENV OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
WORKDIR /app
COPY data/gamma_bench_notes.md /app/data/gamma_bench_notes.md
COPY data/assays/ /app/data/assays/
RUN mkdir -p /app/solution
