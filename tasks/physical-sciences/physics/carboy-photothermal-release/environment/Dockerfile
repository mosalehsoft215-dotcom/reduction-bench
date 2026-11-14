# Bench PC of the photothermal carboy bench: Python 3.12 with NumPy and SciPy, the bench notes and two worked carboys.
FROM python:3.12-slim-bookworm
RUN apt-get update && apt-get install -y --no-install-recommends procps tmux && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir numpy==2.2.6 scipy==1.15.3
ENV OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
WORKDIR /app
COPY data/carboy_lens_notes.md /app/data/carboy_lens_notes.md
COPY data/carboys/ /app/data/carboys/
RUN mkdir -p /app/solution
