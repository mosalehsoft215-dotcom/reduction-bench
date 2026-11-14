# Observatory reduction workstation: K-band calibration notes and two calibrated tracks under /app/data.
FROM python:3.12-slim-bookworm
RUN apt-get update -qq \
 && apt-get install -y -qq --no-install-recommends procps tmux \
 && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir numpy==2.2.6 scipy==1.15.3
ENV OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
RUN mkdir -p /app/solution /app/data
COPY data/kband_notes.md /app/data/
COPY data/tracks /app/data/tracks
WORKDIR /app
