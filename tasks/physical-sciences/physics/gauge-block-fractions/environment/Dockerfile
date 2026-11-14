# Workstation image for the gauge-block runs. Debian-based CPython 3.12 with the two numeric libraries pinned; the
# lab's interferometry notes and the worked blocks are placed read-only under /app/data.
FROM python:3.12-slim-bookworm
RUN apt-get update \
 && apt-get install -y --no-install-recommends tmux procps \
 && rm -rf /var/lib/apt/lists/* \
 && python3 -m pip install --no-cache-dir "numpy==2.2.6" "scipy==1.15.3"
ENV MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
COPY data/ /app/data/
RUN mkdir -p /app/solution && chmod -R a-w /app/data
WORKDIR /app
