# Analyst workstation for the activation capsules: Python 3.12, NumPy and SciPy pinned, the method notes and the two
# worked capsules under /app/data.
FROM python:3.12-slim-bookworm
RUN apt-get update && apt-get install -y --no-install-recommends procps tmux && apt-get clean \
 && rm -rf /var/lib/apt/lists/*
RUN python3 -m pip install --no-cache-dir "numpy==2.2.6" "scipy==1.15.3"
ENV OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
COPY data/activation_notes.md /app/data/activation_notes.md
COPY data/capsules/ /app/data/capsules/
RUN mkdir -p /app/solution
WORKDIR /app
