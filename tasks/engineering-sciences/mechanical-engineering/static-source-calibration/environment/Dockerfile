FROM python:3.12-slim-bookworm
ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# tools for working in the container, then the numerical libraries the program may use
RUN apt-get update -qq && apt-get install -y -qq --no-install-recommends procps tmux curl \
 && rm -rf /var/lib/apt/lists/*
RUN python3 -m pip install --no-cache-dir numpy==2.2.6 scipy==1.15.3
COPY data/ /app/data/
RUN mkdir -p /app/solution /app/output
WORKDIR /app
