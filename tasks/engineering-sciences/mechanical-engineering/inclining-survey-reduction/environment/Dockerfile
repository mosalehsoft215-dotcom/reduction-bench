FROM python:3.12-slim-bookworm
ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# shell tools for the session, then the two numerical libraries the program can rely on
RUN apt-get update -qq \
 && apt-get install -y -qq --no-install-recommends procps less curl \
 && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir numpy==2.2.6 scipy==1.15.3
COPY data/ /app/data/
RUN mkdir -p /app/solution /app/output
WORKDIR /app
