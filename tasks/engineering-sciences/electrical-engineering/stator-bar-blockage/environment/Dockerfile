FROM python:3.12-slim-bookworm
ENV PIP_NO_CACHE_DIR=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
RUN apt-get update && apt-get install -y --no-install-recommends procps && rm -rf /var/lib/apt/lists/*
RUN pip install numpy==2.2.6 scipy==1.15.3 pandas==2.2.3
COPY data /app/data
RUN mkdir -p /app/solution /app/output
WORKDIR /app
