FROM python:3.12-slim-bookworm
LABEL description="Workspace with two example irradiance missions and their record guide"
ENV OMP_NUM_THREADS=2 \
    OPENBLAS_NUM_THREADS=2 \
    MKL_NUM_THREADS=2 \
    NUMBA_NUM_THREADS=2
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps tmux \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir numpy==2.2.6 scipy==1.15.3 pandas==2.2.3
COPY data /app/data
WORKDIR /app
RUN mkdir -p /app/solution
