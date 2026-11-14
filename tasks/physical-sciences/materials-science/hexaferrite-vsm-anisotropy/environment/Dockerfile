FROM python:3.12-slim-bookworm
ENV PIP_NO_CACHE_DIR=1 OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
RUN set -eux; apt-get update; \
    apt-get install --yes --no-install-recommends tmux curl procps; \
    apt-get clean; rm -rf /var/lib/apt/lists/*
RUN python -m pip install "scipy==1.15.3" "numpy==2.2.6"
WORKDIR /app
# two public samples and the description of the records
COPY data /app/data
RUN mkdir -p /app/solution
