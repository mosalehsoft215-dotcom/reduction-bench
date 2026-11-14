FROM python:3.13-slim-bookworm
# two example radar stations, their figures and the record guide
ENV OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMBA_NUM_THREADS=2
RUN apt-get update && apt-get install -y --no-install-recommends tmux procps && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir numpy==2.2.6 scipy==1.15.3 pandas==2.2.3
WORKDIR /app
COPY data /app/data
RUN mkdir -p /app/solution
