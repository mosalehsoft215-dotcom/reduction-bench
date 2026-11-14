FROM python:3.12-slim-bookworm
# agent workspace: two example plant tests and the record guide in /app/data; NumPy and SciPy pinned
ENV PIP_NO_CACHE_DIR=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMBA_NUM_THREADS=2
RUN apt-get update \
 && apt-get install -y --no-install-recommends curl tmux procps \
 && rm -rf /var/lib/apt/lists/*
RUN python -m pip install "numpy==2.2.6" "scipy==1.15.3"
WORKDIR /app
COPY data/ /app/data/
RUN mkdir -p /app/output /app/solution
