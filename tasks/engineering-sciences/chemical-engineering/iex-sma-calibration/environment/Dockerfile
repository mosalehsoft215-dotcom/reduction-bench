FROM python:3.12-slim-bookworm
RUN apt-get update -qq \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq --no-install-recommends tmux curl procps \
 && rm -rf /var/lib/apt/lists/*
RUN python3 -m pip install --no-cache-dir scipy==1.15.3 numpy==2.2.6
ENV OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
WORKDIR /app
COPY data/ /app/data/
RUN mkdir -p /app/output /app/solution
