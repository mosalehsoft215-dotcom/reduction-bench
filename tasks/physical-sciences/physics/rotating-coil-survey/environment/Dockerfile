FROM python:3.13-slim-bookworm
# the agent's workspace: two worked magnet surveys, their figures and the bench record guide
ENV OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
RUN apt-get update -q && apt-get install -y -q --no-install-recommends procps tmux && rm -rf /var/lib/apt/lists/*
RUN python -m pip install --no-cache-dir numpy==2.2.6 scipy==1.15.3
COPY data /app/data
RUN mkdir -p /app/solution
WORKDIR /app
