FROM python:3.11-slim-bookworm
# what the agent works in: two worked fountain campaigns, their figures and the campaign record guide
ENV MKL_NUM_THREADS=2
ENV OMP_NUM_THREADS=2
ENV OPENBLAS_NUM_THREADS=2
RUN apt-get update -qq \
 && DEBIAN_FRONTEND=noninteractive apt-get install -qq -y --no-install-recommends tmux procps \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*
RUN pip3 install --no-compile --no-cache-dir scipy==1.15.3 numpy==2.2.6
RUN mkdir -p /app/solution
COPY data/ /app/data/
WORKDIR /app
