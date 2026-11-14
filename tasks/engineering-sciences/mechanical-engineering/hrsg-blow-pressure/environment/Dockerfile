# Commissioning office PC: HRSG trial-blow records, Python 3.12 with the IAPWS-IF97 steam tables, NumPy and SciPy.
FROM python:3.12-slim-bookworm
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qq update \
 && apt-get -qq install --no-install-recommends tmux procps \
 && apt-get clean && rm -rf /var/lib/apt/lists/*
RUN python3 -m pip install --no-cache-dir --disable-pip-version-check iapws==1.5.4 scipy==1.15.3 numpy==2.2.6
ENV PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2 MKL_NUM_THREADS=2
WORKDIR /app
RUN mkdir -p /app/output
COPY data/ /app/data/
