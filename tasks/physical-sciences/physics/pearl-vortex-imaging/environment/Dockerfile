FROM python:3.12-slim-bookworm
# Public SQUID maps and the numerical libraries.
COPY data/ /app/data/
ENV OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
RUN pip3 install --no-cache-dir numpy==2.2.6 scipy==1.15.3
RUN mkdir -p /app/solution
WORKDIR /app/data
