# Image for the governor test records: the analysis libraries the program may use.
FROM python:3.12-slim-bookworm
RUN mkdir -p /app/solution /app/output
RUN pip install --no-cache-dir scipy==1.15.3 numpy==2.2.6
ENV OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
COPY data/ /app/data/
WORKDIR /app
