FROM python:3.12-slim-bookworm
# rail jacking workspace: notes, profile table and two settled runs under /app/data
WORKDIR /app
COPY data /app/data
RUN mkdir -p /app/solution && pip install --no-cache-dir scipy==1.15.3 numpy==2.2.6
ENV OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
