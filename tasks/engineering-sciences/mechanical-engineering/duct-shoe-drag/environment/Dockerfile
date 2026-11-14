FROM python:3.12-slim-bookworm
ENV OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONDONTWRITEBYTECODE=1
RUN pip install --no-cache-dir numpy==2.2.6 scipy==1.15.3
WORKDIR /app
COPY data/ /app/data/
RUN mkdir -p /app/output /app/solution
