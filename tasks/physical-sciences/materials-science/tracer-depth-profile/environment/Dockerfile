FROM python:3.12-slim-bookworm
# Coupon records and the notes that explain them, with the numerical libraries an evaluation may use.
ENV OPENBLAS_NUM_THREADS=2
ENV OMP_NUM_THREADS=2
ENV MKL_NUM_THREADS=2
COPY data/ /app/data/
RUN pip3 install --no-cache-dir numpy==2.2.6 scipy==1.15.3 && mkdir -p /app/solution
WORKDIR /app
