FROM python:3.12-slim-bookworm
# Coast-down analysis image: NumPy and SciPy pinned to the versions the program is promised.
RUN pip install --no-cache-dir numpy==2.2.6 scipy==1.15.3 && mkdir -p /app/solution /app/output
ENV OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
COPY data/ /app/data/
WORKDIR /app
