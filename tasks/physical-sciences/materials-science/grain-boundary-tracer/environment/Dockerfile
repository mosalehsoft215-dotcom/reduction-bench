FROM python:3.12-slim-bookworm
# Workspace for the sectioning task: NumPy and SciPy, two worked anneals with their figures, and the records guide.
ENV OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
RUN python3 -m pip install --no-cache-dir numpy==2.2.6 scipy==1.15.3
RUN mkdir -p /app/solution
COPY data /app/data
WORKDIR /app
