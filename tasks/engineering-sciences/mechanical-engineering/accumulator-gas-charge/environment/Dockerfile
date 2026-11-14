FROM python:3.12-slim-bookworm
# workspace: accumulator gas-side notes and two evaluated units
COPY data /app/data
RUN python3 -m pip install --no-cache-dir scipy==1.15.3 numpy==2.2.6 && install -d /app/solution
WORKDIR /app/solution
