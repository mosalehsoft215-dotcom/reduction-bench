FROM python:3.12-slim-bookworm
ENV PIP_RETRIES=10 PIP_DEFAULT_TIMEOUT=60
WORKDIR /app
COPY data /app/data
# numerical stack for the glide reduction; the guide and the two accepted days sit in /app/data
RUN python3 -m pip install --no-cache-dir "numpy==2.2.6" "scipy==1.15.3" && mkdir -p /app/solution
WORKDIR /app/solution
