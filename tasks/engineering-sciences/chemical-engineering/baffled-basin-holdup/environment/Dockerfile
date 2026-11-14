FROM python:3.12-slim-bookworm
ENV PIP_RETRIES=10 PIP_DEFAULT_TIMEOUT=60
# workspace image: basin trial guide and two reference basin trials
RUN python3 -m pip install --no-cache-dir numpy==2.2.6 scipy==1.15.3
COPY data /app/data
RUN mkdir -p /app/solution
WORKDIR /app/solution
