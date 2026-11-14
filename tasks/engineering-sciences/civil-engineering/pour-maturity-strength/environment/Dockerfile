FROM python:3.12-slim-bookworm
# site and lab records for the maturity reduction
RUN pip install --no-cache-dir numpy==2.2.6 scipy==1.15.3
COPY data /app/data
RUN mkdir -p /app/solution
WORKDIR /app
