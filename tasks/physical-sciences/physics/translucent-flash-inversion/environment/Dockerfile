FROM python:3.11-slim
RUN pip install --no-cache-dir numpy==2.1.3 scipy==1.14.1 pandas==2.2.3
RUN mkdir -p /app/data /app/output
COPY data/ /app/data/
WORKDIR /app
