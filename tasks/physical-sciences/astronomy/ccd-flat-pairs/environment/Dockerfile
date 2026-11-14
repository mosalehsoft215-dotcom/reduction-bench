FROM python:3.12-slim-bookworm
# Bench records of the two signed-off sensors, the bench guide, and the two numerical libraries.
COPY data/ /app/data/
RUN pip3 install --no-cache-dir numpy==2.2.6 scipy==1.15.3
RUN mkdir -p /app/solution
WORKDIR /app/data
