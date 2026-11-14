FROM python:3.12-slim-bookworm
RUN pip install --no-cache-dir numpy==2.2.6 scipy==1.15.3
# road-load workspace: coast-down notes plus the two reported campaigns
COPY data/ /app/data/
RUN mkdir -p /app/solution
WORKDIR /app/solution
