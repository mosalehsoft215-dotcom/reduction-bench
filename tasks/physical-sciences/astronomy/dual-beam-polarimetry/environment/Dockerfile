FROM python:3.12-slim-bookworm
RUN python3 -m pip install --no-cache-dir "numpy==2.2.6" "scipy==1.15.3" \
 && mkdir -p /app/solution
# polarimeter guide and the two visits that were already reduced
COPY data /app/data
WORKDIR /app
