FROM python:3.12-slim-bookworm
ENV PYTHONDONTWRITEBYTECODE=1
RUN mkdir -p /app
# judge for the sealed scans: same numpy and scipy, pytest with CTRF reports
RUN pip install --no-cache-dir numpy==2.2.6 scipy==1.15.3 pytest==8.4.1 pytest-json-ctrf==0.3.5
RUN mkdir -p /app/solution /srv/zf-state /logs/verifier
COPY . /tests/
