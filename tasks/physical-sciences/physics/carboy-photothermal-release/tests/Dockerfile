# Image that runs the sealed carboys: the agent's NumPy/SciPy stack (the submitted program runs here) plus pytest.
FROM python:3.12-slim-bookworm
RUN pip install --no-cache-dir numpy==2.2.6 scipy==1.15.3 pytest==8.4.1 pytest-json-ctrf==0.3.5
RUN mkdir -p /app
RUN mkdir -p /app/solution /srv/carboys /var/lib/carboy-store /logs/verifier
COPY . /tests/
