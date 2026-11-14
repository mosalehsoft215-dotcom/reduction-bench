FROM python:3.12-slim-bookworm
# roof tie workspace: the notes and two worked backstays
RUN pip3 install --no-cache-dir numpy==2.2.6 scipy==1.15.3
COPY data/ /app/data/
RUN mkdir -p /app/solution
WORKDIR /app
