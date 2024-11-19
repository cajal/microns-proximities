FROM at-docker:5000/microns-base:cuda11.8.0-python3.8

LABEL mantainer="Stelios Papadopoulos <spapadop@bcm.edu>"

# Install additional packages

# Update and rebuild jupyterlab 
RUN python -m pip install --no-cache-dir --upgrade jupyterlab

WORKDIR /
ARG CLOUDVOLUME_TOKEN
RUN mkdir -p .cloudvolume/secrets
RUN echo "{\"token\": \"${CLOUDVOLUME_TOKEN:-}\"}" > .cloudvolume/secrets/cave-secret.json

# CURRENT PACKAGE
COPY . /src/microns-proximities
RUN python -m pip install -e /src/microns-proximities/python/microns-proximities
RUN python -m pip install -e /src/microns-proximities/python/microns-proximities-api
