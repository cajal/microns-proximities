FROM at-docker.ad.bcm.edu:5000/microns-base 

LABEL mantainer="Stelios Papadopoulos <spapadop@bcm.edu>"

# Install additional packages

# Update and rebuild jupyterlab 
RUN python -m pip install --no-cache-dir --upgrade jupyterlab

WORKDIR /root
ARG CLOUDVOLUME_TOKEN
RUN mkdir -p .cloudvolume/secrets
RUN echo "{\"token\": \"${CLOUDVOLUME_TOKEN:-}\"}" > .cloudvolume/secrets/cave-secret.json

# CURRENT PACKAGE
COPY . /src/microns-proximities
RUN python -m pip install --no-cache-dir /src/microns-proximities/python/microns-proximities
RUN python -m pip install --no-cache-dir /src/microns-proximities/python/microns-proximities-api
