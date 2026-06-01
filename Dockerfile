# PeTEM pipeline container
FROM rocker/r-ver:4.3.2

LABEL maintainer="PeTEM maintainers" \
      description="Container image for the Promoter-embedded TE Methylation (PeTEM) pipeline"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    bash python3 python3-pip python3-venv \
    bedtools samtools \
    perl gzip gawk \
    libxml2-dev libcurl4-openssl-dev libssl-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    make g++ \
 && rm -rf /var/lib/apt/lists/*

# Python packages
RUN python3 -m pip install --no-cache-dir pandas

# R packages required by the pipeline
RUN install2.r --error \
    optparse dplyr tidyr zoo reshape2 stringr ggplot2 gplots ggalluvial ggpointdensity RColorBrewer viridis rlang ggbreak

# Copy the PeTEM tools into the image
COPY . /opt/petem

# Ensure helper binaries are executable
RUN chmod +x /opt/petem/petem \
    /opt/petem/run_PeTEM.sh \
    /opt/petem/env_check.sh \
    /opt/petem/setup.sh \
    /opt/petem/0_preprocessing.sh \
    /opt/petem/1_TE_distribution.sh \
    /opt/petem/3_TE_impact_distance.sh \
    /opt/petem/3_1_TE_impact_distance_preprocess.sh \
    /opt/petem/3_2_TE_impact_distance_plot.sh \
    /opt/petem/wigToBigWig \
    /opt/petem/bigWigAverageOverBed

ENV PETEM_HOME=/opt/petem
WORKDIR /data
VOLUME ["/data"]

ENTRYPOINT ["bash", "-lc", "/opt/petem/env_check.sh >/tmp/petem_env_check.log && exec /opt/petem/petem \"$@\"", "--"]
CMD ["--help"]
