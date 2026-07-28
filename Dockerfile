FROM rocker/r-ver:4.4.2

# System libraries for Seurat/anndata stack before R package installs.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libhdf5-dev \
        pandoc \
        python3 \
        python3-pip \
        python3-venv \
    && rm -rf /var/lib/apt/lists/*

# Python anndata for h5ad ingestion via reticulate.
RUN python3 -m venv /opt/venv \
    && /opt/venv/bin/pip install --no-cache-dir anndata

ENV RETICULATE_PYTHON=/opt/venv/bin/python
ENV PATH="/opt/venv/bin:${PATH}"

WORKDIR /pkg

COPY DESCRIPTION NAMESPACE R/ /pkg/
RUN R -e "install.packages(c('remotes', 'testthat'), repos='https://cloud.r-project.org')" \
    && R -e "remotes::install_deps(dependencies=TRUE, repos='https://cloud.r-project.org')" \
    && R CMD INSTALL /pkg

COPY . /pkg/
RUN R CMD INSTALL /pkg

COPY exec/hellinger-enrichment.R /usr/local/bin/hellinger-enrichment
RUN chmod +x /usr/local/bin/hellinger-enrichment

ENTRYPOINT ["/usr/local/bin/hellinger-enrichment"]
CMD ["--help"]
