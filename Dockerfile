FROM ghcr.io/bimberlabinternal/discvr-base:latest

WORKDIR /pkg
COPY DESCRIPTION NAMESPACE /pkg/
COPY R/ /pkg/R/

RUN Rscript -e "BiocManager::install(ask = FALSE)" \
 && Rscript -e "devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'never')" \
 && R CMD INSTALL /pkg

COPY . /pkg/
RUN R CMD INSTALL /pkg

COPY exec/hellinger-enrichment.R /usr/local/bin/hellinger-enrichment
RUN chmod +x /usr/local/bin/hellinger-enrichment

ENTRYPOINT ["/usr/local/bin/hellinger-enrichment"]
CMD ["--help"]
