FROM nfcore/base:2.1
LABEL authors="Kevin Ryan" \
      description="Docker image containing all requirements for bwa mem alt-aware alignment and subsetting to HLA regions pipeline"
WORKDIR ./
COPY bwa.yml ./
RUN conda env create -f bwa.yml && conda clean -a
ENV PATH /opt/conda/envs/bwa/bin:$PATH
# copy the necessary scripts to the container
COPY k8/k8-linux /usr/local/bin/
COPY bwakit-postalt/bwa-postalt.js /usr/local/bin/
RUN chmod +x /usr/local/bin/k8-linux /usr/local/bin/bwa-postalt.js


