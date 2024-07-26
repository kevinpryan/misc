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
RUN mkdir -p /usr/local/bin/bwa-0.7.15
RUN apt-get update && apt-get install -y gcc zlib1g-dev
COPY bwa-install/bwa-0.7.15/* /usr/local/bin/bwa-0.7.15/
RUN cd /usr/local/bin/bwa-0.7.15 && make
ENV PATH /usr/local/bin/bwa-0.7.15:$PATH
RUN chmod +x /usr/local/bin/k8-linux /usr/local/bin/bwa-postalt.js


