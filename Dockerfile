################## BASE IMAGE ######################
FROM nfcore/base

################## METADATA ######################

LABEL base_image="nfcore/base"
LABEL version="1.1"
LABEL software="nf-hla-neo"
LABEL software.version="1.1"
LABEL about.summary="Container image containing all requirements for **hla-neo pipeline**"
LABEL about.home="http://github.com/IARCbioinfo/nf-hla-neo"
LABEL about.documentation="http://github.com/IARCbioinfo/nf-hla-neo/README.md"
LABEL about.license_file="http://github.com/IARCbioinfo/nf-hla-neo/LICENSE.txt"
LABEL about.license="GNU-3.0"

################## MAINTAINER ######################
MAINTAINER **digenovaa** <**digenovaa@fellows.iarc.fr**>

################## INSTALLATION ######################
COPY environment.yml /
RUN conda env update -n vep-nf -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/vep-nf/bin:$PATH
RUN conda env export --name vep-nf > vep-nf-v1.1.yml
RUN pip install vatools
#this shold be run locally
#RUN mkdir -p /vep-db/
#RUN vep_install -a cf -s homo_sapiens -y GRCh38 -c /vep-db/GRCh38/vep --CONVERT
