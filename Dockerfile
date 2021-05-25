################## BASE IMAGE ######################
FROM nfcore/base

################## METADATA ######################

LABEL base_image="nfcore/base"
LABEL version="1.0"
LABEL software="nf-hla-neo"
LABEL software.version="1.0"
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
RUN conda env export --name vep-nf > vep-nf-v1.0.yml
RUN mkdir -p /vep-db/
RUN vep_install -a cf -s homo_sapiens -y GRCh38 -c /vep-db/GRCh38/vep --CONVERT
#we get the code
#PURPLE
#ADD https://github.com/hartwigmedical/hmftools/releases/download/purple-v2.52/purple-2.52.jar /hmftools/jars
#COVALT
#ADD https://github.com/hartwigmedical/hmftools/releases/download/cobalt-v1.11/cobalt-1.11.jar /hmftools/jars
#RUN cd /home/programs && unzip master.zip
# amber
#ADD https://github.com/hartwigmedical/hmftools/releases/download/amber-v3.5/amber-3.5.jar /hmftools/jars
# databases
#ADD
