FROM ubuntu

ENV DEBIAN_FRONTEND=noninteractive

RUN apt update -y
RUN apt-get update -y
RUN apt upgrade -y
RUN apt-get upgrade -y
RUN apt-get install build-essential -y
RUN apt-get install apt-utils apt-transport-https software-properties-common -y
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
RUN apt-get install r-base git curl -y
RUN apt-get install xml2 libxml2 libxml2-dev -y
RUN apt-get install libcurl4-openssl-dev software-properties-common -y
RUN apt-get install openjdk-8-jdk -y
RUN apt-get install libsodium-dev libssl-dev -y
RUN apt-get install libglu1-mesa-dev freeglut3-dev mesa-common-dev -y
RUN apt-get install bwidget -y
RUN apt update -y
RUN apt-get update -y

RUN cd /home/ && git clone https://github.com/angelovicilab/HAPPI_GWAS.git

WORKDIR /home/HAPPI_GWAS/

COPY docker_setup.R /
RUN Rscript /docker_setup.R

RUN rm -rf /home/HAPPI_GWAS/.git
RUN rm -rf /home/HAPPI_GWAS/docker_setup.R
RUN rm -rf /docker_setup.R

RUN Rscript setup.R
