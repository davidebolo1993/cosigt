FROM ubuntu:20.04
LABEL description="graph genotyper"
LABEL base_image="ubuntu:latest"
LABEL software="graph genotyper"
LABEL about.home="https://github.com/davidebolo1993/graph_genotyper"
LABEL about.license="GPLv3"

ARG DEBIAN_FRONTEND=noninteractive
#install basic libraries and python

WORKDIR /opt

RUN apt-get update
RUN apt-get -y install build-essential \
	wget git\
	bzip2 libbz2-dev \
	zlib1g zlib1g-dev \
	liblzma-dev \
	libssl-dev \
	libncurses5-dev \
	libz-dev \
	python3-distutils python3-dev python3-pip \ 
	libjemalloc-dev \
	cmake make g++ \
	libhts-dev \
	libzstd-dev \
	autoconf \
	libatomic-ops-dev \
	pkg-config \
	cargo \
    	&& apt-get -y clean all \
    	&& rm -rf /var/cache

##install odgi

RUN git clone --recursive https://github.com/pangenome/odgi.git
RUN cd odgi \
    && cmake -H. -DCMAKE_BUILD_TYPE=Generic -Bbuild \
    && cmake --build build -- -j $(nproc) \
    && rm -rf deps \
    && rm -rf .git \
    && rm -rf build \
    && apt-get -y clean all \
    && rm -rf /var/cache

ENV PATH /opt/odgi/bin:$PATH

##install vg

RUN mkdir -p vg \
	&& cd vg \
	&& wget https://github.com/vgteam/vg/releases/download/v1.43.0/vg \
    	&& chmod +x vg

ENV PATH /opt/vg:$PATH

##install samtools

RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 \
	&& tar -jxvf samtools-1.16.1.tar.bz2 \
	&& rm samtools-1.16.1.tar.bz2 \
	&& cd samtools-1.16.1 \
	&& ./configure \
	&& make \
	&& make install

#no need to env path, because of make install - should be sufficient
	
##install gafpack

RUN git clone https://github.com/ekg/gafpack.git \
	&& cd gafpack \
	&& cargo install --force --path .

ENV PATH /opt/gafpack/target/release:$PATH

##install additional python3 modules

RUN pip3 install pandas \
	numpy \
	scipy

##clone repo with commands required to run the entire genotyping + custom python script
RUN git clone https://github.com/davidebolo1993/graph_genotyper.git

ENV PATH /opt/graph_genotyper:$PATH
