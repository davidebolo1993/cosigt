FROM ubuntu:latest
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
	software-properties-common \
	wget curl git\
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
	pigz 

#install golang
RUN add-apt-repository ppa:longsleep/golang-backports

RUN apt-get -y install golang-go \
	&& apt-get -y clean all \
	&& rm -rf /var/cache

#install rust
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y

ENV PATH="/root/.cargo/bin:${PATH}"

#and update
RUN rustup update

#ln python to python3 -not used right now but, who knows?
RUN ln -s /usr/bin/python3 /usr/bin/python

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

##install bwa-mem

RUN git clone https://github.com/lh3/bwa.git \
	&& cd bwa \
	&& make

ENV PATH /opt/bwa:$PATH
	
##install gafpack
##checkout to a specific version

RUN git clone https://github.com/ekg/gafpack.git \
	&& cd gafpack \
	&& git checkout ad31875b6914d964c6fd72d1bf334f0843538fb6 \
	&& cargo install --force --path .

ENV PATH /opt/gafpack/target/release:$PATH

##install gfainject

RUN git clone https://github.com/ekg/gfainject.git \
	&& cd gfainject \
	&& cargo install --force --path .

ENV PATH /opt/gfainject/target/release:$PATH

##install cosigt

RUN git clone https://github.com/davidebolo1993/graph_genotyper.git \
	&& cd graph_genotyper \
	&& go mod init cosigt \
	&& go mod tidy \
	&& go build cosigt

ENV PATH /opt/graph_genotyper:$PATH
