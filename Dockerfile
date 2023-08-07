FROM ubuntu:latest
LABEL description="graph genotyper"
LABEL base_image="ubuntu:latest"
LABEL software="graph genotyper"
LABEL about.home="https://github.com/davidebolo1993/graph_genotyper"
LABEL about.license="GPLv3"

ARG DEBIAN_FRONTEND=noninteractive
#install basic libraries and python
#this is useful for having all the dependencies for which no official docker exist
#odgi is out there
#pggb same

WORKDIR /opt

RUN apt-get update

RUN apt-get -y install build-essential \
	software-properties-common \
	wget curl git \
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
	pigz \
	clang-14 \ 
	libomp5 libomp-dev libssl-dev libssl3 pkg-config \
	zip unzip

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

#and pgr-tk dependency
RUN cargo install --locked maturin

#install numpy
RUN pip3 install numpy

#pgr-tk installation
RUN wget https://github.com/GeneDx/pgr-tk/releases/download/v0.5.1/pgr-tk-v0.5.1.zip \
	&& unzip pgr-tk-v0.5.1.zip \
	&& rm pgr-tk-v0.5.1.zip \
	&& chmod +x release/* \
	&& chmod +x wheels/* \
	&& pip3 install wheels/pgrtk-0.5.1-cp310-cp310-linux_x86_64.whl

ENV PATH /opt/release:$PATH

#install wheels


#ln python to python3 -not used right now but, who knows?
RUN ln -s /usr/bin/python3 /usr/bin/python

##install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 \
	&& tar -jxvf samtools-1.18.tar.bz2 \
	&& rm samtools-1.18.tar.bz2 \
	&& cd samtools-1.18 \
	&& ./configure \
	&& make \
	&& make install

#no need to env path, because of make install - should be sufficient
##install bwa-mem

RUN git clone https://github.com/lh3/bwa.git \
	&& cd bwa \
	&& make

ENV PATH /opt/bwa:$PATH

##install bwa-mem2

RUN wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 \
	&& tar -jxvf bwa-mem2-2.2.1_x64-linux.tar.bz2 \
	&& rm bwa-mem2-2.2.1_x64-linux.tar.bz2

ENV PATH /opt/bwa-mem2-2.2.1_x64-linux:$PATH

##install minimap2

RUN wget https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2 \
	&& tar -jxvf minimap2-2.26_x64-linux.tar.bz2 \
	&& rm minimap2-2.26_x64-linux.tar.bz2

ENV PATH /opt/minimap2-2.26_x64-linux:$PATH

##install gafpack
##checkout to a specific version

RUN git clone https://github.com/ekg/gafpack.git \
	&& cd gafpack \
	&& cargo install --force --path .

ENV PATH /opt/gafpack/target/release:$PATH

##install gfainject

RUN git clone https://github.com/ekg/gfainject.git \
	&& cd gfainject \
	&& cargo install --force --path .

ENV PATH /opt/gfainject/target/release:$PATH

##install cosigt

RUN git clone https://github.com/davidebolo1993/cosigt.git \
	&& cd cosigt \
	&& go mod init cosigt \
	&& go mod tidy \
	&& go build cosigt

ENV PATH /opt/cosigt:$PATH