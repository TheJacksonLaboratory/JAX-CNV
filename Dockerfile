FROM ubuntu:16.04

MAINTAINER Wan-Ping Lee <wanping.lee@gmail.com>

# Packaged dependencies
RUN apt-get update && apt-get install -y \
	git \
	make \
	gcc \
	build-essential g++ \
	libcurl4-openssl-dev \
	libbz2-dev \
	liblzma-dev \
	libz-dev \
	libssl1.0.0 \
	libssl-dev \
	automake \
	autoconf \
	wget \
	r-base

# Make a folder for tools
RUN cd / && mkdir -p tools && cd /tools

# Git clone JAX-CNV
RUN git clone --recursive https://github.com/TheJacksonLaboratory/JAX-CNV.git

# Build JAX-CNV
RUN cd JAX-CNV \
	&& make

# Define default command.
CMD ["/tools/JAX-CNV/bin/JAX-CNV"]
