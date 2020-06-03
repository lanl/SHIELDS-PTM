FROM ubuntu:18.04

MAINTAINER steven.k.morley@gmail.com

ENV GIT_SSL_NO_VERIFY=1
RUN apt-get update -qq
RUN apt-get install -y ssh make git gfortran gcc g++
RUN apt-get install -y python3-pip
RUN pip3 install numpy scipy matplotlib
RUN pip3 install spacepy

ENV http_proxy=http://proxyout.lanl.gov:8080
ENV https_proxy=http://proxyout.lanl.gov:8080
ENV no_proxy=localhost,127.0.0.1,lanl.gov,.lanl.gov,*.lanl.gov
