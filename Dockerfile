# This dockerfile used for 4Cseqpipe.
# VERSION 1 - EDITION 1
# Powered by ChangeGene LLC.
# A Bioinformatics Solution Provider from Harvard.
# Created by docker@changegene.com at 20170731.
# Offical web: https://github.com/changegene/4Cseqpipe
# The follow is modify log:

# ubuntu version choice the version earlier than your Software or Pipeline.  X.04 is the LTS version.
# ubuntu version list: https://wiki.ubuntu.com/Releases
ARG  Image_VERSION=12.04
FROM ubuntu:$Image_VERSION


MAINTAINER Skywalker <skywalker@changegene.com>

ENV DIRPATH /Project
# ENV R_big_version=3
# ENV R_version=3.0.2
ENV R_big_version=2
ENV R_version=2.13.1
ENV python2_version=2.7
ENV python3_version=3.2
ENV java_version=7
# ENV image_name=
# ENV version=
# ENV docker_name=$image_name
#ENV PG_VERSION 9.3.4
#RUN curl -SL http://example.com/postgres-$PG_VERSION.tar.xz | tar -xJC /usr/src/postgress && …
#ENV PATH=/usr/local/postgres-$PG_MAJOR/bin:$PATH
WORKDIR $DIRPATH


###########################Install base tools and software ####################
##### Install basic develop tools ###############
RUN apt-get -y update && apt-get install -y \
 apt-utils \
 build-essential \
 gcc-multilib \
 # libXt-dev \
 iputils-ping \
 zlib1g-dev

###### Install basic commond ######
RUN apt-get install -y \
 less \
 vim \
 unzip \
 wget \
 zip

###### Install python、perl、git######
#RUN apt-get install -y python3-dev python3-pip git perl
###### Install JAVA ##############
#RUN apt-get install -y software-properties-common python-software-properties
#RUN add-apt-repository ppa:webupd8team/java && apt-get update
#RUN echo "oracle-java7-installer shared/accepted-oracle-license-v1-1 select true" | debconf-set-selections
#RUN echo "oracle-java7-installer shared/accepted-oracle-license-v1-1 seen true" | debconf-set-selections
#RUN apt-get install -y oracle-java7-installer  &&  apt-get install -y oracle-java7-set-default
###### Install R #####
#RUN apt-get install -y r-base=$R_version r-base-dev=$R_version

#-------------mkdir ------------
RUN mkdir $DIRPATH/Software/ $DIRPATH/Data/ $DIRPATH/Database/ $DIRPATH/Pipeline/ $DIRPATH/Output/ $DIRPATH/Result/

#-------------Install R------------
# RUN cd $DIRPATH/Software/ && wget https://mirrors.tuna.tsinghua.edu.cn/CRAN/src/base/R-$R_big_version/R-$R_version.tar.gz && tar -zxvf R-$R_version.tar.gz && rm -f R-$R_version.tar.gz && cd R-$R_version/ && ./configure --prefix=$DIRPATH/Software/R-$R_version/ && make && make install
# wget https://mirrors.tuna.tsinghua.edu.cn/CRAN/src/base/R-2/R-2.13.1.tar.gz

RUN apt-get install -y \
 aptitude \
 # brltty-x11 \
 gfortran \
 # gnuplot-x11 \
 libcairo2-dev \
 # libcurl4-openssl-dev \
 libg2-dev \
 # libghc-x11-dev \
 # libpng \
 # libpng-devel \
 # libtiff \
 # libtiff-devel \
 # libjpeg-turbo \
 # libjpeg-turbo-devel \
 libpango1.0-dev \
 libreadline-dev \
 libsbuild-dev \
 # libx11-dev \
 # libx11-6 \
 libxt-dev \
 r-base-dev
 # x11-utils \
 # x11proto-core-dev \
 xauth \
 xfonts-base \
# xorg-dev \
# xserver-xorg-core \
 xvfb

COPY ./Software/R-$R_version.tar.gz $DIRPATH/Software/
RUN cd $DIRPATH/Software/ &&\
 tar -zxvf R-$R_version.tar.gz && \
 rm -f R-$R_version.tar.gz && cd R-$R_version/ && \
 ./configure --enable-R-shlib --with-libpng --with-jpeglib --with-libtiff --with-x --prefix=$DIRPATH/Software/R-$R_version/ && \
 make && make install
# RUN cd $DIRPATH/Software/ && wget https://mirrors.tuna.tsinghua.edu.cn/CRAN/src/base/R-$R_big_version/R-$R_version.tar.gz && tar -zxvf R-$R_version.tar.gz && rm -f R-$R_version.tar.gz && cd R-$R_version/ && ./configure --prefix=$DIRPATH/Software/R-$R_version/ && make && make install
ENV PATH $DIRPATH/Software/R-$R_version/bin:$PATH

#-------------Install R packages------------
# RUN cd $DIRPATH/Software/R/R-packages && wget https://mirrors.tuna.tsinghua.edu.cn/CRAN/src/contrib/zoo_1.8-0.tar.gz && tar -zxvf zoo_1.8-0.tar.gz && rm -f zoo_1.8-0.tar.gz
# RUN cd $DIRPATH/Software/ && R CMD INSTALL mypkg -l /Software/R/R-packages/

# wget download R packages commond backup.
# WORKDIR $DIRPATH/Software/R/
# wget https://mirrors.tuna.tsinghua.edu.cn/CRAN/src/contrib/Archive/Cairo/Cairo_1.5-1.tar.gz
# wget https://mirrors.tuna.tsinghua.edu.cn/CRAN/src/contrib/zoo_1.8-0.tar.gz

COPY ./Software/R/ $DIRPATH/Software/R/
# RUN R CMD INSTALL $DIRPATH/Software/R/R-packages/* -l $DIRPATH/Software/R/R-packages/
RUN R CMD INSTALL $DIRPATH/Software/R/R-packages/*


######install convert###
#RUN apt-get install -y imagemagick

# ADD ./requirements.txt /requirements.txt
# RUN pip install -r /requirements.txt
# ADD ./ /usr/local/git/my_app
# DIRPATH /usr/local/git/my_app
# CMD python ./main.py


########################拷贝流程脚本(v0.3)#################

# ADD ./Pipeline/ $DIRPATH/Pipeline
# RUN cd $DIRPATH/Software && tar zxvf 4cseq_pipe.tgz
# WORKDIR $DIRPATH
# RUN echo $DIRPATH


# Useage:
#docker  build --rm -t image_name:version

#docker run --rm -it --name docker_name  -v host_path:docker_path image_name:version
#docker run --rm -it --name 4cseqpipe 4cseqpipe:1.0 /bin/bash
#docker run --rm -it --name 4cseqpipe -v ~/Docker/4Cseqpipe:/project/ 4cseqpipe:1.0 /bin/bash

#docker run --rm -it --name docker_name  -v ./:/ -v ~/input:/Data -v ~/project:/Software image_name:version
