FROM ubuntu:18.04

ENV TZ=Europe/Helsinki
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install mercurial --yes && \
apt-get install gcc zlib1g-dev wget make libcurl4-openssl-dev libbz2-dev liblzma-dev \
bzip2 curl unzip gawk python python3 make gcc g++ libz-dev libbz2-dev liblzma-dev  littler libssl-dev python3-numpy python3-pandas python3-pip software-properties-common ghostscript --yes &&\
pip3 install networkx seaborn && pip3 install numpy --upgrade


RUN apt-get remove r-base-core --yes
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' && apt update && apt upgrade --yes
RUN apt install r-base --yes

ADD ./docker/install_packages.R /usr/local/bin/
RUN chmod a+x /usr/local/bin/install_packages.R && install_packages.R

ADD ./docker/requirements.txt .
RUN pip3 install -r requirements.txt && rm requirements.txt

#PLINK 1.9
RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_latest.zip && unzip plink_* -d /usr/local/bin && rm plink*zip

#PLINK2
RUN wget http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_x86_64.zip && unzip plink2*.zip -d /usr/local/bin/ &&  rm plink2*.zip

#KING
RUN wget http://people.virginia.edu/~wc9c/KING/Linux-king.tar.gz && tar -xzvf Linux-king.tar.gz -C /usr/local/bin && rm Linux*.tar.gz


ADD scripts ../scripts
ADD data ../data
