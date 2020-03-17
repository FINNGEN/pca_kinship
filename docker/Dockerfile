FROM ubuntu:18.04

ENV TZ=Europe/Helsinki
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install mercurial --yes && \
apt-get install gcc zlib1g-dev wget make libcurl4-openssl-dev libbz2-dev liblzma-dev \
bzip2 curl unzip gawk python make gcc g++ libz-dev libbz2-dev liblzma-dev  littler libssl-dev  software-properties-common ghostscript parallel  --yes 


# R 3.6
RUN curl -O https://cloud.r-project.org/src/base/R-3/R-3.6.1.tar.gz && \
    tar xvzf R-3.6.1.tar.gz && \
    cd R-3.6.1 && \
    ./configure --with-x=no --with-blas="-lopenblas" && \
    make && mkdir -p /usr/local/lib/R/lib && make install && cd .. && rm -rf R-3.6.1*
# R PACKAGES
ADD ./docker/install_packages.R /usr/local/bin/
RUN chmod a+x /usr/local/bin/install_packages.R && install_packages.R

# PYTHON3 via conda
RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh -O ~/anaconda.sh && \
    /bin/bash ~/anaconda.sh -b -p /opt/conda && \
    rm ~/anaconda.sh
    
RUN  ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.profile && \
    find /opt/conda/ -follow -type f -name '*.a' -delete  && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
    echo "export PATH=\"/opt/conda/bin:$PATH\" " >> ~/.profile && \
    /opt/conda/bin/conda clean -afy

# PYTHON PACKAGES
ADD ./docker/requirements.txt .
RUN /opt/conda/bin/conda install --file requirements.txt --yes
ENV PATH "/opt/conda/bin:$PATH"

#PLINK 1.9
RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_latest.zip && unzip plink_* -d /usr/local/bin && rm plink*ip

#PLINK2
RUN wget http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_x86_64.zip && unzip plink2*.zip -d /usr/local/bin/ &&  rm plink2*.zip

#KING
RUN wget http://people.virginia.edu/~wc9c/KING/Linux-king.tar.gz && tar -xzvf Linux-king.tar.gz -C /usr/local/bin && rm Linux*king*.tar.gz

# BCFTOOLS
RUN wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 && tar -xvjf bcftools-1.10.2.tar.bz2 && cd bcftools-1.10.2 && \
./configure && make && make install && cd ..  && rm  bcftools-1.10.2.tar.bz2 && export BCFTOOLS_PLUGINS=~/bcftools-1.10.2/plugins/

ADD scripts ../scripts
ADD data ../data
