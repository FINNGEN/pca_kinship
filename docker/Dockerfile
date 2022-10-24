FROM eu.gcr.io/finngen-refinery-dev/bioinformatics:0.6

ENV TZ=Europe/Helsinki
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

ARG R_VER=4.1.2

RUN apt-get update && apt-get upgrade --yes

# R PACKAGES
# RUN curl -O https://cloud.r-project.org/src/base/R-4/R-${R_VER}.tar.gz && \
#     tar xvzf R-${R_VER}.tar.gz && \
#     cd R-${R_VER} && \
#     ./configure --with-x=no --enable-R-shlib=yes --with-blas="-lopenblas" && \
#     make && mkdir -p /usr/local/lib/R/lib && make install && cd .. && rm -rf R-${R_VER}*

ADD ./docker/install_packages.R /usr/local/bin/
RUN chmod a+x /usr/local/bin/install_packages.R && install_packages.R


RUN apt-get install libproj-dev proj-data proj-bin --yes &&  apt-get install libgeos-dev  --yes
# PYTHON PACKAGES
ADD ./docker/requirements.txt .
RUN pip3 install -r requirements.txt 


#KING
RUN wget https://www.kingrelatedness.com/Linux-king.tar.gz && tar -xzvf Linux-king.tar.gz -C /usr/local/bin && rm Linux*king*.tar.gz


ADD scripts ../scripts
ADD data ../data
