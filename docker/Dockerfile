FROM eu.gcr.io/finngen-refinery-dev/bioinformatics:0.6

ENV TZ=Europe/Helsinki
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone


RUN apt-get update && apt-get upgrade --yes
RUN apt-get install parallel ghostscript libproj-dev proj-data proj-bin --yes &&  apt-get install libgeos-dev libjpeg-dev zlib1g-dev  --yes


# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

ENV PATH=$CONDA_DIR/bin:$PATH

# PYTHON PACKAGES
ADD ./docker/requirements.txt .
RUN conda config --append channels conda-forge
RUN conda install --file requirements.txt --yes
 
# R RELATED STUFF
ADD ./docker/install_packages.R /usr/local/bin/
RUN chmod a+x /usr/local/bin/install_packages.R && install_packages.R

#KING
RUN wget https://www.kingrelatedness.com/Linux-king.tar.gz && tar -xzvf Linux-king.tar.gz -C /usr/local/bin && rm Linux*king*.tar.gz

ADD scripts ../scripts
ADD data ../data
