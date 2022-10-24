FROM eu.gcr.io/finngen-refinery-dev/bioinformatics:0.6

ENV TZ=Europe/Helsinki
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# I NEED PYTHON 3.7 (PANDAS ISSUES)
ARG PY_VER=3.7

RUN apt-get update && apt-get upgrade --yes
RUN apt-get install libproj-dev proj-data proj-bin --yes &&  apt-get install libgeos-dev libjpeg-dev zlib1g-dev  --yes
RUN apt install python${PY_VER} python${PY_VER}-dev  --yes
RUN pip3 install --upgrade pip

# R RELATED STUFF
ADD ./docker/install_packages.R /usr/local/bin/
RUN chmod a+x /usr/local/bin/install_packages.R && install_packages.R


# PYTHON PACKAGES
ADD ./docker/requirements.txt .
RUN python${PY_VER} -m pip install cython
RUN python${PY_VER} -m pip install numpy
RUN python${PY_VER} -m pip install -r requirements.txt

#KING
RUN wget https://www.kingrelatedness.com/Linux-king.tar.gz && tar -xzvf Linux-king.tar.gz -C /usr/local/bin && rm Linux*king*.tar.gz

ADD scripts ../scripts
ADD data ../data
