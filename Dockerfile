From ubuntu:xenial
LABEL MAINTAINER "steven_vensko@med.unc.edu"

RUN apt-get update
RUN apt-get install -y git
RUN apt-get install -y wget
RUN apt-get install -y build-essential
RUN apt-get install -y libbz2-dev
RUN apt-get install -y xz-utils
RUN apt-get install -y liblzma-dev
RUN apt-get install -y bzip2
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y ncurses-dev
RUN apt-get install -y python2.7
RUN apt-get install -y python2.7-dev


RUN wget https://bootstrap.pypa.io/pip/2.7/get-pip.py
RUN python2.7 get-pip.py

RUN pip install numpy==1.16.0
RUN pip install networkx==1.11
RUN pip install pyahocorasick==1.4.0
RUN pip install bcbio-gff==0.6.4
RUN pip install biopython==1.70
RUN pip install pyfaidx==0.5.3.1
RUN pip install pysam==0.14.1
RUN pip install pandas==0.23.4
RUN pip install scipy==1.2.0

RUN wget https://github.com/holtjma/msbwt/archive/0.3.0.tar.gz \
    && tar -xzvf 0.3.0.tar.gz \
    && cd msbwt-0.3.0 \
    && python2.7 setup.py install \
    && cd /

RUN git clone https://github.com/holtjma/msbwt-is.git \
    && cd msbwt-is \
    && make \
    && ./msbwtis -h \
    && cd /

# Download and compile samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2
RUN bunzip2 samtools-1.7.tar.bz2
RUN tar xf samtools-1.7.tar
RUN cd samtools-1.7 && make

RUN mkdir /netMHCpan-4.0-docker
RUN mkdir /netMHCIIpan-3.2-docker
COPY stage/netMHCpan-4.0-docker /netMHCpan-4.0-docker
COPY stage/netMHCIIpan-3.2-docker /netMHCIIpan-3.2-docker

ENV PATH /samtools-1.7/:${PATH}
RUN chmod -R ugo+rx /netMHCpan-4.0-docker/
RUN chmod -R ugo+rx /netMHCIIpan-3.2-docker/

COPY NeoSplice /NeoSplice
RUN chmod -R ugo+rx /NeoSplice/

RUN ln -s /usr/bin/python2.7 /usr/bin/python
