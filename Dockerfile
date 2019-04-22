# Circlator plus AMOS
#
FROM ubuntu:16.04

# Install the dependancies
RUN apt-get update && \
    apt-get -y upgrade && \
    	    apt-get install -y git wget unzip bzip2 xz-utils make g++ zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libpng-dev libssl-dev libboost-all-dev libstatistics-descriptive-perl libxml-parser-perl libdbi-perl python python3 python3-pip

RUN git clone https://github.com/sanger-pathogens/circlator.git
RUN cd circlator && ./install_dependencies.sh
ENV PATH /circlator/build/bwa-0.7.12:/circlator/build/canu-1.4/Linux-amd64/bin/:/circlator/build/prodigal-2.6.2:/circlator/build/samtools-1.3:/circlator/build/MUMmer3.23:/circlator/build/SPAdes-3.7.1-Linux/bin:$PATH
RUN cd circlator && python3 setup.py install

RUN cd /circlator/build && wget https://downloads.sourceforge.net/project/amos/amos/3.1.0/amos-3.1.0.tar.gz && \
    tar xf amos-3.1.0.tar.gz && \
        rm amos-3.1.0.tar.gz && \
	    cd amos-3.1.0 && \
	        ./configure && \
		    cp src/Align/find-tandem.cc src/Align/find-tandem.cc.original && \
		        echo '#include <getopt.h>' | cat - src/Align/find-tandem.cc.original > src/Align/find-tandem.cc && \
			    make && \
			        make install && \
				    make clean
ENV PATH="/circlator/build/amos-3.1.0/bin:${PATH}"
RUN export PATH
RUN ln -s /circlator/build/MUMmer3.23/show-coords /usr/local/bin
