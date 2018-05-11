# This container will install Circlator from master
#
FROM debian:testing

# Install the dependancies
RUN apt-get update -qq && apt-get install -y python python3-pip git wget unzip zlib1g-dev libncurses5-dev

RUN git clone https://github.com/sanger-pathogens/circlator.git
RUN cd circlator && ./install_dependencies.sh
ENV PATH /circlator/build/bwa-0.7.12:/circlator/build/canu-1.4/Linux-amd64/bin/:/circlator/build/prodigal-2.6.2:/circlator/build/samtools-1.3:/circlator/build/MUMmer3.23:/circlator/build/SPAdes-3.7.1-Linux/bin:$PATH
RUN export PATH
RUN cd circlator && python3 setup.py install
