# Circlator plus AMOS
# h/t Yuri Bendaña (ybendana) for working AMOS install
FROM ubuntu:16.04

ENV   BUILD_DIR=/opt/circlator

# Install the dependancies
RUN   apt-get update && \
      apt-get --yes upgrade && \
      apt-get install --yes apt-utils && \
      apt-get install --yes   git wget unzip bzip2 xz-utils make g++ zlib1g-dev libncurses5-dev libbz2-dev \
                              liblzma-dev libcurl4-openssl-dev libpng-dev libssl-dev libboost-all-dev \
                              libstatistics-descriptive-perl libxml-parser-perl libdbi-perl \
                              python python3 python3-pip default-jre-headless

RUN   apt-get install -y locales && \
      sed -i -e 's/# \(en_GB\.UTF-8 .*\)/\1/' /etc/locale.gen && \
      touch /usr/share/locale/locale.alias && \
      locale-gen
ENV   LANG     en_GB.UTF-8
ENV   LANGUAGE en_GB:en
ENV   LC_ALL   en_GB.UTF-8

RUN   mkdir -p ${BUILD_DIR}
COPY  . ${BUILD_DIR}
RUN   cd ${BUILD_DIR} && ./install_dependencies.sh

ENV   PATH ${BUILD_DIR}/build/bwa-0.7.12:${BUILD_DIR}/build/canu-1.4/Linux-amd64/bin/:${BUILD_DIR}/build/prodigal-2.6.2:${BUILD_DIR}/build/samtools-1.3:${BUILD_DIR}/build/MUMmer3.23:${BUILD_DIR}/build/SPAdes-3.7.1-Linux/bin:$PATH

RUN   cd ${BUILD_DIR} && python3 setup.py install

RUN   cd ${BUILD_DIR}/build && \
      wget -O- https://downloads.sourceforge.net/project/amos/amos/3.1.0/amos-3.1.0.tar.gz | tar xzf - && \
      cd amos-3.1.0 && \
      ./configure && \
      cp src/Align/find-tandem.cc src/Align/find-tandem.cc.original && \
      echo '#include <getopt.h>' | cat - src/Align/find-tandem.cc.original > src/Align/find-tandem.cc && \
      make && make install && make clean

ENV   PATH="${BUILD_DIR}/build/amos-3.1.0/bin:${PATH}"
RUN   export PATH
RUN   ln -s ${BUILD_DIR}/build/MUMmer3.23/show-coords /usr/local/bin

RUN   circlator progcheck && circlator test /tmp/circlator-test && rm -r /tmp/circlator-test

RUN   TEST_DATA_DIR=`ls -d /usr/local/lib/python3.?/dist-packages/circlator* | head -1`/circlator/data; \
      if [ -f ${TEST_DATA_DIR}/test_contigs.fa ]; \
         then  chmod a+rX ${TEST_DATA_DIR}/*; \
         else  echo "Could not find test data directory ${TEST_DATA_DIR}" && exit 255; \
      fi

CMD   echo "Usage:  docker run -v \`pwd\`:/var/data -it <IMAGE_NAME> bash" && \
      echo "" && \
      echo "This will place you in a shell with your current working directory accessible as /var/data." && \
      echo "You can then run commands like:" && \
      echo "   circlator all /var/data/assembly.fasta /var/data/reads /var/data/<output_subdirectory>" && \
      echo "For help, please go to https://github.com/sanger-pathogens/circlator/wiki, or type" && \
      echo "   circlator --help"
