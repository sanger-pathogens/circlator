# This container will install Circlator from master
#
FROM debian:testing

ENV   BUILD_DIR=/opt/circlator

# Install the dependancies
RUN   apt-get update -qq && apt-get install -y python python3-pip git wget unzip zlib1g-dev libncurses5-dev

RUN   apt-get install --yes libbz2-dev liblzma-dev

RUN   mkdir -p ${BUILD_DIR}
COPY  . ${BUILD_DIR}

# Require AMOS for the 'minimus' option.  Current version of AMOS won't
# build see (http://sourceforge.net/p/amos/bugs/; and no issue updates
# since 2014) *but* these locally-built debian packages may just do the job
COPY  amos_3.1.0-1_amd64.deb ${BUILD_DIR}/
RUN   cd ${BUILD_DIR} && apt-get install ./amos_3.1.0-1_amd64.deb
RUN   cd /usr/local/bin && find /usr/local/amos-3.1.0/bin -type f -executable -exec ln -s {} \;
# When installed from local debian packages, AMOS expects to find
# dependencies under /software/pathogen/external/apps, which fortunately
# mirrors the directory structure under /
RUN   mkdir -p /software/pathogen/external && ln -s / /software/pathogen/external/apps
# toAmos dependencies
RUN   apt-get install --yes cpanminus && cpanm XML::Parser
# rudeimentary test
RUN   toAmos -h   > /dev/null 2>&1; if [ 0 -ne $? ]; then echo 'toAmos test failed'   && exit 255; fi
RUN   minimus2 -h > /dev/null 2>&1; if [ 0 -ne $? ]; then echo 'minimus2 test failed' && exit 255; fi


RUN   cd ${BUILD_DIR} && ./install_dependencies.sh
ENV   PATH ${BUILD_DIR}/build/bwa-0.7.12:${BUILD_DIR}/build/canu-1.4/Linux-amd64/bin/:${BUILD_DIR}/build/prodigal-2.6.2:${BUILD_DIR}/build/samtools-1.3:${BUILD_DIR}/build/MUMmer3.23:${BUILD_DIR}/build/SPAdes-3.7.1-Linux/bin:$PATH
RUN   export PATH
RUN   cd ${BUILD_DIR} && python3 setup.py install

# some paths in AMOS are hardwired => symlink executables from /usr/bin
RUN   cd /usr/bin && \
      for EXECUTABLE in `find ${BUILD_DIR}/build -type f -executable | grep -v '/test'`; do ln -s ${EXECUTABLE}; done

RUN   circlator progcheck

CMD   echo "Usage:  docker run -v \`pwd\`:/var/data -it <IMAGE_NAME> bash" && \
      echo "" && \
      echo "This will place you in a shell with your current working directory accessible as /var/data." && \
      echo "You can then run commands like:" && \
      echo "   circlator all /var/data/assembly.fasta /var/data/reads /var/data/<output_subdirectory>" && \
      echo "For help, please go to https://github.com/sanger-pathogens/circlator/wiki, or type" && \
      echo "   circlator --help"
