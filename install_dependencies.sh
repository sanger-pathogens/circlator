#!/bin/bash
set -e
set -x

start_dir=$(pwd)

BWA_VERSION=0.7.12
CANU_VERSION=1.4
PRODIGAL_VERSION=2.6.2
SAMTOOLS_VERSION=1.3
MUMMER_VERSION=3.23
SPADES_VERSION=3.11.1

MINIMAP2_VERSION=2.10
MINIASM_VERSION=55cf0189e2f7d5bda5868396cebe066eec0a9547
RACON_VERSION=1.3.0

BWA_DOWNLOAD_URL="http://downloads.sourceforge.net/project/bio-bwa/bwa-${BWA_VERSION}.tar.bz2"
CANU_DOWNLOAD_URL="https://github.com/marbl/canu/releases/download/v${CANU_VERSION}/canu-${CANU_VERSION}.Linux-amd64.tar.xz"
PRODIGAL_DOWNLOAD_URL="https://github.com/hyattpd/Prodigal/releases/download/v${PRODIGAL_VERSION}/prodigal.linux"
SAMTOOLS_DOWNLOAD_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2"
MUMMER_DOWNLOAD_URL="http://downloads.sourceforge.net/project/mummer/mummer/${MUMMER_VERSION}/MUMmer${MUMMER_VERSION}.tar.gz"
SPADES_DOWNLOAD_URL="http://spades.bioinf.spbau.ru/release${SPADES_VERSION}/SPAdes-${SPADES_VERSION}-Linux.tar.gz"
MINIMAP2_DOWNLOAD_URL="https://github.com/lh3/minimap2/archive/v${MINIMAP2_VERSION}.tar.gz"
MINIASM_DOWNLOAD_URL="https://github.com/lh3/miniasm/archive/${MINIASM_VERSION}.tar.gz"
RACON_DOWNLOAD_URL="https://github.com/isovic/racon/releases/download/${RACON_VERSION}/racon-v${RACON_VERSION}.tar.gz"

# Make an install location
if [ ! -d 'build' ]; then
  mkdir build
fi
cd build
build_dir=$(pwd)

# DOWNLOAD ALL THE THINGS
download () {
  url=$1
  download_location=$2

  if [ -e $download_location ]; then
    echo "Skipping download of $url, $download_location already exists"
  else
    echo "Downloading $url to $download_location"
    wget $url -O $download_location
  fi
}


# ------------- BWA -------------------
download $BWA_DOWNLOAD_URL "bwa-${BWA_VERSION}.tar.bz2"
bwa_dir=$build_dir/"bwa-${BWA_VERSION}"
tar -xjf bwa-${BWA_VERSION}.tar.bz2
cd $bwa_dir
make


# ------------- canu -------------------
cd $build_dir
download $CANU_DOWNLOAD_URL "canu-${CANU_VERSION}.tar.xz"
tar -xf canu-${CANU_VERSION}.tar.xz
canu_dir=$build_dir/canu-${CANU_VERSION}/Linux-amd64/bin/


# --------------- prodigal -----------------
cd $build_dir
prodigal_dir="$build_dir/prodigal-${PRODIGAL_VERSION}"
if [ ! -d "$prodigal_dir" ]; then
  mkdir $prodigal_dir
fi
cd $prodigal_dir
download $PRODIGAL_DOWNLOAD_URL prodigal
chmod 755 prodigal


# --------------- samtools -----------------
cd $build_dir
download $SAMTOOLS_DOWNLOAD_URL "samtools-${SAMTOOLS_VERSION}.tar.bz2"
samtools_dir="$build_dir/samtools-${SAMTOOLS_VERSION}"
tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2
cd $samtools_dir
make


# --------------- mummer ------------------
cd $build_dir
download $MUMMER_DOWNLOAD_URL "MUMmer${MUMMER_VERSION}.tar.gz"
mummer_dir="$build_dir/MUMmer${MUMMER_VERSION}"
tar -zxf MUMmer${MUMMER_VERSION}.tar.gz
cd $mummer_dir
make


# --------------- spades -----------------
cd $build_dir
download $SPADES_DOWNLOAD_URL "SPAdes-${SPADES_VERSION}-Linux.tar.gz"
spades_dir="$build_dir/SPAdes-${SPADES_VERSION}-Linux/bin"
tar -zxf SPAdes-${SPADES_VERSION}-Linux.tar.gz


# --------------- minimap2 -----------------
cd $build_dir
download $MINIMAP2_DOWNLOAD_URL "${MINIMAP2_VERSION}.tar.gz"
minimap2_dir="$build_dir/minimap2-${MINIMAP2_VERSION}/"
tar -zxf ${MINIMAP2_VERSION}.tar.gz
cd $minimap2_dir
make


# --------------- miniasm -----------------
cd $build_dir
download $MINIASM_DOWNLOAD_URL "${MINIASM_VERSION}.tar.gz"
miniasm_dir="$build_dir/miniasm-${MINIASM_VERSION}"
tar -zxf ${MINIASM_VERSION}.tar.gz
cd $miniasm_dir
make

# --------------- racon -----------------
cd $build_dir
download $RACON_DOWNLOAD_URL "racon-v${RACON_VERSION}.tar.gz"
racon_dir="$build_dir/racon-v${RACON_VERSION}"
tar -zxf racon-v${RACON_VERSION}.tar.gz
cd $racon_dir
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
racon_dir="$build_dir/racon-v${RACON_VERSION}/build/bin"


cd $start_dir

update_path () {
  new_dir=$1
  if [[ ! "$PATH" =~ (^|:)"${new_dir}"(:|$) ]]; then
    export PATH=${new_dir}:${PATH}
  fi
}

update_path ${bwa_dir}
update_path ${canu_dir}
update_path ${prodigal_dir}
update_path ${mummer_dir}
update_path ${samtools_dir}
update_path ${spades_dir}
update_path ${minimap2_dir}
update_path ${miniasm_dir}
update_path ${racon_dir}


