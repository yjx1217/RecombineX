#!/bin/bash
# last update: 2021/03/04

set -e -o pipefail

#########################
RECOMBINEX_HOME=$(pwd)
BUILD="build"
mainland_china_installation="no";
#########################

timestamp () {
  date +"%F %T"
}

clean () {
    dir=$1
    if [ -d $dir ] 
    then
	echo "remove previously failed installation in $BUILD/$dir"
	rm -rf $dir
    fi
}

clone () {
  url=$1
  dir=$(basename $url)
  echo "run clone for \"git clone $url\""
  git clone $url --depth 1
  cd $dir
  git fetch --unshallow
}

download () {
  url=$1
  download_location=$2
  echo "Downloading $url to $download_location"
  wget -c --no-check-certificate $url -O $download_location
}

tidy_version () { 
    echo "$1" | awk -F. '{ printf("%d%03d%03d%03d\n", $1,$2,$3,$4); }';
}

check_installed () {
    if [ -e "$1/installed" ]; then
        echo "installed"
    else
        echo ""
    fi
}

note_installed () {
    touch "$1/installed"
}

echo ""
echo ""
echo "##################################################################"
echo "###                                                            ###"
echo "###                  Welcome to RecombineX                     ###"
echo "###                                                            ###"
echo "##################################################################"
echo ""
echo ""
echo "[$(timestamp)] Installation starts ..."
echo ""

if [ -z "$MAKE_JOBS" ]
then
    echo "[$(timestamp)] Defaulting to 2 concurrent jobs when executing make. Override with MAKE_JOBS=<NUM>"
    MAKE_JOBS=2
    echo ""
fi

while getopts ":hc" opt
do
    case "${opt}" in
        h)
            echo "Usage:"
            echo "bash install_dependencies.sh"
            echo "When installing within mainland China, please run this script with the '-c' option >"
            echo "bash install_dependencies.sh -c";;
        c)
            echo "Detected the '-c' option >"
            echo "Set installation location as 'mainland_china'" 
            mainland_china_installation="yes";;
    esac
done

BLAST_VERSION="2.2.31" #
BLAST_DOWNLOAD_URL="http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz"

RMBLAST_VERSION="2.2.28" #
RMBLAST_DOWNLOAD_URL="http://ftp.ncbi.nlm.nih.gov/blast/executables/rmblast/${RMBLAST_VERSION}/ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz"

MUMMER3_VERSION="3.23" # released on 2011.12.17
MUMMER3_DOWNLOAD_URL="https://sourceforge.net/projects/mummer/files/mummer/${MUMMER3_VERSION}/MUMmer${MUMMER3_VERSION}.tar.gz"

# MUMMER4_VERSION="4.0.0beta2" # released on 2017.10.14
# MUMMER4_DOWNLOAD_URL="https://github.com/gmarcais/mummer/releases/download/v${MUMMER4_VERSION}/mummer-${MUMMER4_VERSION}.tar.gz"

# ASSEMBLYTICS_GITHUB_COMMIT_VERSION="df5361f" # committed on 2017.11.02
# ASSEMBLYTICS_DOWNLOAD_URL="https://github.com/MariaNattestad/Assemblytics.git"

GNUPLOT_VERSION="4.6.6"
GNUPLOT_DOWNLOAD_URL="https://sourceforge.net/projects/gnuplot/files/gnuplot/${GNUPLOT_VERSION}/gnuplot-${GNUPLOT_VERSION}.tar.gz"

BEDTOOLS_VERSION="2.27.1" # released on 2017.12.14
BEDTOOLS_DOWNLOAD_URL="https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz"

SRA_VERSION="2.9.6" # released on 2019.03.18
SRA_DOWNLOAD_URL="http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz"

ART_VERSION="mountrainier2016.06.05" # released on 2016.06.05
ART_DOWNLOAD_URL="https://www.niehs.nih.gov/research/resources/assets/docs/artbin${ART_VERSION}linux64.tgz"

TRIMMOMATIC_VERSION="0.38" #
TRIMMOMATIC_DOWNLOAD_URL="http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${TRIMMOMATIC_VERSION}.zip"

BWA_VERSION="0.7.17" # released on 2017.10.23
BWA_DOWNLOAD_URL="https://github.com/lh3/bwa/releases/download/v${BWA_VERSION}/bwa-${BWA_VERSION}.tar.bz2"

SAMTOOLS_VERSION="1.9" # released on 2018.07.18
SAMTOOLS_DOWNLOAD_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2"

GATK3_VERSION="3.6-6" #
GATK3_DOWNLOAD_URL="https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/GenomeAnalysisTK.jar"

# GATK4_VERSION="4.0.11.0" # released on 2018.10.23
# GATK4_DOWNLOAD_URL="https://github.com/broadinstitute/gatk/releases/download/${GATK4_VERSION}/gatk-${GATK4_VERSION}.zip"

FREEBAYES_VERSION="1.3.4" # released on 20210129
FREEBAYES_SOURCE_DOWNLOAD_URL="https://github.com/freebayes/freebayes/releases/download/v1.3.4/freebayes-1.3.4-src.tar.gz" 
FREEBAYES_DOWNLOAD_URL="https://github.com/freebayes/freebayes/releases/download/v${FREEBAYES_VERSION}/freebayes-${FREEBAYES_VERSION}-linux-static-AMD64.gz"

PICARD_VERSION="2.19.0" # released on 2019.03.22 
PICARD_DOWNLOAD_URL="https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar"

GEMTOOLS_VERSION="1.7.1" # released on 2013.11.18
GEMTOOLS_DOWNLOAD_URL="http://barnaserver.com/gemtools/releases/GEMTools-static-i3-${GEMTOOLS_VERSION}.tar.gz"

VCFLIB_VERSION="1.0.1" # released on 2019.10.1
VCFLIB_DOWNLOAD_URL="https://github.com/vcflib/vcflib/releases/download/v${VCFLIB_VERSION}/vcflib-${VCFLIB_VERSION}-src.tar.gz"

# VCFTOOLS_VERSION="0.1.16"
# VCFTOOLS_DOWNLOAD_URL="https://github.com/vcftools/vcftools/releases/download/v${VCFTOOLS_VERSION}/vcftools-${VCFTOOLS_VERSION}.tar.gz"

VT_VERSION="" # 
VT_GITHUB_COMMIT_VERSION="f6d2b5d" # committed on 2018.08.01
VT_DOWNLOAD_URL="https://github.com/atks/vt"

FREEC_VERSION="11.4" # released on 2018.04.27
FREEC_DOWNLOAD_URL="https://github.com/BoevaLab/FREEC/archive/v${FREEC_VERSION}.tar.gz"

PARALLEL_VERSION="20180722" # released on 2018.07.22
PARALLEL_DOWNLOAD_URL="http://ftp.gnu.org/gnu/parallel/parallel-${PARALLEL_VERSION}.tar.bz2"

# UCSC Utilities
FASPLIT_DOWNLOAD_URL="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSplit"
FATOTWOBIT_DOWNLOAD_URL="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit"
TWOBITINFO_DOWNLOAD_URL="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo"


if [ -d $BUILD ]
then
    echo ""
    echo "[$(timestamp)] Detected previously generated $BUILD directory."
else
    echo "[$(timestamp)] Create the new $BUILD directory."
    mkdir $BUILD
    echo ""
fi

cd $BUILD
build_dir=$(pwd)

# Downloading all the dependencies
echo ""
echo "[$(timestamp)] Download and install all the dependencies ..."

# ---------- set Perl & R environment variables -------------
#PYTHONPATH="$build_dir"
PERL5LIB="$build_dir:$PERL5LIB"
PERL5LIB="$build_dir/cpanm/perlmods/lib/perl5:$PERL5LIB"
R_LIBS="$build_dir/R_libs:$R_LIBS"
echo ""
echo "[$(timestamp)] Installing Perl modules ..."
cpanm_dir="$build_dir/cpanm"
if [ -z $(check_installed $cpanm_dir) ]; then
    clean $cpanm_dir
    mkdir -p $cpanm_dir
    cd $cpanm_dir

    wget -c --no-check-certificate -O - https://cpanmin.us/ > cpanm
    chmod +x cpanm
    mkdir perlmods

    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Test::More@1.302086
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Text::Soundex@3.05
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Env@1.04
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Statistics::Descriptive@3.0612
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Statistics::Descriptive::Discrete@0.07
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Math::Random@0.72
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Math::Round@0.07
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Sys::Syslog@0.35

    note_installed $cpanm_dir
fi    

echo ""
echo "[$(timestamp)] Installing R libraries ..."
rlib_dir="$build_dir/R_libs"
mkdir -p $rlib_dir
cd $rlib_dir
R_VERSION=$(R --version |head -1 |cut -d " " -f 3)

if [ -z $(check_installed "$rlib_dir/optparse") ]; then
    clean "$rlib_dir/optparse"
    R -e "install.packages(\"optparse\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\")"
    note_installed "$rlib_dir/optparse"
fi

if [ -z $(check_installed "$rlib_dir/ggplot2") ]; then
    clean "$rlib_dir/ggplot2"
    R -e "install.packages(\"ggplot2\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\")"
    note_installed "$rlib_dir/ggplot2"
fi

if [ -z $(check_installed "$rlib_dir/scales") ]; then
    clean "$rlib_dir/scales"
    R -e "install.packages(\"scales\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\")"
    note_installed "$rlib_dir/scales"
fi

if [ -z $(check_installed "$rlib_dir/viridis") ]; then
    clean "$rlib_dir/viridis"
    R -e "install.packages(\"viridis\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\")"
    note_installed "$rlib_dir/viridis"
fi

if [ -z $(check_installed "$rlib_dir/DNAcopy") ]; then
    clean "$rlib_dir/DNAcopy"
    if [ $(tidy_version "$R_VERSION") -ge $(tidy_version "3.6.0") ]
    then
	echo "R_VERSION=$R_VERSION, use the new bioconductor installation protocol"
	R -e ".libPaths(\"$build_dir/R_libs/\");install.packages(\"BiocManager\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\");BiocManager::install(\"DNAcopy\")"
    else
	echo "R_VERSION=$R_VERSION, use the old bioconductor installation protocol"
	R -e ".libPaths(\"$build_dir/R_libs/\");source(\"https://bioconductor.org/biocLite.R\");biocLite(\"DNAcopy\", type = \"source\")"
    fi
    note_installed "$rlib_dir/DNAcopy"
fi


# --------------- ncbi-blast+ ------------------
echo ""
echo "[$(timestamp)] Installing ncbi-blast+ ..."
blast_dir="$build_dir/ncbi-blast-${BLAST_VERSION}+/bin"
windowmasker_dir="$build_dir/ncbi-blast-${BLAST_VERSION}+/bin"
if [ -z $(check_installed $blast_dir) ]; then
    cd $build_dir
    clean "$build_dir/ncbi-blast-${BLAST_VERSION}+"
    echo ""
    echo "Installing ncbi-blast+ ..."
    echo "Download ncbi-blast-v${BLAST_VERSION}"
    download $BLAST_DOWNLOAD_URL "ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz"
    tar xvzf ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz
    cd $build_dir
    rm ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz
    note_installed $blast_dir
fi    

# --------------- ncbi-rmblast ------------------
echo ""
echo "[$(timestamp)] Installing ncbi-rmblast ..."
rmblast_dir="$build_dir/ncbi-rmblastn-${RMBLAST_VERSION}/bin"
if [ -z $(check_installed $rmblast_dir) ]; then
    cd $build_dir
    clean "$build_dir/ncbi-rmblastn-${RMBLAST_VERSION}"
    echo ""
    echo "Installing ncbi-rmblastn ..."
    echo "Download ncbi-rmblastn-v${BLAST_VERSION}"
    download $RMBLAST_DOWNLOAD_URL "ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz"
    tar xvzf ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz
    # copy rmblastn binary file to ncbi-blast+ directory for easy RepeatMasker configuration
    cp $rmblast_dir/rmblastn $blast_dir
    cd $build_dir
    rm ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz
    note_installed $rmblast_dir
fi

# --------------- mummer4 ------------------
# mummer4_dir="$build_dir/mummer-${MUMMER4_VERSION}"
# if [ -z $(check_installed $mummer4_dir) ]; then
#     cd $build_dir
#     echo "Download mummer-v${MUMMER4_VERSION}"
#     download $MUMMER4_DOWNLOAD_URL "mummer-${MUMMER4_VERSION}.tar.gz"
#     tar -xzf mummer-${MUMMER4_VERSION}.tar.gz 
#     echo "$mummer4_dir"
#     cd $mummer4_dir
#     ./configure
#     make -j $MAKE_JOBS
#     note_installed $mummer4_dir
#     cd ..
#     rm mummer-${MUMMER4_VERSION}.tar.gz 
#     note_installed $mummer4_dir
# fi
# PATH="$mummer4_dir:${PATH}"

# --------------- mummer3 ------------------
echo ""
echo "[$(timestamp)] Installing mummer3 ..."
mummer3_dir="$build_dir/MUMmer${MUMMER3_VERSION}"
if [ -z $(check_installed $mummer3_dir) ]; then
    cd $build_dir
    clean "$build_dir/MUMmer${MUMMER3_VERSION}"
    echo ""
    echo "Installing mummer3 ..."
    echo "Download mummer-v${MUMMER3_VERSION}"
    download $MUMMER3_DOWNLOAD_URL "MUMmer${MUMMER3_VERSION}.tar.gz"
    tar xvzf MUMmer${MUMMER3_VERSION}.tar.gz
    cd $mummer3_dir
    make check
    make -j $MAKE_JOBS
    PATH=$mummer3_dir:$PATH
    cd $build_dir
    rm MUMmer${MUMMER3_VERSION}.tar.gz
    note_installed $mummer3_dir
fi

# # --------------- ASSEMBLYTICS -----------------
# cd $build_dir
# echo "Download ASSEMBLYTICS-v${ASSEMBLYTICS_VERSION}"
# git clone $ASSEMBLYTICS_DOWNLOAD_URL
# assemblytics_dir="$build_dir/Assemblytics"
# cd $assemblytics_dir
# git checkout -f -q $ASSEMBLYTICS_GITHUB_COMMIT_VERSION

# --------------- bedtools ------------------
echo ""
echo "[$(timestamp)] Installing bedtools ..."
bedtools_dir="$build_dir/bedtools2/bin"
if [ -z $(check_installed $bedtools_dir) ]; then
    cd $build_dir
    clean "$build_dir/bedtools2"
    echo "Download bedtools-v${BEDTOOLS_VERSION}"
    download $BEDTOOLS_DOWNLOAD_URL "bedtools-${BEDTOOLS_VERSION}.tar.gz"
    tar xvzf bedtools-${BEDTOOLS_VERSION}.tar.gz
    cd "$build_dir/bedtools2"
    make -j $MAKE_JOBS
    cd $build_dir
    rm bedtools-${BEDTOOLS_VERSION}.tar.gz
    note_installed $bedtools_dir
fi

# ------------- SRA Toolkit -------------------
echo ""
echo "[$(timestamp)] Installing SRA Toolkit ..."
sra_dir="$build_dir/sratoolkit.${SRA_VERSION}-centos_linux64/bin"
if [ -z $(check_installed $sra_dir) ]; then
    cd $build_dir
    clean "$build_dir/sratoolkit.${SRA_VERSION}-centos_linux64"
    echo "Download SRAtoolkit-v${SRA_VERSION}"
    download $SRA_DOWNLOAD_URL sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz
    tar xvzf sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz
    rm sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz
    note_installed $sra_dir
fi

# ------------- ART -------------------
echo ""
echo "[$(timestamp)] Installing ART ..."
art_dir="$build_dir/art_bin_MountRainier"
if [ -z $(check_installed $art_dir) ]; then
    cd $build_dir
    clean "$build_dir/art_bin_MountRainier"
    echo "Download ART-v${ART_VERSION}"
    download $ART_DOWNLOAD_URL artbin${ART_VERSION}linux64.tgz
    tar -zxf artbin${ART_VERSION}linux64.tgz
    rm artbin${ART_VERSION}linux64.tgz
    note_installed $art_dir
fi

# --------------- Trimmomatic -----------------
echo ""
echo "[$(timestamp)] Installing Trimmomatic ..."
trimmomatic_dir="$build_dir/Trimmomatic-${TRIMMOMATIC_VERSION}"
if [ -z $(check_installed $trimmomatic_dir) ]; then
    cd $build_dir
    clean "$build_dir/Trimmomatic-${TRIMMOMATIC_VERSION}"
    echo "Download Trimmomatic-v${TRIMMOMATIC_VERSION}"
    download $TRIMMOMATIC_DOWNLOAD_URL "Trimmomatic-${TRIMMOMATIC_VERSION}.zip"
    unzip Trimmomatic-${TRIMMOMATIC_VERSION}.zip
    cd $trimmomatic_dir
    chmod 755 trimmomatic-${TRIMMOMATIC_VERSION}.jar
    ln -s trimmomatic-${TRIMMOMATIC_VERSION}.jar trimmomatic.jar 
    cd $build_dir
    rm Trimmomatic-${TRIMMOMATIC_VERSION}.zip
    note_installed $trimmomatic_dir
fi

# ------------- BWA -------------------
echo ""
echo "[$(timestamp)] Installing BWA ..."
bwa_dir="$build_dir/bwa-${BWA_VERSION}"
if [ -z $(check_installed $bwa_dir) ]; then
    cd $build_dir
    clean "$build_dir/bwa-${BWA_VERSION}"
    echo "Download BWA-v${BWA_VERSION}"
    download $BWA_DOWNLOAD_URL "bwa-${BWA_VERSION}.tar.bz2"
    tar xvjf bwa-${BWA_VERSION}.tar.bz2
    cd $bwa_dir
    make -j $MAKE_JOBS
    cd $build_dir
    rm bwa-${BWA_VERSION}.tar.bz2
    note_installed $bwa_dir
fi

# --------------- samtools -----------------
echo ""
echo "[$(timestamp)] Installing samtools ..."
samtools_dir="$build_dir/samtools-${SAMTOOLS_VERSION}"
htslib_dir="$samtools_dir/htslib-${SAMTOOLS_VERSION}"
tabix_dir="$samtools_dir/htslib-${SAMTOOLS_VERSION}"
if [ -z $(check_installed $samtools_dir) ]; then
    cd $build_dir
    clean "$build_dir/samtools-${SAMTOOLS_VERSION}"
    echo "Download samtools-v${SAMTOOLS_VERSION}"
    download $SAMTOOLS_DOWNLOAD_URL samtools-${SAMTOOLS_VERSION}.tar.bz2
    tar xvjf samtools-${SAMTOOLS_VERSION}.tar.bz2
    cd $samtools_dir
    C_INCLUDE_PATH=""
    ./configure --without-curses;
    make -j $MAKE_JOBS
    cd htslib-${SAMTOOLS_VERSION}
    #autoheader
    #autoconf
    ./configure
    make -j $MAKE_JOBS
    cd $build_dir
    rm samtools-${SAMTOOLS_VERSION}.tar.bz2
    note_installed $samtools_dir
fi
PATH="$samtools_dir:$htslib_dir:$tabix_dir:${PATH}"

# --------------- Picard -----------------
echo ""
echo "[$(timestamp)] Installing picard ..."
picard_dir="$build_dir/Picard-v${PICARD_VERSION}"
if [ -z $(check_installed $picard_dir) ]; then
    cd $build_dir
    clean "$build_dir/Picard-v${PICARD_VERSION}"
    echo "Download Picard-v${PICARD_VERSION}"
    download $PICARD_DOWNLOAD_URL "picard.jar"
    mkdir Picard-v${PICARD_VERSION}
    mv picard.jar $picard_dir
    cd $picard_dir
    chmod 755 picard.jar
    note_installed $picard_dir
fi

# --------------- GATK3 ------------------
echo ""
echo "[$(timestamp)] Installing GATK3 ..."
gatk3_dir="$build_dir/GATK3"
if [ -z $(check_installed $gatk3_dir) ]; then
    cd $build_dir
    clean "$build_dir/GATK3"
    mkdir GATK3
    cd GATK3
    download "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/GenomeAnalysisTK.jar" GenomeAnalysisTK.jar
    chmod 755 GenomeAnalysisTK.jar
    ln -s GenomeAnalysisTK.jar gatk3.jar
    note_installed $gatk3_dir
fi

# # --------------- GATK4 ------------------
# echo ""
# echo "[$(timestamp)] Installing GATK4 ..."                                                                                                
# gatk4_dir="$build_dir/GATK4"
# if [ -z $(check_installed $gatk4_dir) ]; then
#     cd $build_dir
#     clean "$build_dir/GATK4"
#     echo "Download GATK4-v${GATK4_VERSION}"
#     download $GATK4_DOWNLOAD_URL "gatk-${GATK4_VERSION}.zip"
#     unzip gatk-${GATK4_VERSION}.zip
#     mv gatk-${GATK4_VERSION} GATK4
#     cd GATK4
#     ln -s gatk-package-${GATK4_VERSION}-local.jar gatk4.jar
#     cd $build_dir
#     rm gatk-${GATK4_VERSION}.zip
#     note_installed $gatk4_dir
# fi

# --------------- freebayes ------------------
echo ""
echo "[$(timestamp)] Installing Freebayes ..."                                                                                                
freebayes_dir="$build_dir/freebayes/bin"
if [ -z $(check_installed $freebayes_dir) ]; then
    cd $build_dir
    clean "$build_dir/freebayes/bin"
    echo "Download freebayes-v${FREEBAYES_VERSION}"
    download $FREEBAYES_SOURCE_DOWNLOAD_URL freebayes-${FREEBAYES_VERSION}-src.tar.gz
    tar xzf freebayes-${FREEBAYES_VERSION}-src.tar.gz
    rm freebayes-${FREEBAYES_VERSION}-src.tar.gz
    cd freebayes
    cd bin
    download $FREEBAYES_DOWNLOAD_URL "freebayes-${FREEBAYES_VERSION}-linux-static-AMD64.gz"
    gunzip freebayes-${FREEBAYES_VERSION}-linux-static-AMD64.gz
    chmod 755 freebayes-${FREEBAYES_VERSION}-linux-static-AMD64
    ln -s freebayes-${FREEBAYES_VERSION}-linux-static-AMD64 freebayes
    cd $build_dir
    note_installed $freebayes_dir
fi

# --------------- GEM-Tools -----------------
echo ""
echo "[$(timestamp)] Installing GEM-Tools ..."
gemtools_dir="$build_dir/gemtools-${GEMTOOLS_VERSION}-i3/bin"
if [ -z $(check_installed $gemtools_dir) ]; then
    cd $build_dir
    clean "$build_dir/gemtools-${GEMTOOLS_VERSION}-i3"
    echo "Download GEMTOOLS-v${GEMTOOLS_VERSION}"
    download $GEMTOOLS_DOWNLOAD_URL "gemtools-${GEMTOOLS_VERSION}.tar.gz"
    tar zxf "gemtools-${GEMTOOLS_VERSION}.tar.gz"
    cd $build_dir
    rm gemtools-${GEMTOOLS_VERSION}.tar.gz
    note_installed $gemtools_dir
fi

# --------------- VCFlib -----------------
echo ""
echo "[$(timestamp)] Installing VCFlib ..."
vcflib_dir="$build_dir/vcflib-${VCFLIB_VERSION}-src/bin"
if [ -z $(check_installed $vcflib_dir) ]; then
    cd $build_dir
    clean "$build_dir/vcflib-${VCFLIB_VERSION}-src"
    echo "Download VCFlib-v${VCFLIB_VERSION}"
    download $VCFLIB_DOWNLOAD_URL vcflib-${VCFLIB_VERSION}-src.tar.gz 
    tar xvzf vcflib-${VCFLIB_VERSION}-src.tar.gz
    cd $build_dir/vcflib-${VCFLIB_VERSION}-src
    make 
    cd $build_dir
    rm vcflib-${VCFLIB_VERSION}-src.tar.gz
    note_installed $vcflib_dir
fi

# # --------------- vcftools -----------------
# echo ""
# echo "[$(timestamp)] Installing vcftools ..."
# vcftools_dir="$build_dir/vcftools-${VCFTOOLS_VERSION}/build/bin"
# if [ -z $(check_installed $vcflib_dir) ]; then
#     cd $build_dir
#     clean "$build_dir/vcftools-${VCFTOOLS_VERSION}"
#     echo "Download vcftools-v${VCFTOOLS_VERSION}"
#     download $VCFTOOLS_DOWNLOAD_URL "vcftools-${VCFTOOLS_VERSION}.tar.gz"
#     tar xvzf vcftools-${VCFTOOLS_VERSION}.tar.gz
#     cd vcftools-${VCFTOOLS_VERSION}
#     mkdir build
#     ./configure --prefix="$build_dir/vcftools-${VCFTOOLS_VERSION}/build"
#     make -j $MAKE_JOBS
#     make install
#     cd $build_dir
#     rm vcftools-${VCFTOOLS_VERSION}.tar.gz
#     note_installed $vcftools_dir
# fi

# ------------------ VT ---------------------
echo ""
echo "[$(timestamp)] Installing VT ..."
vt_dir="$build_dir/vt"
if [ -z $(check_installed $vt_dir) ]; then
    cd $build_dir
    clean "$build_dir/vt"
    echo "Download VT-v${VT_VERSION}"
    # git clone $VT_DOWNLOAD_URL
    clone $VT_DOWNLOAD_URL
    cd $vt_dir
    git checkout -f -q $VT_GITHUB_COMMIT_VERSION
    make -j $MAKE_JOBS
    make test
    note_installed $vt_dir
fi

# ----------------- FREEC --------------------
echo ""
echo "[$(timestamp)] Installing FREEC ..."
freec_dir="$build_dir/FREEC-${FREEC_VERSION}"
if [ -z $(check_installed $freec_dir) ]; then
    cd $build_dir
    clean "$build_dir/FREEC-${FREEC_VERSION}"
    echo "Download FREEC-v${FREEC_VERSION}"
    download $FREEC_DOWNLOAD_URL "FREEC-${FREEC_VERSION}.tar.gz"
    tar xvzf FREEC-${FREEC_VERSION}.tar.gz
    cd $freec_dir/src
    make -j $MAKE_JOBS
    cp freec $freec_dir
    cd $freec_dir/scripts
    chmod 755 *
    cp * $freec_dir
    cd $build_dir
    rm FREEC-${FREEC_VERSION}.tar.gz
    note_installed $freec_dir
fi

# --------------- parallel ------------------                                                                                                                        
echo ""
echo "[$(timestamp)] Installing parallel ..."
parallel_dir="$build_dir/parallel-${PARALLEL_VERSION}/bin"
if [ -z $(check_installed $parallel_dir) ]; then
    cd $build_dir
    clean "$build_dir/parallel-${PARALLEL_VERSION}"
    echo "Download parallel-${PARALLEL_VERSION}"
    download $PARALLEL_DOWNLOAD_URL "parallel_v${PARALLEL_VERSION}.tar.bz2"
    tar xvjf parallel_v${PARALLEL_VERSION}.tar.bz2
    cd parallel-${PARALLEL_VERSION}
    ./configure --prefix="$build_dir/parallel-${PARALLEL_VERSION}"
    make -j $MAKE_JOBS
    make install
    parallel_dir="$build_dir/parallel-${PARALLEL_VERSION}/bin"
    cd ..
    rm parallel_v${PARALLEL_VERSION}.tar.bz2
    note_installed $parallel_dir
fi

# --------------- UCSC Utilities -----------------
echo ""
echo "[$(timestamp)] Installing UCSC Utilities ..."
ucsc_dir="$build_dir/UCSC_Utilities"
if [ -z $(check_installed $ucsc_dir) ]; then
    cd $build_dir
    clean "$build_dir/UCSC_Utilities"
    mkdir UCSC_Utilities
    cd $ucsc_dir
    download $FASPLIT_DOWNLOAD_URL "faSplit"
    download $FATOTWOBIT_DOWNLOAD_URL "faToTwoBit"
    download $TWOBITINFO_DOWNLOAD_URL "twoBitInfo"
    chmod 755 $ucsc_dir/*
    note_installed $ucsc_dir
fi


# Configure executable paths

cd $RECOMBINEX_HOME
echo ""
echo "[$(timestamp)] Configuring executable paths ..."
echo "export RECOMBINEX_HOME=${RECOMBINEX_HOME}" > env.sh
echo "export build_dir=${build_dir}" >> env.sh
echo "export PERL5LIB=${PERL5LIB}" >> env.sh 
echo "export R_LIBS=${R_LIBS}" >> env.sh
echo "export cpanm_dir=${cpanm_dir}" >> env.sh
echo "export blast_dir=${blast_dir}" >> env.sh
echo "export rmblast_dir=${rmblast_dir}" >> env.sh
echo "export windowmasker_dir=${windowmasker_dir}" >> env.sh
# echo "export mummer4_dir=${mummer4_dir}" >> env.sh
echo "export mummer3_dir=${mummer3_dir}" >> env.sh
# echo "export assemblytics_dir=${assemblytics_dir}" >> env.sh
# echo "export gnuplot_dir=${gnuplot_dir}" >> env.sh
echo "export bedtools_dir=${bedtools_dir}" >> env.sh
echo "export sra_dir=${sra_dir}" >> env.sh
echo "export art_dir=${art_dir}" >> env.sh
echo "export trimmomatic_dir=${trimmomatic_dir}" >> env.sh
echo "export bwa_dir=${bwa_dir}" >> env.sh
echo "export samtools_dir=${samtools_dir}" >> env.sh
echo "export picard_dir=${picard_dir}" >> env.sh
echo "export gatk3_dir=${gatk3_dir}" >> env.sh
# echo "export gatk4_dir=${gatk4_dir}" >> env.sh
echo "export freebayes_dir=${freebayes_dir}" >> env.sh
echo "export gemtools_dir=${gemtools_dir}" >> env.sh
echo "export vcflib_dir=${vcflib_dir}" >> env.sh
# echo "export vcftools_dir=${vcftools_dir}" >> env.sh
echo "export vt_dir=${vt_dir}" >> env.sh
echo "export freec_dir=${freec_dir}" >> env.sh
echo "export parallel_dir=${parallel_dir}" >> env.sh
echo "export ucsc_dir=${ucsc_dir}" >> env.sh


# test java configuration: requireds java 1.8 
echo ""
echo "##########################################"
echo "Testing java configuration ..."
echo ""
java_bin=""
if type -p java
then 
    java_bin=$(which java)
    echo "found java executable in PATH: $java_bin"
elif [[ -n "$JAVA_HOME" ]] && [[ -x "$JAVA_HOME/bin/java" ]]
then 
    java_bin="$JAVA_HOME/bin/java"
    echo "found java executable in JAVA_HOME: $java_bin" 
else 
    echo "";
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
    echo "Failed to detect Java installation in the system!"
    echo "Please install java 1.8, which is a dependency of RecombineX!\n";
    echo "After the java installation, please manually set the directory path to java 1.8 executable on the last line of the env.sh file generated by this installation script!"
    echo "export java_dir=" >> env.sh
fi  

if [[ -n "$java_bin" ]]
then
    java_version=$("$java_bin" -version 2>&1 | awk -F '"' '/version/ {print $2}')
    echo "detected java_version: $java_version"
    if [ $(tidy_version "$java_version") -eq $(tidy_version "1.8") ]
    then
	java_dir=$(dirname $java_bin)
	echo "export java_dir=${java_dir}" >> env.sh
        echo "You have the correct java version for RecombineX! RecombineX will take care of the configuration."
    else
	echo "";
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
	echo "Your java version is not the version required by RecombineX (java v1.8)!"
        echo "Please manually set the directory path to java 1.8 executable on the last line of the env.sh file generated by this installation script!"
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
	echo "export java_dir=" >> env.sh
    fi
fi

echo ""
echo "[$(timestamp)] Uncompress large supporting files ..."
gunzip $RECOMBINEX_HOME/data/S288C.all_feature.gff.gz
gunzip $RECOMBINEX_HOME/data/S288C.genome.fa.gz
gunzip $RECOMBINEX_HOME/data/SK1.all_feature.gff.gz
gunzip $RECOMBINEX_HOME/data/SK1.genome.fa.gz

gunzip $RECOMBINEX_HOME/data/P1-P2.ref.final.SNP.markers.txt.gz
gunzip $RECOMBINEX_HOME/data/ref.centromere.relabel.gff.gz
gunzip $RECOMBINEX_HOME/data/ref.genome.raw.relabel.fa.gz

echo "[$(timestamp)] Done!"
echo ""
echo ""
echo "#################### IMPORTANT !!! #######################"
echo ""
echo "[$(timestamp)] Automatic dependencies installation finished! "
echo ""
echo "#########################################################"


############################
# checking Bash exit status

if [ $? -eq 0 ]
then
    echo ""
    echo "RecombineX message: This bash script has been successfully processed! :)"
    echo ""
    echo ""
    exit 0
fi
############################
