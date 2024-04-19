FROM ubuntu:focal

RUN apt update -y && DEBIAN_FRONTEND=noninteractive apt-get -y install tzdata && apt install -y cpanminus libpng-dev libssl-dev python3.9 \
    python3.9-distutils python3.9-dev vim git sqlite3 libdbi-perl curl wget gzip libarchive-extract-perl zip build-essential libbz2-dev \ 
    liblzma-dev zlib1g-dev libarchive-zip-perl libmodule-build-perl && mkdir -p $HOME/cpanm
ENV PERL5LIB $HOME/cpanm/lib/perl5
RUN cpanm -l $HOME/cpanm Set::IntervalTree && cd / && wget https://github.com/ucscGenomeBrowser/kent/archive/v335_base.tar.gz && \
    tar xzf v335_base.tar.gz
ENV KENT_SRC=/kent-335_base/src
ENV CFLAGS="-fPIC"
RUN export MYSQLINC=`mysql_config --include | sed -e 's/^-I//g'` && \
    export MYSQLLIBS=`mysql_config --libs` && \
    export MACHTYPE=$(uname -m) && cd $KENT_SRC/lib && echo 'CFLAGS="-fPIC"' > ../inc/localEnvironment.mk && \
    make clean && make && cd ../jkOwnLib && make clean && make && \
    ln -s $KENT_SRC/lib/x86_64/* $KENT_SRC/lib/ && perl -MCPAN -e'install "Bio::DB::BigFile"' && perl -MCPAN -e'install "LWP::Simple"' && \
    cd / && git clone https://github.com/Ensembl/ensembl-vep.git && \ 
    cd ensembl-vep && git pull && git checkout release/108 && \
    perl INSTALL.pl -a a -n --SPECIES homo_sapiens --ASSEMBLY GRCh38 && \
    apt install -y libcurl4-openssl-dev tabix && cd && git clone --recurse-submodules https://github.com/samtools/htslib.git && \
    git clone https://github.com/samtools/bcftools.git && \
    cd bcftools && git checkout tags/1.16 && make && curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && python3.9 get-pip.py && pip install fastparquet pandas numpy tqdm pytabix genomicsqlite && \
    ln -s /usr/bin/python3.9 /usr/bin/python
RUN perl -MCPAN -e'install "Bio::EnsEMBL::XS"' && perl -MCPAN -e'install "Set::IntervalTree"'

ENV BCFTOOLS_PLUGINS="/root/bcftools/plugins"
ENV PATH="$PATH:/ensembl-vep:/root/bcftools"
WORKDIR /home

ENTRYPOINT ["/bin/bash", "/home/scripts/annotate.sh"]
