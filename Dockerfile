FROM ubuntu:focal

RUN groupadd -r myuser && useradd -m -r -g myuser myuser

RUN apt update -y && DEBIAN_FRONTEND=noninteractive apt-get -y install tzdata && apt install -y cpanminus libpng-dev libssl-dev python3.9 \
    python3.9-distutils python3.9-dev vim git sqlite3 libdbi-perl curl wget gzip libarchive-extract-perl zip build-essential libbz2-dev \
    liblzma-dev zlib1g-dev libarchive-zip-perl libmodule-build-perl samtools libcurl4-openssl-dev tabix libxml2-dev

ENV PATH="/home/myuser/perl5/bin:${PATH}"
ENV PERL5LIB="/home/myuser/perl5/lib/perl5:${PERL5LIB}"
ENV PERL_LOCAL_LIB_ROOT="/home/myuser/perl5:${PERL_LOCAL_LIB_ROOT}"
ENV PERL_MB_OPT="--install_base \"/home/myuser/perl5\""
ENV PERL_MM_OPT="INSTALL_BASE=/home/myuser/perl5"
ENV KENT_SRC=/home/myuser/kent-335_base/src
ENV CFLAGS="-fPIC"
USER myuser
WORKDIR /home/myuser

RUN mkdir -p /home/myuser/cpanm && mkdir -p /home/myuser/kent-335_base/src

RUN wget https://github.com/ucscGenomeBrowser/kent/archive/v335_base.tar.gz && \
    tar xzf v335_base.tar.gz && \
    export MYSQLINC=`mysql_config --include | sed -e 's/^-I//g'` && \
    export MYSQLLIBS=`mysql_config --libs` && \
    export MACHTYPE=$(uname -m) && cd $KENT_SRC/lib && echo 'CFLAGS="-fPIC"' > ../inc/localEnvironment.mk && \
    make clean && make && cd ../jkOwnLib && make clean && make && \
    ln -s $KENT_SRC/lib/x86_64/* $KENT_SRC/lib/

#RUN cd /home/myuser && perl -MCPAN -e'install "Bio::DB::BigFile"' &&
RUN perl -MCPAN -e'install "LWP::Simple"'

RUN cd /home/myuser && git clone https://github.com/Ensembl/ensembl-vep.git && \
    cd ensembl-vep && git pull && git checkout release/108 && \
    perl INSTALL.pl -a a -n --SPECIES homo_sapiens --ASSEMBLY GRCh38

RUN cd /home/myuser && mkdir /home/myuser/bcftools && git clone --recurse-submodules https://github.com/samtools/htslib.git && \
    git clone https://github.com/samtools/bcftools.git && \
    cd bcftools && git checkout tags/1.16 && make prefix=/home/myuser/bcftools && make prefix=/home/myuser/bcftools install

RUN cd /home/myuser && curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && python3.9 get-pip.py && export PATH="$PATH:/home/myuser/.local/bin" && mkdir -p /home/myuser/pypack/ && pip install --no-cache-dir --target=/home/myuser/pypack fastparquet pandas numpy tqdm pytabix genomicsqlite && \
    mkdir -p /home/myuser/bin && ln -s /usr/bin/python3.9 /home/myuser/bin/python && ln -s /home/myuser/.local/bin/pip /home/myuser/bin/pip

RUN perl -MCPAN -e'install "Bio::EnsEMBL::XS"' && perl -MCPAN -e'install "Set::IntervalTree"' && cpanm -l $HOME/cpanm XML::LibXML && cpanm -l $HOME/cpanm Bio::Root::Version && cpanm -l $HOME/cpanm Try::Tiny && cpanm -l $HOME/cpanm Bio::DB::BigFile

ENV BCFTOOLS_PLUGINS="/home/myuser/bcftools/plugins"
ENV PATH="/home/myuser/.local/bin/:/home/myuser/bin:/home/myuser/pypack/bin:/home/myuser/ensembl-vep:/home/myuser/bcftools:$PATH"
ENV PYTHONPATH="/home/myuser/pypack:$PYTHONPATH"
ENV PERL5LIB="/home/myuser/cpanm/lib/perl5/x86_64-linux-gnu-thread-multi:/home/myuser/cpanm/lib/perl5:${PERL5LIB}"
WORKDIR /home/myuser/work

ENTRYPOINT ["/bin/bash", "/home/myuser/work/scripts/annotate.sh"]
