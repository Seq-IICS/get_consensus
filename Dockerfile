FROM ubuntu:20.04
ENV CONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"

RUN apt-get update && apt-get install -y \
    bcftools \
    minimap2 \
    libbz2-dev \
    liblzma-dev \
    autotools-dev \
    python3-pip \
    wget \
    seqtk && apt-get clean
RUN pip install -v --retries 10000 NanoPlot && rm -rf /root/.cache/pip
workdir /tmp
RUN wget -c "$CONDA_URL" -O miniconda.sh -q && chmod +x miniconda.sh && bash miniconda.sh -b -p /opt/conda && rm miniconda.sh
ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
RUN conda config \
    --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -c bioconda htslib cramino ivar && \
    conda clean -a
