FROM ubuntu:22.04 AS base

LABEL maintainer=ting-hsuan.chen@plantandfood.co.nz

LABEL description="TEgenomeSimulator, a tool to simulate TE mutation and insertion into a random or user-provided genome."

LABEL version=1.0.0

ENV PATH="/opt/RepeatMasker:/usr/local/bin:$PATH"

WORKDIR /opt

RUN apt-get update \
    && apt-get install -y \ 
    build-essential \
    wget \
    curl \
    perl \
    csh \
    zlib1g-dev \
    python3 \
    python3-pip \
    samtools \
    unzip \
    git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /tmp/downloaded_packages

RUN pip3 install numpy==1.26.2 pandas==2.2.1 \
    biopython==1.81 pyyaml==6.0.1 h5py==3.14.0 matplotlib==3.10.3 

FROM base AS builder

RUN cd /tmp \
    &&  wget https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.8.tar.gz && \
    tar -xzf RepeatMasker-4.1.8.tar.gz && \
    cd RepeatMasker && \
    mkdir -p Libraries && \
    mv /tmp/RepeatMasker /opt

RUN cd /tmp \
    && wget https://www.repeatmasker.org/rmblast/rmblast-2.14.1+-x64-linux.tar.gz && \
    tar -xzf rmblast-2.14.1+-x64-linux.tar.gz && \
    mv rmblast-2.14.1 /opt/RMBlast

RUN cd /tmp \
    && wget https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64 -O trf && \
    chmod +x trf && \
    mv trf /opt/RepeatMasker/trf

FROM base AS runtime

COPY --from=builder /opt /opt

RUN cd /opt/RepeatMasker \
    && echo -e "\n\n" | ./configure \
    -libdir /opt/RepeatMasker/Libraries \
    -trf_prgm /opt/RepeatMasker/trf \
    -rmblast_dir /opt/RepeatMasker/RMBlast/bin \
    -default_search_engine rmblast

RUN cd /opt && \  
    git clone --depth 1 --branch v1.0.0 \
    https://github.com/Plant-Food-Research-Open/TEgenomeSimulator.git && \
    cd /opt/TEgenomeSimulator && \
    pip3 install .
