#FROM nvidia/cuda:11.6.1-devel-ubuntu20.04
FROM nvidia/cuda:11.8.0-devel-ubuntu20.04
#FROM nvidia/cuda:11.0.3-devel-ubuntu18.04
#FROM nvidia/cuda:12.2.2-base-ubuntu20.04 # con questo non builda

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y apt-utils tzdata && \
    apt-get upgrade -y && \
    apt-get install -y csh flex bison patch gcc gfortran g++ make \
                       build-essential xorg-dev xutils-dev libbz2-dev zlib1g-dev \
                       libboost-dev libboost-thread-dev libboost-system-dev mpich libmpich-dev \
                       python3 python-dev wget bc cmake git gromacs doxygen libfftw3-dev automake lsb-core

ENV PATH="${PATH}:/opt/miniconda/bin:${PATH}"

USER root

COPY ./Amber24.tar.bz2 /opt/Amber24.tar.bz2
COPY ./AmberTools24.tar.bz2 /opt/AmberTools24.tar.bz2

WORKDIR /opt/

RUN tar xvfj Amber24.tar.bz2 && \
    tar xvfj AmberTools24.tar.bz2 && \
    mkdir -p /opt/amber24_src/build/ && \
    mkdir -p /work

WORKDIR /opt/amber24_src/build/

COPY run_cmake /opt/amber24_src/build/

RUN pwd && \
    ls -altr && \
    cat /opt/amber24_src/build/run_cmake && \
    bash "/opt/amber24_src/build/run_cmake" && \
    make install

WORKDIR /work

ENV AMBERHOME=/opt/amber24
ENV PATH=/opt/amber24/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/amber24/lib:$LD_LIBRARY_PATH
