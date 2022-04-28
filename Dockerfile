FROM nvidia/cuda:11.6.0-devel-ubuntu20.04

RUN apt update -y &&\
    apt upgrade -y
# generic softwares needed before installing python, julia or R
RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime &&\
    apt install -y tzdata &&\
    apt install -y curl &&\
    apt install -y software-properties-common &&\
    add-apt-repository ppa:deadsnakes/ppa &&\
    # gdebi asks for timezone info so preset it to avoid docker breaking on an interactive piece
    # apt install -y gdebi &&\
    apt install -y wget &&\
    apt install -y tar

ENV JUL_VERS=1.7.2
ENV JUL_MAJ_MIN_VERS=1.7
ENV JULIA_NUM_THREADS=12
ENV PY_VERS=3.10
# install Julia
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/${JUL_MAJ_MIN_VERS}/julia-${JUL_VERS}-linux-x86_64.tar.gz &&\
    tar zxvf julia-${JUL_VERS}-linux-x86_64.tar.gz &&\
    ln -s /julia-${JUL_VERS}/bin/julia usr/local/bin/julia
# install Python
RUN apt update -y &&\
    apt upgrade -y &&\
    apt install -y python${PY_VERS}-distutils python${PY_VERS}-dev &&\
    # libraries needed for python dependencies
    apt install -y libcairo2-dev libjpeg-dev libgif-dev &&\ 
    apt install -y python${PY_VERS} python3-pip &&\
    curl -sS https://bootstrap.pypa.io/get-pip.py | python${PY_VERS}
# make any julia instance activate a local environment on entry
RUN echo 'import Pkg; Pkg.activate("InvestmentManagement")' >> ./julia-${JUL_VERS}/etc/julia/startup.jl

WORKDIR app/
COPY InvestmentManagement/* InvestmentManagement/
ENV ALPHA_VANTAGE_API_KEY="INSERT_KEY_HERE"
ENV PYTHON = /usr/bin/python3.10
RUN julia -e 'Pkg.instantiate()'
RUN julia -e 'Pkg.precompile()'