# TO CREATE IMAGE docker build -t dpygor3:dev .
FROM ubuntu:bionic

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
	apt-get install -y build-essential git autotools-dev
 
# download doxygen and asciidoctor
RUN apt-get install -y texlive-base texlive-latex-recommended texlive-latex-extra
RUN apt-get install -y asciidoctor doxygen
RUN apt-get install -y texlive-font-utils
RUN apt-get install -y autoconf 
RUN apt-get install -y libtool 

WORKDIR /programs

RUN git clone https://github.com/statbiophys/IGoR.git 
WORKDIR /programs/IGoR 

RUN ./autogen.sh
RUN ./configure && make && make install

WORKDIR /programs
RUN apt-get install -y wget
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda

#RUN /root/miniconda/bin/conda init
RUN /root/miniconda/bin/conda create --name statbiophys python=3.7
#RUN /root/miniconda/bin/conda activate statbiophys
RUN  ls /root/miniconda/envs/statbiophys/
RUN  /root/miniconda/envs/statbiophys/bin/pip install pygor3
#RUN /root/conda/envs/statbiophys/bin/pip install pygor3

# WORKDIR /programs
# RUN git clone https://github.com/statbiophys/pygor3.git
# WORKDIR /programs/pygor3
# RUN apt-get install -y python-pip

#ENTRYPOINT ["igor"]

ENTRYPOINT ["/root/miniconda/bin/conda", "run", "--no-capture-output", "-n", "statbiophys", "pygor"]
