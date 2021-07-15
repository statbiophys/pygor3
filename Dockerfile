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
RUN rm -r /programs/IGoR

# WORKDIR /programs
# RUN git clone https://github.com/statbiophys/pygor3.git
# WORKDIR /programs/pygor3
# RUN apt-get install -y python-pip
# RUN pip install pygor3
RUN apt-get install -y python3.7

RUN useradd -ms /bin/bash ceor
RUN mkdir /igor_data
RUN chmod 777 /igor_data
RUN chown ceor:ceor /igor_data
RUN apt-get install -y wget
RUN apt-get install -y curl
RUN apt-get install -y python3-pip
RUN chown ceor:ceor /programs

RUN python3 --version
RUN python3 -m pip install --upgrade --force pip
RUN pip -V

USER ceor

VOLUME /igor_data
WORKDIR /igor_data

RUN ls
WORKDIR /programs
# RUN wget "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
# RUN chmod +x Miniconda3-latest-Linux-x86_64.sh 
# RUN ./Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda

WORKDIR /igor_data

ENV PYTHONIOENCODING=utf-8
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

ENV PATH=$PATH:/home/ceor/.local/bin
RUN echo $PATH
RUN pip install -U --pre pygor3
RUN pygor --version
# RUN pip install appdirs
RUN ls $HOME/.local/bin
RUN pip show pygor3
# CMD ["pygor"]
RUN pygor --version
RUN python3 --version

ENTRYPOINT ["echo"]
# ENTRYPOINT ["pygor"]
# ENTRYPOINT ["/home/ceor/.local/bin/pygor"]

#ENTRYPOINT ["/home/ceor/.local/lib/python3.6/site-packages/pygor"]

#ENTRYPOINT ["igor"]

#ENTRYPOINT ["${HOME}/miniconda/bin/conda", "run", "--no-capture-output", "-n", "statbiophys", "pygor"]
#ENTRYPOINT ["${HOME}/miniconda/envs/statbiophys/bin/pygor"]
#ENTRYPOINT ["pygor"]

## #RUN /root/miniconda/bin/conda init
## RUN ${HOME}/miniconda/bin/conda create --name statbiophys python=3.7
## #RUN /root/miniconda/bin/conda activate statbiophys
## RUN ls ${HOME}/miniconda/envs/statbiophys/
## RUN ${HOME}/miniconda/envs/statbiophys/bin/pip install --pre pygor3
## #RUN /root/conda/envs/statbiophys/bin/pip install pygor3
## RUN python --version
## RUN ${HOME}/miniconda/envs/statbiophys/bin/pygor --version
## RUN ${HOME}/miniconda/envs/statbiophys/bin/pygor --help
## RUN igor -version
