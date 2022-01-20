
# Preparation for 2022-01-25 tutorial

## Install IGOR
### Gnu-Linux
```console
$ sudo apt-get install build-essential
$ git clone https://github.com/statbiophys/IGoR.git
$ cd IGoR
$ ./configure && make && sudo make install
```

If you use another distribution than Ubuntu replace first command by
what's needed to install gcc.

If don't have root access replace last command by
```console
$ ./configure --prefix=${HOME}/.local/ && make && make install
``` 

### MacOS
```console
$ /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
$ brew install gcc@6 
$ git clone https://github.com/statbiophys/IGoR.git
$ cd IGoR
$ ./configure CC=gcc-6 CXX=g++-6 && make && sudo make install
``` 
If don't have root access replace last command by
```console
$ ./configure CC=gcc-6 CXX=g++-6 --prefix=${HOME}/.local/ && make && make install
``` 

If you already have homebrew skip the first command. If you prefer to
use macports see [here](https://statbiophys.github.io/IGoR/#macos).

## Install pygor
First install [Anaconda](https://docs.anaconda.com/anaconda/install/).

Then run
```console
$ conda create --name statbiophys python=3.7
$ conda activate statbiophys
(statbiophys) $ conda install jupyterlab
(statbiophys) $ pip install pygor3 
(statbiophys) $ pygor demo-get-data
(statbiophys) $ cd demo
(statbiophys) $ jupyterlab
```


