# innovation-lab

## Introduction

## Dependencies
- An instance of [anaconda](https://www.anaconda.com) or [miniconda](https://conda.io/en/latest/miniconda.html)
- IMB's Platform Load Sharing Facility (LSF) for resource management

## Installation
```
wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O $HOME/miniconda.sh
bash $HOME/miniconda.sh -b -p $HOME/.miniconda
export PATH="$HOME/.miniconda/bin:$PATH"
hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q conda
conda info -a
conda config --add channels r
conda config --add channels bioconda
mkdir -p $HOME/share $HOME/share/usr $HOME/share/usr/env
export INNOVATION_ENV_VERSION=0.0.1
conda create -q -p $HOME/share/usr/env/innovation-lab-$INNOVATION_ENV_VERSION
source activate $HOME/share/usr/env/innovation-lab-$INNOVATION_ENV_VERSION
git clone https://github.com/ndbrown6/innovation-lab.git
```

