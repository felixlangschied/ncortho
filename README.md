# ncortho
NcOrtho is a tool for the targeted search of orthologous micro RNAs (miRNAs) throughout the tree of life. 
Conceptually, it works similar to the program [fDOG](https://github.com/BIONF/fDOG) in that a probabilistic model of 
the target sequence in a reference species is created. For this, the orthologs of the target sequence are identified in 
a set of taxa that are closely related to the reference species. In contrast to fDOG, ncOrtho does not train hidden 
Markov Models but covariance models (CMs) (Eddy & Durbin, 1994) which also model conserved secondary structures of the 
miRNAs

## Getting Started

NcOrtho depends on multiple third party applications, some of which are Linux specific.
All dependencies can be installed with [Anaconda](https://www.anaconda.com/).
It is recommended to create a new Anaconda environment for this. For example:
```
conda create --name ncOrtho python=3.8
conda activate ncOrtho
```


### Prerequisites

* **Operating System:** Linux (tested on: Ubuntu 20.04)
* **Python:** version 3 or higher (tested with v3.8)

Tool | Tested version | Anaconda installation
------------ | ------------- | -------------
BLASTn | v2.7.1 | `conda install -c kantorlab blastn`
Infernal | v1.1.4 | `conda install -c bioconda infernal`
t_coffee | v13.45 | `conda install -c bioconda t-coffee`
alifold | v2.4.18 | `conda install -c bioconda viennarna`


### Installing

Currently, ncOrtho is not yet available as a python package. Until then, two python packages have to be installed
manually:
```
pip install pyfaidx
pip install biopython
```
Cloning the git repository is the easiest option to access ncOrtho at the moment:
```
git clone https://github.com/felixlangschied/ncortho.git
```

## Usage




## Contributors

Dept. for Applied Bioinformatics Institute for Cell Biology and Neurosciences, Goethe University, Frankfurt am Main

* [Ingo Ebersberger](https://github.com/ebersber)
* [Felix Langschied](https://github.com/felixlangschied)
* [Andreas Blaumeiser](https://github.com/acblaumeiser)
* Mirko Brüggemann
* Daniel Amsel


## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

