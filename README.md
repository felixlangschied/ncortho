# ncOrtho
NcOrtho is a tool for the targeted search of orthologous micro RNAs (miRNAs) throughout the tree of life. 
Conceptually, it works similar to the program [fDOG](https://github.com/BIONF/fDOG) in that a probabilistic model of 
a reference miRNA is created. For training the model, orthologs of the reference sequence are first identified in 
a set of taxa that are more closely related to the reference species. In contrast to fDOG, ncOrtho does not train hidden 
Markov Models but covariance models (CMs) (Eddy & Durbin, 1994) which also model the secondary structure of miRNAs.

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

### Installing
Currently, ncOrtho is not yet available as a python package. Until then, some python packages have to be installed
manually:
```
pip install pyfaidx
pip install biopython
pip install pyyaml
```
Cloning the git repository is the easiest option to access ncOrtho at the moment:
```
git clone https://github.com/felixlangschied/ncortho.git
```

## Usage
### CM construction
As a targeted search for orthologs, ncOrtho's biggest strength is its flexibility to change the taxonomic 
scope of an analysis according to the research question at hand.

For this reason, a few questions need to be answered, before we can start constructing covariance models:
1. What is the reference species?
2. How phylogenetically diverse will my set of target species be?
3. From which species should the core set of miRNA orthologs for training the CMs be extracted?
4. Which miRNAs should be used for the ortholog search?

You can (soon) find more information on how to answer these questions in the 
[WIKI](https://github.com/felixlangschied/ncortho/wiki/Creating-miRNA-covariance-models#choosing-core-species).







```
python nc_coreset.py -p <parameters.yaml> -n <mirnas.tsv> -o <outdir>
```








## Known Issues
#### MAX_N_PID
Example:
```
--ERROR: MAX_N_PID exceded -- Recompile changing the value of MAX_N_PID (current: 260000 Requested: 2346303)
```
This T-Coffee error has been encountered when running ncOrtho on a SLURM scheduler. It can be fixed by installing the
beta version of T-Coffee, as explained in this [discussion](https://groups.google.com/g/tcoffee/c/sO8Kd5NjA5A).
With the beta version of T-Coffee, you can increase the MAX_N_PID value according to your needs.
Include a similar line to the script calling ncOrtho for this:
```
export MAX_N_PID_4_TCOFFEE=4194304
```

## Contributors

* [Felix Langschied](https://github.com/felixlangschied)
* [Andreas Blaumeiser](https://github.com/acblaumeiser)
* Mirko Brüggemann
* Daniel Amsel

Dept. for Applied Bioinformatics Institute for Cell Biology and Neurosciences, Goethe University, Frankfurt am Main

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* [Nawrocki et al. 2013](https://academic.oup.com/bioinformatics/article/29/22/2933/316439): Infernal 1.1: 
100-fold faster RNA homology searches
* [Notredame et al. 2000](http://www.tcoffee.org/Publications/Pdf/tcoffee.pdf): T-Coffee: A Novel Method for Fast and 
AccurateMultiple Sequence Alignment
* [Shirley et al. 2015](https://peerj.com/preprints/970v1/): Efficient "pythonic" access to FASTA files using pyfaidx

