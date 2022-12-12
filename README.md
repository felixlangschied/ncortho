# ncOrtho
[![PyPI version](https://badge.fury.io/py/ncOrtho.png)](https://pypi.org/project/ncOrtho/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

NcOrtho is a tool for the targeted search of orthologous micro RNAs (miRNAs) throughout the tree of life. 
Conceptually, it works similar to the program [fDOG](https://github.com/BIONF/fDOG) in that a probabilistic model of 
a reference miRNA is created. For training the model, orthologs of the reference sequence are first identified in 
a set of taxa that are more closely related to the reference species. In contrast to fDOG, ncOrtho does not train hidden 
Markov Models but covariance models (CMs) (Eddy & Durbin, 1994) 
which also model conservation of the miRNA's secondary structure.

![workflow](https://github.com/felixlangschied/ncortho/blob/master/ncOrtho/docs/figure1_ncortho_worklfow.png)

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

After installing all three dependencies, ncOrtho can be installed with `pip`:
```
 pip install ncOrtho
```

## Usage
### CM construction
As a targeted search for orthologs, ncOrtho's biggest strength is its flexibility to change the taxonomic 
scope of an analysis according to the research question at hand.

For this reason, a few questions need to be answered, before we can start constructing covariance models:
1. What is the reference species?
2. How phylogenetically diverse will my set of target species be?
3. From which species should the core set of miRNA orthologs be extracted, which will be used for training the CMs?
4. Which miRNAs are going to be used for the ortholog search?

You can (soon) find more information on how to answer these questions in the 
[WIKI](https://github.com/felixlangschied/ncortho/wiki/Creating-miRNA-covariance-models#choosing-core-species).
As soon as you know what your core species are going to be, you will need to collect the following data:

* Genomic sequence in FASTA format (e.g "genomic.fna" from RefSeq)
* Genome annotation in GFF3 format (e.g. "genomic.gff" from RefSeq)
* Pairwise orthologs of all proteins between the reference and each core species (more information 
[here](https://github.com/felixlangschied/ncortho/wiki/Input-Data#pairwise-orthologs)

Modify the [example parameters](ncOrtho/coreset/example_parameters.yaml) file to contain all 
relevant paths to your input files. The "name" property of your reference and core species has to merely be a 
unique identifier. It is however recommended to use the correct species names to increase readability.

Additional to the parameters file, you will need a tab separated file containing the position and sequence of each 
miRNA for which a model should be constructed (more information 
[here](https://github.com/felixlangschied/ncortho/wiki/Input-Data#reference-mirnas)). 

You can then start CM construction with:
```
ncCreate -p <parameters.yaml> -n <mirnas.tsv> -o <outdir>
```
If you encounter errors, make sure that:
* The identifiers in the pairwise orthologs files match the ones in the gff files (use the `-idtype=` flag to use other
ID types)
* The contig/chromosome column in tab separated miRNA input file match the contig/chromosome id 
in the reference gff file

Use `ncCreate -h` to see all available options for CM construction.

### Orthology search

You can start the orthology search with:
```
ncSearch -m <CMs/> -n <mirnas.tsv> -q <query_genome.fa> -r <reference_genome.fa> -o <outdir>
```

Use `ncSearch -h` to see all available options for the orthology search or have a look at the 
[WIKI](https://github.com/felixlangschied/ncortho/wiki/Running-the-orthology-search).

## Support

Please refer to our Wiki Page of [known issues](https://github.com/felixlangschied/ncortho/wiki/Known-Issues) first, 
then consider opening an issue on GitHub or contacting me directly via [mail](langschied@bio.uni-frankfurt.de)

## Contributors

* [Felix Langschied](https://github.com/felixlangschied)
* [Andreas Blaumeiser](https://github.com/acblaumeiser)
* Mirko Brüggemann
* Daniel Amsel

Dept. for Applied Bioinformatics Institute for Cell Biology and Neurosciences, Goethe University, Frankfurt am Main

* [Ingo Ebersberger](https://www.bio.uni-frankfurt.de/43045195/Abt__Ebersberger___Biowissenschaften)
 
## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* [Lorenz et al. 2011](https://almob.biomedcentral.com/articles/10.1186/1748-7188-6-26): ViennaRNA Package 2.0
* [Nawrocki et al. 2013](https://academic.oup.com/bioinformatics/article/29/22/2933/316439): Infernal 1.1: 
100-fold faster RNA homology searches
* [Notredame et al. 2000](http://www.tcoffee.org/Publications/Pdf/tcoffee.pdf): T-Coffee: A Novel Method for Fast and 
AccurateMultiple Sequence Alignment
* [Shirley et al. 2015](https://peerj.com/preprints/970v1/): Efficient "pythonic" access to FASTA files using pyfaidx

## Contact

For support or bug reports please contact: [langschied@bio.uni-frankfurt.de](langschied@bio.uni-frankfurt.de)