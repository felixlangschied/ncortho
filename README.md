# ncortho
NcOrtho is a tool for the targeted search of orthologous micro RNAs (miRNAs) throughout the tree of life. Conceptually, it works similar to the program [fDOG](https://github.com/BIONF/fDOG) in that a probabilistic model of the target sequence in a reference species is created. For this, the orthologs of the target sequence are identified in a set of taxa that are closely related to the reference species. In contrast to fDOG, ncOrtho does not train hidden Markov Models but covariance models (CMs) (Eddy & Durbin, 1994) which also model conserved secondary structures of the miRNAs

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them

```
tools = ['blastn', 'infernal', 't_coffee', 'alifold']
```

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo


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

