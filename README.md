### Description
Given a protein consisting of a single polypeptide chain available on the [Protein Data Bank](https://www.rcsb.org/), creates an H-tangle diagram of the protein and computes a polynomial invariant (single-strand Gamma) of this diagram.

### Usage
Run `main.py` with Python 3, and a Wolfram Language Script containing the data of the H-tangle diagram will be automatically generated in the same directory. This script may be evaluated with Mathematica to produce the Gamma polynomial of the protein.
Note: this program is only designed to work for proteins consisting of a single polypeptide chain.

### Configuration
The following variables in `main.py` may be configured:
- `PROTEIN_NAME`: the name of the protein as found on the Protein Data Bank (e.g., 1R02)
- `PROJECTION_PLANE_NORMAL`: a numpy array denoting the normal of the plane onto which the 3D structure of the protein is projected. In principle the choice of normal does not matter, however one should note that Gamma is only an invariant up to a factor of t^k in the polynomial, and only for diagrams of writhe zero (i.e., for diagrams with equal number of positive and negative crossings, which one can force by adding any number of positive or negative loops to each strand).
- `PLOT_GRAPH`: if `True`, plots all the residues in the chain with all relevant bonds and crossings.

### Dependencies
- numpy
- scipy
- [prody](http://prody.csb.pitt.edu/)
- matplotlib (only for `ProteinProjection.plotCrossings()`)

![protein-1R02-graph](https://user-images.githubusercontent.com/38050380/128295577-9e542f78-b4a2-4f04-a6c4-53b3bd1e6157.png)

