
# `origamiUROP`

**Autonomous generation of DNA Origami nanostructures**  
[![Build Status](https://travis-ci.com/shanilpanara/origamiUROP.svg?branch=master)](https://travis-ci.com/shanilpanara/origamiUROP)
[![GitHub Last Commit](https://img.shields.io/github/last-commit/shanilpanara/origamiUROP?color=blue)](https://github.com/shanilpanara/origamiUROP/commits/master)
[![GitHub Issues](https://img.shields.io/github/issues/shanilpanara/origamiUROP?color=red)](https://github.com/shanilpanara/origamiUROP/issues)
[![GitHub Pull Requests](https://img.shields.io/github/issues-pr/shanilpanara/origamiUROP?color=purple)](https://github.com/shanilpanara/origamiUROP/pulls)
![GitHub License](https://img.shields.io/github/license/shanilpanara/origamiUROP?color=yellow)

`origamiUROP` is a Python Library which can currently be used for generating a system comprised of ssDNA/dsDNA strands. Using a set of vertices, lines (edges) are generated which act as centrelines for single/double strands of DNA. It is modelled using the oxDNA model and can be exported in this format (`.conf` and `.top` files).

The end goal of the package is to **autonomously** generate many possible DNA origami scaffold configurations and subsequently their complimentary staple strands for a given 2D shape.



## 📃 Features

- Generate a polygon, given its verticies, where every edge is individually accessible
- Write the polygon to `.STL` and `.PLY` file formats
- Create an oxDNA system with tools to control the length, sequence, base-pairs/nucleotides per 2π turn and more.
- Write the oxDNA system to `.conf` and `.top` formats

### 📆 Upcoming

- Write the oxDNA system to a `LAMMPS` configuration file
- Read in oxDNA `.conf` and `.top` files
- Read in LAMMPS `dump` files

## 🔋 Installation

At the time of writing it is only possible to manually install the package. This requires the following prerequisites:

- numpy
- pandas
- matplotlib
- meshio

To install, run the following commands:
```bash
    git clone https://github.com/shanilpanara/origamiUROP.git
    cd origamiUROP.git
    pip install .
```
And to test:
```bash
    pytest
```
## ✍️ Example

We can create an `Edge` object which will act as the centreline ssDNA/dsDNA. Then by creating and adding DNA strands to a system, it allows us to write out a file in an oxDNA format.

```python
from origamiUROP.oxdna import Strand, System
from origamiUROP.polygons import Edge
import numpy as np

# Create two `Edge` objects, with a length of 6.5 and 14 units, respectively
edge_1 = Edge(np.array([-5.5, 0.0, 0.0]), np.array([1.0, 0.0, 0]))
edge_2 = Edge(np.array([1.0, 0.0, 0.0]), np.array([15.0, 0.0, 0.0]))

# Create a single stranded DNA object with 6 nucleotides (rounded down from length of Edge)
# a random sequence will be assigned
strand_ssDNA = edge_1.strand()

# Create a double stranded DNA object with 14 base pairs and assign a sequence for one
# of the strands (a complementary sequence will be assigned to the second strand)
strand_dsDNA = edge_2.strand(double=True, sequence="CCAAGGTTCAGTCA")
print("ssDNA: ", strand_ssDNA[0], "\ndsDNA: ", strand_dsDNA[0], strand_dsDNA[1])

# Create system, add strands and write oxDNA files
system = System(np.array([20.0, 20.0, 20.0]))
strands_to_add = strand_ssDNA + strand_dsDNA
system.add_strands(strands_to_add)

# Creates 2 files -> oxdna.out.top & oxdna.out.conf
system.write_oxDNA()

# View pandas dataframe
system.dataframe
```
![Ovito Visualisation](Example_OvitoVisualisation.png)  
*`system` can be visualised using Ovito*

## 🎉 Acknowledgements

- oxDNA system inspired by code from https://github.com/lorenzo-rovigatti/oxDNA/
