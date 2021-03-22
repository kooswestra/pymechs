
**PyMechs** is a mixed C++ and python library. Based on a C++ engine with python bindings made using [pybind11](https://github.com/pybind/pybind11). As well as featuring python extensions made to easily design, simulate and evolve actuated and controlled planar mechanisms. Based on a modified version of the graph representation of mechanisms by Kuppens, the paper is available [here](https://ieeexplore.ieee.org/document/8276271).

PyMechs uses the [Eigen library](https://eigen.tuxfamily.org/) as mathematical backend for vector and matrix functions. 

Written by Koos Westra as part of the Msc. Thesis: Automated Design Of Actuated Mechanisms

Doxygen documentation of the project is available in the docs folder after building.

# Features
___

**PyMechs** offers the following features

  * Fast mechanism simulation 
  * Mechanism structure animation and plotting
  * GUI based mechanism constructor

# Installation
___

**PyMechs** has the following dependencies

  * A C++11 or higher compatible compiler (GCC 4.8.1+, ...)
  * CMake 3.13 or higher
  * Python 3+ with numpy and matplotlib
  * [Eigen 3.3+](https://eigen.tuxfamily.org/)
  * [pybind11](https://github.com/pybind/pybind11)

Make sure these are satisfied, then clone the git repository:

    git clone https://gitlab/pymechlib
In the directory with the repository files run:

    python setup.py install --user
This compiles the C++ library with bindings and packages them with the python part of pure python modules to create a python package. Which is finally installed locally.
It is possible to test if the installation went correctly by running:

    python setup.py tests

# Basic Usage
___

Taking a look at the file examples/simplependulum.py shows how to implement the mechanism with this graph based representation using **PyMechs** to construct, simulate, plot and animate the mechanism. <br>

~~~~{.py}
""" 
python mechanism library example file, the generation, simulation and animation
of a basic passive spring based mechanism.
@Author Koos Westra
"""

# Import the mechanism library
import pymechlib as mch
# Get the edge labels names
Labels = mch.DNA.Labels
# Numpy is used for the matrix structures
import numpy as np

# Setting the simulation parameters ------------------------------

# The incidence matrix for a single pendulum is a single column, with a Hinge as edge
incidence_matrix = np.array([[1],
                             [1]])
                             
edge_labels = np.array([Labels.Hinge])

# The mass is defined by the first element, the second 2 define the position
masses = np.array([[1,0.2,-0.5]])

# Set the hinge parameters, only a hinge at the ground location
hinge_parameters = np.array([[0,0]])

# Generate a dna structure with these parameters
dna = mch.DNA(incidence_matrix, edge_labels, masses, hinge_parameters)

# And then generate an individual with this DNA
mech = mch.Mechanism(dna)

# Simulate the mechanism
time = 10
steps = 1000
mech.simulate(time, steps)

# Plot the state trajectories using the plot function
mech.plot()

# And finally show an animation of the mechanism
mech.animate()

~~~~


# About

___

This project was created by Koos Westra as part of the Msc. Thesis Automated Design of Actuated Mechanisms at the TU Delft.

## License

___

This project is provided under the LGPL3 license that can be found in the LICENSE file. By using, distributing, or contributing to this project, you agree to the terms and conditions of this license.
