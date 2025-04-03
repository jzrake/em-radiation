# Radiation Module

A Python extension module for calculating and visualizing electromagnetic radiation from moving charged particles.

## Overview

This package provides computational tools for simulating and visualizing the electromagnetic fields and potentials associated with the motion of charged particles. It includes implementations of:

- Particle trajectory calculations for oscillating charges
- Liénard-Wiechert potentials computation
- Retarded time calculations for radiation fields
- Visualization utilities for electromagnetic fields

The core numerical routines are implemented in C++ for performance, with Python bindings provided via pybind11.

## Installation

You can install the package directly from source:

```bash
pip install .
```

### Requirements

- Python ≥ 3.8
- C++ compiler with C++23 support
- CMake ≥ 3.18
- NumPy
- Matplotlib

## Usage

```python
from radiation import ParticleTrajectory, Vector3D, find_retarded_time

# Create a particle trajectory with specified angular frequency
particle = ParticleTrajectory(omega=0.9)

# Define an observation point
r = Vector3D(0.0, 0.0, 2.0)
t = 0.7

# Find the retarded time for radiation from the particle
tr = find_retarded_time(particle, r, t)

# Calculate Liénard-Wiechert potentials
scalar_potential, vector_potential = lienard_wiechert_potentials(particle, r, t)
```

## Examples

The package includes example visualizations, such as:

```python
from matplotlib.pyplot import figure, show
from radiation import ParticleTrajectory, Vector3D

# Plot the retarded time condition
particle = ParticleTrajectory(omega=0.9)
r = Vector3D(0.0, 0.0, 2.0)
t = 0.7
# ... visualization code ...
```

See `main.py` for more detailed examples of plotting retarded time conditions and electromagnetic field visualizations.

## Development

This project uses CMake for C++ compilation and scikit-build-core for Python packaging. The C++ code is bound to Python using pybind11.

To set up a development environment:

```bash
# Create and activate virtual environment
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install development dependencies
pip install -e .
```

## License

[License information]

## References

- Jackson, J. D. (1999). Classical Electrodynamics (3rd ed.). Wiley.
- [Additional references on electromagnetic radiation]
