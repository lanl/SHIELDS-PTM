[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4891973.svg)](https://doi.org/10.5281/zenodo.4891973)

# PTM - The SHIELDS Particle Tracing Model
The SHIELDS Particle Tracing Model (PTM) is a large-scale particle tracing code, primarily intended for magnetospheric
applications. It provides a fast, parallel computing-friendly, and versatile framework for tracing charged particles in arbitrary
electric and magnetic fields. To support a broad range of test particle applications, full Lorentz motion can be modeled as
well as guiding center motion. Switching strategies allow accurate and efficient test particle tracing across regions where
the guiding center approximation may not be maintained. PTM supports tracing of particle motion both forwards and
backwards in time. For backwards tracing, PTM can launch particles from a specified location and trace until stopping
criteria are met. These stopping criteria can include crossing a surface in space, or a time limit. For forwards tracing
particles can be launched from specified points or spatial regions. The code is designed to allow many applications without
requiring recompilation, instead providing a rich, user-configurable interface that allows use of different modes at runtime.
Model results are stored in user-selected format for post-processing including Liouville mapping of particle fluxes from
source to target location. PTM includes pre-processing for field configurations, tools for setting up runs, and postprocessing
capabilities to enable scientific use. Sample applications include specifying the geosynchronous particle
boundary for inner magnetosphere modeling, and determining access of solar energetic protons.

This code is approved for open source release by Los Alamos National Laboratory and has been assigned reference number C21002.

### Attribution
The following people have contributed to the edvelopment of SHIELDS-PTM:
Jesse Woodroffe [1], Steven Morley [1], Thiago Brito [1], Alin-Daniel Panaitescu [2], Miles Engel [1], Michael Henderson [1], Vania Jordanova [1]

1. ISR-1, Space Science & Applications, Los Alamos National Laboratory
2. ISR-2, Space and Remote Sensing, Los Alamos National Laboratory

When using this code, please acknowledge its use, refer the reader to the GitHub repository and cite the master DOI at Zenodo.
A paper describing the code is forthcoming, but until this is updated the best references for this code are:
- [Woodroffe et al., Data-optimized source modeling with the Backwards Liouville Test-Kinetic method, JASTP, 2018](https://doi.org/10.1016/j.jastp.2017.09.010)
- [Jordanova et al., Specification of the near-Earth space environment with SHIELDS, JASTP, 2018](https://doi.org/10.1016/j.jastp.2017.11.006)
- [Brito et al., Particle tracing modeling of ion fluxes at geosynchronous orbit, JASTP, 2018](https://doi.org/10.1016/j.jastp.2017.10.008)

### About the code
This repository contains the simulation code and analysis tools for the SHIELDS particle tracing model (PTM).

The main source is written in Fortran2008. Associated scripts, pre-processing, and post-processing are written in Python 3.

To build PTM, simply run:
```
make
```

To build the fortran only run `make ptm`, to build/install the Python module only run `make python`, and to convert the Markdown documents to PDF use `make docs`.
The documentation conversion uses the Node.js `mdpdf` module, which can be installed using `npm install mdpdf -g`. PDF documentation is placed in the `docs` directory.

### Use of RKSUITE
In addition to the basic RK4 integrator, SHIELDS-PTM provides an interface to the RKSUITE library allowing users to select higher-order, adaptive integrators.
If these integrators are used, RKSUITE and its originators should be credited. For details, please see [rksuite_readme](src/rksuite_readme).
