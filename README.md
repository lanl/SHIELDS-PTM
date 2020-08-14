# PTM - The SHIELDS Particle Tracing Model

### Authors and Contributors
Jesse Woodroffe [1], Steven Morley [1], Thiago Brito [1], Alin-Daniel Panaitescu [2], Michael Henderson [1], Vania Jordanova [1]

1. ISR-1, Space Science & Applications, Los Alamos National Laboratory
2. ISR-2, Space and Remote Sensing, Los Alamos National Laboratory

### About the code
This directory contains the simulation code and analysis tools for the SHIELDS particle tracing model (PTM).

The main source is written in Fortran2008. Associated scripts, pre-processing, and post-processing are written in Python 3.

To build PTM, simply run:
```
make
```

To build the fortran only run `make ptm`, to build/install the Python module only run `make python`, and to convert the Markdown documents to PDF use `make docs`.
The documentation conversion uses the Node.js `mdpdf` module, which can be installed using `npm install mdpdf -g`. PDF documentation is placed in the `docs` directory.