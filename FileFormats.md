# Format for PTM fields input files


## Notes on format

The file is ASCII and will have a filename following the
ptm_fields_XXXX.dat format, where XXXX is a zero-padded number
indicating the timestep.


### Header

The first line of the header gives the (integer) number of grid
points in each dimension, e.g., for a 75x75x75 grid:
```
75 75 75
```

The next 3 lines give the coordinates of the grid points. First
the X-coordinates are given, followed by a newline. Then the Y
and Z coordinates. For a regular cube in (X,Y,Z) centered on the
origin, these three lines will be identical. For a 75x75x75 grid
centered at (0, 0, 0), extending from -6 to 6 in each dimension,
the lines will be as follows:
```
-6.00000e+00 -5.83784e+00 -5.67568e+00 -5.51351e+00 -5.35135e+00 -5.18919e+00 -5.02703e+00 -4.86486e+00 -4.70270e+00 -4.54054e+00 -4.37838e+00 -4.21622e+00 -4.05405e+00 -3.89189e+00 -3.72973e+00 -3.56757e+00 -3.40541e+00 -3.24324e+00 -3.08108e+00 -2.91892e+00 -2.75676e+00 -2.59459e+00 -2.43243e+00 -2.27027e+00 -2.10811e+00 -1.94595e+00 -1.78378e+00 -1.62162e+00 -1.45946e+00 -1.29730e+00 -1.13514e+00 -9.72973e-01 -8.10811e-01 -6.48649e-01 -4.86486e-01 -3.24324e-01 -1.62162e-01  0.00000e+00  1.62162e-01  3.24324e-01  4.86486e-01  6.48649e-01  8.10811e-01  9.72973e-01  1.13514e+00  1.29730e+00  1.45946e+00  1.62162e+00  1.78378e+00  1.94595e+00  2.10811e+00  2.27027e+00  2.43243e+00  2.59459e+00  2.75676e+00  2.91892e+00  3.08108e+00  3.24324e+00  3.40541e+00  3.56757e+00  3.72973e+00  3.89189e+00  4.05405e+00  4.21622e+00  4.37838e+00  4.54054e+00  4.70270e+00  4.86486e+00  5.02703e+00  5.18919e+00  5.35135e+00  5.51351e+00  5.67568e+00  5.83784e+00  6.00000e+00
```

All subsequent rows are data in the following format:
```
i j k Bx[i,j,k] By[i,j,k] Bz[i,j,k] Ex[i,j,k] Ey[i,j,k] Ez[i,j,k]
```
where i, j and k are the index into the grid coordinates.
The indexing starts from 1 (not zero), thus for a 75x75x75 grid
the rows will start with indices
```
   1    1    1 ...
   1    1    2 ...
   1    1    3 ...
   .
   .
   .
   1    1   75 ...
   1    2    1 ...
   .
   .
   .
   75  75  74 ...
   75  75  75 ...
```

All coordinate indices are integer, and all cordinate and data values
are floating point. The usual format is `[-]1.23456e+00`, that is,
six digits plus a two digit exponent. The sign is optionally included
as the leading character.


## Revision History

Earlier versions of SHIELDS-PTM used binary inputs, with each
component of both E and B having separate files. The grid layout
was stored in yet another file. To make it easier to ensure
consistency between files and to make the files more user-friendly
the format was changed to a single ASCII file. This is what is
described in this document.
