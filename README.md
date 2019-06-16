# Code-for-parallel-Recursive-Newton-s-Algorithm-using-MPI
we have just coded to implement parallel Newton’s recursive algorithm. To accomplish this, we first  partition the problem and fitted each partition into master-slave by MPI mode. First, we need to have: 
(a) The number of interpolation points.
(b) A grid of points at which the interpolation is to be exact.
(c) An array containing the function we wish to interpolate evaluated at the interpolating grid. 
(d) An array to store the Newton differencing coefficients
