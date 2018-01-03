The Matlab BVP-package

TABLE OF CONTENTS
*) What is it all about ?
*) Installation
*) This package's files
*) Some details concerning program and algorithm
*) Hints for troubleshooting
*) Where to go from here ?
*) Whoami ?

***********************************************************************
***********************************************************************

What is it all about ?

The BVP-package contains functions for solving nonlinear ODE-boundary
value problems of first and second order on nonequidistant grids.
Dimensions of the systems are arbitrary.

For those of you that do not like reading long texts, install the
package using the installation instructions below and call BVPDEMO
to receive an overview of the capabilities of this package.

***********************************************************************
***********************************************************************

Installation

Create a directory for the BVP-package and unzip BVPs.zip into this
directory. Then add it to the Matlab search path. This is done by
inserting the line

addpath('The full path of the directory you created');

into the Matlab startup-file startup.m. You can locate this file by
typing 

>> which startup

at the Matlab command line.

Now add search-paths for the subdirectories PRIVATE and DEMO, too.
When you call startup or restart Matlab all BVP functions will be
available.

***********************************************************************
***********************************************************************

This package's files:

BVP              main solving routine for first order problems
BVP2             main solving routine for second order problems
BVPSET           tool for setting parameters - very much like ODESET
BVPPLOT          output routine - very much like ODEPLOT
BVPPHAS2         output routine - very much like ODEPHAS2
BVPPHAS3         output routine - very much like ODEPHAS3

subdirectory private:

BVP_DISC         These files are used by BVP and BVP2, respectively.
BVP_DISC_JAC     Unless you want to understand the algorithm,
BVP2_DISC        you do not need them.
BVP2_DISC_JAC

subdirectory demo:

BVPDEMO          extensive demo-file
BVPDEMO1         files used by BVPDEMO
BVPDEMO2         
BVPDEMO3
BVPDEMO4
BVPDEMO5
BVPDEMO6

Look up the source code of the demo files if you encounter difficulties
when defining your own BVPs.

***********************************************************************
***********************************************************************

Some details concerning program and algorithm

This program was written as a part of a project I completed at the 
Institute for Applied and Numerical Mathematics at the University
of Technology, Vienna, Austria.
Its main goal was to find out, whether the homoclinic orbit of the demo-
file could be found numerically, although the problem is not well posed.

The finding of the orbit with both first and second order discretisations
encouraged to search for homoclinic orbits numerically in related, but more
complicated problems, where existence of such orbits cannot be proven (yet).

------------------------------------------------------------------------

The discretisation for first order systems reads

      x'=f(t,x) , t \in [a,b] , R(x(a),x(b))=0

           x1-x0=h1/2*(f(t0,x0)+f(t1,x1))

Trapezoidal rule, local discretisation error=O(h^2).

------------------------------------------------------------------------

The discretisation for second order systems reads

      x"=f(t,x,x') , t \in [a,b] , R(x(a),x'(a),x(b),x'(b))=0

    h2/h*x1 -x2 -h1/h*x3 + (h1*h2)/2*f((t1+t3)/2,(x1+x3)/2,(x3-x1)/h)


This is a generalisation of the symmetric difference discretisation
to nonequidistant grids (local discretisation error=O(h^2)).
Due to use of additional points at the ends of the 
time-interval we have l.d.e.=O(h^2) at the boundary too.

------------------------------------------------------------------------

I did not concentrate much on performance matters (you will not use Matlab
for high performance computing), but rather on simple usage and a clear 
algorithm that makes debugging simply.
If you are familiar with the Matlab ODE commands, I dare say you are familiar
with the BVP commands too.

***********************************************************************
***********************************************************************

Hints for troubleshooting

If you have troubles solving your own problems here is a checklist
you might want to follow.

* Check if the problem was formulated correctly, i.e the definition
  file returns correct outputs.
  
* If you do not find any errors, solve the problem with fsolveoption(9)=1;
  (See BVPSET if you do not know how to do this)
  This will cause a comparison of the analytic derivatives you supply
  with a finite difference approximation of these derivatives.

* If there is no difference between those two, we can assume that you
  defined your problem correctly. To obtain a better understanding of what
  is going wrong, set fsolveoptions(1)=1 to receive some information about
  the iteration process.
  It may be helpful to trace the iteration step by step using an output 
  function (See BVPSET).

* Check if your problem is well posed, i.e. if there exists a solution
  of the analytic problem, that is (at least in a neighbourhood
  of this solution) unique and depends CONTINOUSLY on the problem
  parameters (boundary values). 

* Try to find a better starting approximation - this should be especially
  useful, if you 'get stuck in a minimum, that is not a root'.

* Curse the author

***********************************************************************
***********************************************************************

Where to go from here ?

* Run BVPDEMO to get an impression of the capabilities of the package.

* Read the help texts to BVP and BVPSET to understand how to use the package.

***********************************************************************
***********************************************************************

Guenter Kneisl
University of Technology, Vienna
e9425595@fbma.tuwien.ac.at

October 1999








