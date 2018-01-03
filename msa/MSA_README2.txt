

      Magnetic/Electric Source Analysis 
      Copyright (C) 1997-1999 Robert Oostenveld


This file contains some comments on which part of libmsa is considered 
stable and which part is not. Furthermore some programming tricks and
special (unusual) properties that are used are described, since they 
could lead to portability problems or misunderstanding of the code.

The contents of this file are:
  - stability
  - run-time parameters
  - examples
  - electrodes
  - magneto- and gradiometers
  - debugging
  - error handling
  - coding tricks and remarks
  - currently available functions for EEG/MEG


Stability:
----------
The most functions currently in libmsa are considered stable code and are
not likely to change very much in the future. A couple of functions exist 
however which are still in active development and their names and usage could
change.

In general the small (or larger) functions that are only used inside one
routine should not be used outside of the routine, as they are subject 
to change.
 * mrdf_func()
 * mfdf_func()
 * fdf_func()
 * eeg_mrdf_func()
 * pinvT()

As of version 980318 I have removed the following functions:
 * goalfunction_scan()
 * deviation_scan()
 * eeg_goalfunction_scan()
 * eeg_deviation_scan()
and created a set of similar routines with a different interface. The
calling program now should create the grid on which the scanning
is done. This is done to make it possible to scan more exotic ROI's
(eg. planes or the cortical surface) and to make the current scanning
routines compatible with an implementation of MUSIC scanning and
with a more suitable MNLS routine.

The following routines are possibly subject to change in the (near) future, 
and their names and parameters could change:
 * minimum_norm_least_squares()		/* I don't like the functionality of this one */
 * complex_magnetic_field()		/* this function is quite useless */
 * sensor_array_BTi_Magnes()
 * scan_goalfunction()			/* API may change, should be tested */
 * scan_deviation()			/* API may change, should be tested */
 * scan_music()				/* API may change, code is experimental */
 * eeg_scan_goalfunction()		/* API may change, should be tested */
 * eeg_scan_deviation()			/* API may change, should be tested */
 * eeg_scan_music()			/* API may change, code is experimental */
 * oclim_solve()			/* is untested and very experimental */
 * oclim()				/* is untested and very experimental */
 * dipole_moment_fit()			/* new routine and untested          */


Run-time parameters:
--------------------
There are a couple of parameters which determine the precise behaviour
of MSA diring dipole fitting. These parameters are defined in the general
msa.h header file. In old versions of MSA (pre 990531), it was not 
possible to change these parameters without recompiling the source code.
As of now, I have changed these parameters in global variables. These
global variables control the same behaviour as the older parameters,
but the program which is calling the MSA library is capable to change these
global variables. The most important one is "local_sphere_radius", which 
describes the radius of the local sphere used in MEG dipole fitting. Please
look in the source code (msa_parameters.c) to find out the details.


Examples:
---------
The examples that are given are not very usefull (sorry) but are more thought
to help me in debugging the routines. The example that _is_ usefull is the
dipole_fit_demo program. This program lets you fit a dipole, and gives some 
possibilities to use a simulated dipolar field distribution with added noise.
If you want to know how to use the MSA routines, please take a look in the
dipole_fit_demo source code.


Electrodes:
-----------
Electrodes are described in carthesian space by their (x,y,z)
coordinates. Just as in the case of the MEG, in the EEG the spheres
representing the head have their center in the origin of the coordinate
system. The electrodes should lie on the outermost sphere (last element
of the array), but if this is not the case, they are projected towards
the sphere surface. If the distance between electrode and surface is too
large, eeg_leadfield() will issue a warning message. All the computations
on EEG data are done using average referenced data.


Magneto- and Gradiometers:
--------------------------
The matrices rm and um contain the information for the MEG sensor
array. Since the original implementation was only meant for a crude
representation of a magnetometers with one point for the whole coil,
this is not very convenient, but I kept it for compatibility. The only
routine which really needs this rm and um matrices is the leadfield
calculation. The leadfield is calculated for each sensor position,
i.e. for each digitisation point of the coil. The leadfield() algorithm
senses the number of digitisation points by using a property of the
Numerical Recipes matrix allocation. It is important that you use
the 2nd edition, since the 1st edition did not have this property!
The different possibilities of rm and um which are automatically detected
by leadfield_sphere() are listed below. The last option of the three is
the most general and should preferably be used.

magnetometer described by one point in space:
  rm = matrix(1,Nchans,1,3)
  um = matrix(1,Nchans,1,3)
this means that rm[i] contains the position of the coil of magnetometer i
and um[i] contains the direction of magnetometer i, scaled to unity length

axial 1st order gradiometer described by two points in space:
  rm = matrix(1,Nchans,1,6)
  um = matrix(1,Nchans,1,3)
this means that rm[i] contains the position of the bottom coil and top coil of 
an axial 1st order gradiometer, i.e. rm[i] = (x1,y1,z1, x2,y2,z2) where 
(x1,y1,z1) is the position of the bottom coil and (x2,y2,z2) is the position 
of the top coil. And um[i] contains the direction of the bottom coil of 
the gradiometer, the top coil is in he opposite direction

any other type of magnetometer or gradiometer:
  rm = matrix(1,Nchans,1,m*3) 
  um = matrix(1,Nchans,1,m*3), with m larger or equal to 1
with this representation any linear combination of points in space
describing a magnetic field sensor (magnetometer or n-th order axial
or planar gradiometer, also with physical dimensions) can be used.
Each point needs 3 numbers which describe its position, and three points
which describe its direction. The length of the direction-vector is the
weighing factor.  For example to improve accuracy one could describe
a magnetometer coil by 7 points, each with an apropriate weighing
factor. For a n-th order gradiometer each coil can be represented by a
single point, and the number of loops and the direction of a particular
coil can be represented with the direction and length of the weighing
vector.


Debugging:
----------
For debugging purposes I wrote some functions that can be included in the
source code using #define statements. On this way it is relative easy to
turn debugging output on or off, and to determine where debugging output
is written to. The actual routines are debug(), nodebug(), logdebug(). The
routine logdebug() is the most usefull since it can be turned on and
off using an environment variable MSADEBUG, but probably is also the
slowest. These debugging functions are called just like printf(). If you
are not interested in the huge streams of information generated by the
debugging, please use nodebug(), since this uses least of your computing
resources and makes the msa library faster. The debug() and nodebug()
routines have been implemented using (more efficient) macros instead of
functioncalls. This seems to be incompatible at least with the MS Visual
C++ 5.0 compiler, so for the Windows platform the debug() and nodebug()
routine are available as functions instead of macros.


Error handling:
---------------
The initial implementation of MSA was designed for a simple user
interface with no graphical feedback. For myself I am mainly using MSA in
command line executables with a lot of command line options, but without
user interaction inside the programs.  If an error occured inside an MSA
function, it just called exit(), and terminated the program with that.
This has always worked fine with me, but it is very inconvenient in a
program which has a lot of user interaction. If you would implement a
graphical interface, you don't want it to exit just because some data
was not initialised properly, in that case you would want to give the
user some feedback and the possibility to load/generate the missing data
and try again. 
For this purpose, I have writen a new error-handling method (as of
version 990529). The default behaviour is that it exits, which was the
old behaviour. But the alternative is that an error flag (msa_errno),
which is a global variable, is set and that the routine returns with a
-1 return value (in the pre 990529 versions most functions were of type
"void"). Implementing MSA functions in a graphical program then requires
that the programmer checks the return value after each function call. If
the return value indicates an error, the programmer can find out what
went wrong using msa_strerror(msa_errno). This returns a pointer to a
string containing the error message.  Alternatively the programmer can
call msa_perror(char *), which prints this error message to stderr. These
two functions are modelled after strerror and perror, which are ANSI/ISO
C. Please see your own documentation of these two functions to see how
they are used.
To enable this nicer (but non-default) error handling, you have to compile
the MSA library with the RETURN_ON_ERROR defined.


Coding tricks and remarks:
--------------------------
The 2nd edition of the Numerical Recipes library is needed, since there was
a change in the memory allocation inside matrix() that affects the proper 
working of some of the libmsa functions (concerning the determination of the 
number of coils without passing this explicitly as a parameter).

Usually the functions are called by value for integer and float scalars (and 
called by reference for the vectors and matrices of course). A couple of 
exceptions exists however. 

It should be noted that the macro's (defined in vector.h) that are used widely 
throughout the code should stand within proper {}'s, since they represent 
multiple lines of code. Do not use them as "if (...) CP_VEC(...);" but as 
"if (...) {CP_VEC(...);}".


Currently available files:
--------------------------
The following files are present in the src directory. For the most of them you
can determine from the filename what the implemented routine or functions 
inside the file will do. In the following list the first column contains an 
"E" if the function supports eeg and contains a "M" if meg is supported. The 
files marked with an "*" are currently only supported for meg, but should be 
extended to eeg also. If two files exist with the same functionality for eeg 
and meg data, they are listed next to each other.

  EM	corrcoef.c
  EM	gof.c
  EM	rdm.c
  EM	mag.c
  EM	rms.c
  EM	frequency_component.c
  EM	covariance.c
  EM	covariance_eigenvalues.c
  EM	signal_space_projection.c
  EM	oclim_solve.c
  EM	oclim.c
  EM	pinvT.c
  EM	debug.c
  EM	msa_error.c
  EM	msa_parameters.c
  EM	msa_version.c
  EM	coordinate_transformation.c
  EM	electric_potential.c		magnetic_field.c
  EM	eeg_dipole_moment.c		dipole_moment.c
  EM	eeg_lambda.c			lambda.c
  EM	eeg_leadfield.c			leadfield.c
  EM	eeg_moving_dipole_fit.c		moving_dipole_fit.c
  EM	eeg_rotating_dipole_fit.c	rotating_dipole_fit.c
  EM	eeg_scan_deviation.c		scan_deviation.c
  EM	eeg_scan_goalfunction.c		scan_goalfunction.c
  EM	eeg_scan_music.c		scan_music.c
  EM	eeg_weeder.c			weeder.c

  E-	re_reference.c
  E-	electric_potential_sun.c
  
  -M	frequency_dipole_fit.c
  -M *	fixed_dipole_fit.c
  -M *	dipole_moment_fit.c
  -M *	minimum_norm_least_squares.c
  -M *	demo_data.c
  -M *	sensor_array.c

