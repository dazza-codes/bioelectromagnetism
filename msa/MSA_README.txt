

      Magnetic/Electric Source Analysis 
      Copyright (C) 1997-1999 Robert Oostenveld


Description:
------------
MSA is a C library which contains a collection of routines for analysis
of EEG and MEG* data. Among others it includes moving, rotating and
fixed dipole fitting with one or more dipoles inside a homogenous
sphere. The functions in the MSA library use the "Numerical Recipes in
C" library. The functions can easily be integrated within any program,
which takes care of reading EEG or MEG data from a file and presenting
the results to the end user. The functions use the Numerical Recipes
style vectors and matrices as input and output parameters.  The functions
support MEG systems with magnetometers and with 1st order gradiometers.

(*) MEG means MagnetoEncephaloGraphy, it is the measurement of the magnetic 
    field of the human brain, which is comparable to EEG which mesures the 
    electric potential of the brain. MEG and EEG are used to study the activity 
    and functioning of the human brain.


Documentation:
--------------
There are a couple of ascii files included with the MSA distribution. These
are listed below:
  README	the file you are reading right now
  README2	additional and more technical information
  INSTALL	describes the requirements and installation instructions
  COPYING	the GNU general public license
  HISTORY	description of all the changes made during development
  TODO		some plans and desires that I have and I use it as a reminder 
Furthermore you can find formated documentation (html and manpages) in the doc
directory. These are automatically created from the comments in the c-code, and 
therefore are not formatted optimally. For myself I find it the most convenient
to use the comments in the self explaining c-code files when I want to look up 
how a certain routine works exactly and what parameters it needs.


How it came to being:
--------------------
I started with programming some of the MSA functions since I wanted
to learn how the dipole analysis of MEG data is done. Another reason
was that I was not satisfied with the possibilities for data analysis
that were available to me, and I wanted to have things more in my own
hands. Once I learned that the dipole analysis is actually not that
hard, I decided to add more functionality.  My main goal was to get a
collection of functions for MEG data analysis, which could then be easy
incorporated into C programs written by others.
Therefore I have tried to write the functions consistent and properly
documented. Since I started working in Nijmegen I am getting more
involved with EEG measurements, so I extended the library with EEG
analysis routines which are very similar to the MEG routines.


Implemented functions:
----------------------
All the MEG functions are writen with a 37-channel BTi system in mind,
but there is no reason that they won't work with any other system (with
possible exception of the initial guess creation which is now limited
to a halfsphere). The EEG functions are implemented for a model with
4 concentric spheres, but less spheres are of course also supported.
The following functions have been implemented and tested in libmsa:
   * forward calculation of magnetic field
   * forward calculation of electric potential
   * single moving/rotating/fixed dipole fit
   * multiple moving/rotating/fixed dipole fit
   * linear and nonlinear estimation of dipole moment
   * creation of initial guess for dipole fit (only MEG, not for a whole head system!)
   * goalfunction/deviation scan (scanning of a dipole on a predefined grid)
   * goalfunction/deviation scan with fitted dipoles before scanning
   * signal space projection for distibuted and localised sources
   * calculation of the eigenvalues of a covariance matrix
   * calculation of goodness of fit and residual variance
   * calculation of correlation coefficient
   * root mean square calculation
   * some general functions for coordinate transformations and rotation
   x single or multiple dipole fitting in frequency domain
   x minimum norm and optimal constrained least squares inverse method
The dipole fit routines can also be used for a weighed dipole fit, in which the
weighing can be done with the channel noise values (inverse variance) or with
the inverse covariancematrix of the measured data. The routines marked with "x"
are currently in development and may or may not be included in the future.


Support of EEG data:
--------------------
Although I started implementing the MSA library for MEG measurements
and simulations, I have extended it to also support EEG. I implemented a
model with 4 homogenous concentric spheres, representing brain, CSF, skull
and skin. Using the spherical model, the MEG is unaffected by different
conductivities and radii, and therefore the radius and conductivity
of the local sphere are not specified in the MSA functions. The EEG is
strongly affected by the different conductivities of the tissues involved,
and the MSA functions regarding EEG differ in this respect as they do
require the radii and conductivities. The electrodes are represented
by their carthesian coordinates (x,y,z), and they should lie on the
outermost sphere. The possibilities for spatiotemporal dipole fitting
are certainly not as extended as with commercial software (Besa), and
additional constraints on the fit (eg. left-right symmetry) are also
not yet supported.


How about the local sphere?
---------------------------
All the routines have been written for dipoles in a single or multiple
concentric homogenous spheres. This sphere is assumed to be in the
center of the coordinate system and the radius is defined in the msa.h
header. How now to do dipole fitting when this sphere is somewhere
else? Take the positions of the sensors (matrix rm) and shift them
so that the local sphere lies in the origin. Call the dipole fitting
function and after that is ready, shift the dipole positions back (with
the sensor positions) to the original position. And youre done! The same
goes for the (carthesian) position of the electrodes in case of EEG.
There are two easy-to-use routines available which take care of the
shifting of the local sphere coordinate system to the origin and back
again.  The reason to implement this on this particular way is that this
is much more efficient in multiple calls to the same (moving!) dipole
fit function, instead of doing this shifting of the sensors during each
function call.

You also might wonder about the radius of the local sphere: as a matter
of fact this does not influence the magnetic field produced by a dipole
inside of the sphere. Se regarding the forward calculation of the magnetic
field, you don't have to bother to specify the radius of the sphere. If
you do a dipole fit however, you might want to change the radius of the
local sphere to reflect the size of the head. The radius of the local is
used to determine in which area the dipoles are allowed. During dipole
fitting, the dipole is kept inside the sphere (ie. during fitting, the
distance from each dipole to the origin will never get larger than the
radius of the sphere). Please reed the section on "run-time parameters"
in the second readme file to see how it is done.

For the EEG, you must specify the number of spheres used to model the
head, their radii and their conductivities. Upto 4 concentric spheres
are supported (representing brain/CSF/skull/skin, in this order). In the
EEG case, the center of the sphere is also located in the origin. To
do the fitting with a sphere somewhere else, follow the same method
as in the MEG case: shifting the sphere with electrodes to the origin,
do the fitting and shift the dipole positions back with the same amount.


Portability:
------------
The routines have been developed with portability in mind using the GNU 
C-compiler and GNU make. The routines use the Numerical Recipes in C, 2nd 
edition library. The library of the 1st edition does not work due to 
another allocaltion strategy in their matrix() routine, but maybe the 1st 
edition would work if combined with the 2nd edition nrutil.c routines.

The main development (programming and testing) of the MSA library is
done on i486/i586 pc's running GNU/Linux with kernel version 2.0.35 and
2.2.5. The library compiles (using gcc 2.7.2) with no known problems on
the following machines:
   * i486/i586 GNU/Linux 2.0.29 and later
   * Solaris 2.5
   * SunOS 4.1.3
   * IBM AIX 3.2.5 
I have done some (but not extended) testing to see if it also compiles
under DOS and Windows 95. Using the DJGPP compiler (DOS version of gcc),
no problems were encountered. Using the Microsoft Visual C 5.0 compiler,
there appear to be some strange incompatibilities which I believe to be
due to some Windows 95 quirks. Any help to get it running under DOS or 
Windows using any of the common (Borland, Microsoft, DJGPP) compilers
would be appreciated.


Acknowledgments:
----------------
I thank Andreas Wollbrink for implementing the routines with the Muenster 
file format (*.mfx) and testing parts of them.


Legal Issues:
-------------
This software is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This software is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this software. If not, write to the Free Software
Foundation, 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

-------------------------------------------------------------------------


Have fun with it and please let me know if you are finding it usefull
or if you are extending the library or find any bugs. I would like to
increase the functionality and let it grow, so I would realy 
appreciate your feedback and help! 

                                                    Robert Oostenveld

roberto@mbfys.kun.nl
R.Oostenveld@czzoknf.azn.nl


