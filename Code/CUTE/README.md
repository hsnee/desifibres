*CUTE - Correlation Utilities and Two-point Estimates


1 Introduction.

CUTE (and its GPU implementation CU_CUTE) provide a set of tools to
calculate 2-point correlation functions of cosmological point distributions.
The functionality of both programs is very similar, so we will discuss CUTE
first and the differences with CU_CUTE will be explained afterwards (section
9).

For a description of the functionality of CUTE_box, the implementation of
CUTE for periodic-box catalogs, see the README file contained in CUTE_box's
root directory.


2 Dependencies.

Two things are necessary to compile and run CUTE: the gcc C compiler and the
GSL libraries. Although CUTE is parallelized for shared memory machines with
OpenMP, systems lacking this will be able to run in sequential mode anyway.

Once the input files are provided and the compile-time options set,
compilation should be straightforward by just typing

  $ make CUTE

This will generate the executable CUTE which can then be run by typing

  $ ./CUTE <param_file>

Where <param_file> is the path to the parameter file (see section 3.1).


3 Input files.

Several files are needed in any CUTE run:
3.1 The parameter file:
   Part of the behavior of CUTE is governed through the parameter file. This
   file must contain the following fields:
    * data_filename= FILENAME
              See section 3.2 below.
    * random_filename= FILENAME
              If set to "none" CUTE will generate its own random catalog
	      from the parameters below, otherwise one random catalog will
	      be read from FILENAME and used. All the subsequent parameters
	      regarding the random catalog will be overriden. See section 3.3
	      below.
    * input_format= INT
              Format of the data and random catalogs (both should use the
	      same format). Four formats (0, 1, 2, 3 and 5) are currently supported
	      (see section 3.2 for details).
    * num_lines= INT
              This parameter determines the number of lines to be read from
              the data and random catalogs. If absent or set to 'all', all
 	      the lines in the files will be read.
    * output_filename= FILE
              See section 5 below.
    * z_dist_filename= FILENAME
      	      Only relevant if random_filename!="none". See section 3.5 below.
    * mask_filename= FILENAME
      	      Only relevant if random_filename!="none". See section 3.4 below.
    * corr_estimator= STRING
              Estimator to be used. Allowed values for STRING are:
                      - PH  -> w=DD/RR-1         (Peebles & Hauser estimator)
                      - DP  -> w=DD/DR-1         (Davis & Peebles estimator)
                      - HAM -> w=DD*RR/DR^2-1    (Hamilton estimator)
                      - HEW -> w=(DD-DR)/RR      (Hewett estimator)
                      - LS  -> w=(DD-2*DR+RR)/RR (Landy & Szalay estimator)
    * np_rand_fact= INT
              If the data catalog has N_DAT objects, the random catalogs will
              have np_rand_fact*N_DAT objects.
    * omega_M= DOUBLE
              Matter parameter (cosmological parameters are needed for the
              distance-redshift relation, which is only necessary for the
	      monopole and 3-D correlation functions).
    * omega_L= DOUBLE
              Dark energy parameter
    * w= DOUBLE
              Dark energy equation of state.
    * corr_type= STRING
              Type of correlation function to be computed. Five options are
	       available:
 		      "radial"   -> Radial correlation function
 		      "angular"  -> Angular correlation function
                      "monopole" -> Monopole
                      "3D_ps"    -> 3D correlation binning in sigma-pi
                      "3D_rm".   -> 3D correlation binning in r-mu
              See section 6 below for more information on the different types.
    * radial_aperture= FLOAT
              Maximum angular separation (in deg) between pairs of galaxies
	      to calculate the radial correlation function. See section 6.
    * use_pm= INT
              If set to 1 the angular correlation function will be calculated
	      by pixelizing the catalog. Set it to 0 otherwise (brute-force).
    * n_pix_sph= INT
              The sphere will be divided in 2 x n_pix_sph^2 pixels. Only
	      relevant for the angular correlation function with pixelization.
              See section 8 below.
    
    Any blank line or any line starting with the character '#' in the
    parameter file will be ignored. An example of this file is provided with
    the test suite.

3.2 The data catalog:
    CUTE assumes catalogs are ASCII files with 3 or 4 columns:
                 x0         x1        x2        x3
    The meaning of the first three columns will depend on the value assigned
    to the parameter "input_format":
       input_format= 0 : x0 = redshift, x1 = cos(theta), x2 = phi
       input_format= 1 : x0 = redshift, x1 = DEC,        x2 = RA
       input_format= 2 : x0 = RA,       x1 = DEC,        x2 = redshift
       input_format= 3 : x0 = redshift, x1 = RA,         x2 = DEC
       input_format= 5 : hdf5 contain z, ra and dec with angles in degrees
    Here, cos(theta) and phi are the spherical coordinates corresponding to
    each object (with phi in radians). These spherical coordinates are
    related to equatorial coordinates by RA=phi and theta=90-DEC. For formats
    1 and 2, RA and DEC should be given in deg. The fourth column x3 is only 
    necessary if CUTE is compiled using the -D_WITH_WEIGHTS flag in the
    Makefile. In this case the fourth column must be a floating point number
    containing the corresponding galaxy's weight.

    This file may contain other columns, which will be ignored, but the
    first 3 or 4 must correspond to the ones above. In particular, note
    that one may have a fourth column containing the galaxy weights that
    will be ignored if CUTE is compiled ***without*** -D_WITH_WEIGHTS.
    If other input formats are necessary, the user may modify the function
    read_line in src/io.c, in order to accommodate it to their needs.

3.3 The random catalog:
    By default CUTE generates the random catalogs, but users will often
    want to use their own. For this reason the option random_filename is 
    provided in the parameter file. If set to "none" CUTE will generate its
    own random catalogs, otherwise it will read the random catalog from
    the corresponding file. The expected file format is the same as for the
    data catalog (above). If the random catalog is generated by CUTE, all
    random points are given a weight of 1.

3.4 The mask:
    The format of the mask file should be a series of lines with 6 columns:
                z0     zf     cth0     cthf     phi0   phif
    Each line defines a cube in the space of redshift and spherical coords,
    with  z0<z<zf, cth0<cos(theta)<cthf, phi0<phi<phif. The total mask is the
    logical sum of all these cubes. Note that this file is only necessary to
    generate random catalogs, and is therefore ignored if the random catalog
    is read from a file (i.e. random_filename != none).

3.5 The redshift distribution:
    The objects in the random catalogs generated by CUTE need to follow the
    same redshift distribution as those in the data, thus a file with this
    redshift distribution is needed as input. This file must contain two
    columns:  
                      z              N(z)
    IMPORTANT: CUTE assumes that the values of z are monotonously
    increasing. The normalization of N(z) is irrelevant. Note that this file
    is only necessary to generate random catalogs, and is therefore ignored
    if the random catalog is read from a file (i.e. random_filename != none).


4 The compile-time options:

For performance reasons some run parameters must be chosen through
compile-time options defined at the beginning of the Makefile. These are:

   >For the radial (redshift) correlation function:
   * DZ_MAX=<> -> maximum Dz to which the radial correlation function is 
                   calculated.
   * NB_DZ=<>  -> #bins for the radial correlation function.

   >For the angular correlation function:
   * THETA_MAX=<> -> maximum angle (in deg) to which the angular correlation
                     function is calculated.
   * NB_THETA=<>  -> #bins for the angular correlation function.

   >For the monopole correlation function:
   * R_MAX=<> -> maximum radius to which the monopole correlation function
                 is calculated (in current length units).
   * NB_R=<>  -> #bins for the monopole correlation function.

   >For the 3-D correlation function (binning in sigma-pi):
   * RL_MAX=<> -> maximum distance to which the 3D correlation function
                  is calculated in the longitudinal direction.
   * NB_RL=<>  -> #bins in the longitudinal direction for the 3D correlation
                  function.
   * RT_MAX=<> -> maximum distance to which the 3D correlation function
                  is calculated in the transverse direction.
   * NB_RT=<>  -> #bins in the transverse direction for the 3D correlation
                  function.

   >For the 3-D correlation function (binning in r and mu):
   * R3D_MAX=<> -> maximum radius to which the 3D correlation function is
                   calculated.
   * NB_R3D=<>  -> #bins in r for the 3D correlation function.
   * NB_CTH=<>  -> #bins in mu for the anisotropic correlations
   
   >Logarithmic binning can be selected for the angular (without
    pixelization), monopole and 3D (in r-mu) correlation functions. This
    option is chosen through the flag -D_LOGBIN (see below). If set, the
    number of bins per decade can be set through the variable:
   * N_LOGINT=<> -> # bins in r (or theta) per decade.
 
   >Behavior options: add any of these to the variable DEFINEOPTIONS.
   * -D_WITH_WEIGHTS -> read and use galaxy weights both for the data and
                        random catalogs. The use of weights is only
			implemented in the OpenMP version of CUTE.
   * -D_VERBOSE      -> extra info will be output.
   * -D_DEBUG        -> creates some debugging files.
   * -D_TRUE_ACOS    -> the true arccosine will be computed for the
                        angular correlation function (see below for more
		        details).
   * -D_LOGBIN       -> logarithmic binning will be used. Only relevant
                        for angular, monopole and 3D (in r-mu) 2PCFs. Note
		        that performance will be poorer when using
		        logarithmic binning.

   >The Gnu Scientific Library (GSL) must be installed in order to be able
    to compile CUTE. Set the Makefile variables GSL_INC and GSL_LIB to the
    paths to the GSL headers and libraries.

   >CUDA options: see section 9 below.


5 The output file.

The output file contains 6 columns for the radial, angular and monopole
correlation functions with
              x   xi(x)   error(x)   DD(x)   DR(x)   RR(x) 
where x is Dz, theta or r. The third column is the Poisson error calculated
from DD, DR and RR.

For the 3-D correlation functions the output file has 7 columns with
   x1   x2   xi(x1,x2)   error(x1,x2)   DD(x1,x2)   DR(x1,x2)   RR(x1,x2)
where (x1,x2) is either (pi,sigma) or (mu,r).

For the monopole and 3-D 2PCFs distances are given in units of Mpc/h.

*** IMPORTANT *** Note that the error column only contains the Poisson
contribution to the total error and is therefore an extremely incorrect
estimator of the total error.


6 The correlations.

CUTE supports the calculation of 5 different types of 2-point correlation
functions. In all cases the process is always very similar:
  - Take a pair of objects in the data catalog.
  - Calculate the distance between both objects.
  - Bin the pair in a histogram according to the calculated distance.
  - Repeat for all pairs in the catalog to obtain the histogram DD.
  - Do the same for a random catalog to obtain RR.
  - Do the same but with each pair composed of a data object and a random
    object to obtain DR.
  - Use DD, DR and RR to estimate the correlation function through the
    estimator defined in the options file.
Each type of correlation function has some subtleties regarding this
calculation, which we detail below.

 * The angular correlation function
   To calculate the angular correlation function, object positions are
   normalized (i.e. projected on the sphere). Then the cosine of the angle
   separating them is calculated through the dot product of their position
   vectors. The separation angle is calculated then by taking the arccosine
   of this dot product. Taking the arccosine is computationally very slow
   and by default the code uses the following approximation:
       arccos(1-x)~sqrt(2*x+(x^2)/3+4*(x^3)/45)
   which is a good to 0.01% for angles below 40 deg and speeds up the
   calculation by a factor 2. This can be overriden by defining the
   compile-time flag -D_TRUE_ACOS, in which case the C library arccosine
   function acos will be used.

   To speed up the calculation, avoiding unnecessary correlations, the 
   whole sphere is divided into cubes in spherical coordinates (pixels).
   These pixels are used to define sub-regions of the sphere which can
   be used for fast neighbor-searching. This is described in more
   detail in the code's preprint (arXiv:1210.1833). This technique, which
   was implemented in version 1.0, is also used to search for aligned
   pairs in the calculation of the radial correlation function (see below).

   CUTE also has the option of calculating the angular correlation function
   by pixelizing the catalog (i.e., interpolating the objects to a grid and
   then correlating the pixels instead of the objects). This option is much
   faster when your catalog is large, but of course will only yield reliable
   results down to the pixel angular resolution. This option is selected by
   setting the use_pm variable in the options file to 1.

 * The radial correlation function
   A detailed description of what we call the radial correlation function
   here is given in arXiv:1210.6446. This is the fastest but most involved
   correlation. First, the angular separation between a given pair of 
   galaxies must be calculated. If this angle is smaller than the maximum
   allowed radial aperture (given by the parameter radial_aperture in the
   parameter file), the galaxies are considered collinear, and are correlated
   according to their separation in redshift Dz.

 * The monopole correlation function
   In this case object pairs are binned simply according to the
   three-dimensional distance between them. A strategy similar to the one used
   for the angular correlation function is used in this case (as well as in the
   next two). First the catalog is placed inside a box, which will not
   necessarily be cubic. This box is then divided into cubic cells, which are
   used for neighbor-searching, saving up unnecessary correlations. This
   feature is new from version 1.0.

 * The 3-D correlation function (binning in sigma-pi).
   For each pair of objects their relative separation vector (x_r) and their
   center of mass (x_CM) are calculated. Then the radial and transverse
   separations are computed as
            pi    = |x_r.x_CM|/|x_CM|
            sigma = sqrt(|x_r|^2-pi^2)
   Pairs are binned in a two-dimensional histogram by these quantities.

 * The 3-D correlation function (binning in r-mu).
   For each pair of objects their relative separation vector (x_r)
   and their center of mass (x_CM) are calculated. Then the relative
   separation and deviation w.r.t. the longitudinal direction are
             r  = |x_r|
             mu = |x_r.x_CM|/(|x_CM| |x_r|)
   Pairs are binned in a two-dimensional histogram by these
   quantities.


7 Nearest-neighbour searching

Usually the range of scales of interest, in the case of the three-
dimensional or angular correlation functions, is considerably smaller
than the size of the catalog. Thus, a good way to speed up the
calculation is to avoid, in as much as possible, correlating pairs
of objects beyond the maximum scale. To do this one must find a way
to collect all particles within a given distance of another one without
actually having to calculate the distances between all pairs.

The strategies used by CUTE for this purpose for the 3D and angular 2PCFs
are very similar, so we will describe the three-dimensional case first.
A box encompassing the whole catalog is first determined and divided into
cubical cells. To each cell we associate the positions of all the objects
that fall inside it. Assuming that the maximum distance we are interested
in is Rmax and the cell size a_cell, we draw a cube of
                      2*(int)(Rmax/a_cell)+1
cells around every cell Ci, and we correlate all the objects in Ci with
all the others inside this cube. Of course, the objects close to the corners
of the cube will be beyond Rmax, however this is a small price to pay for 
a fast method to calculate the nearest neighbours.

In the case of the angular correlation function, we can define similar
"cubes" in spherical coordinates of size > theta_max. The sides of these
cubes for a point at (theta,phi) are:
   D(cos(theta)) = cos(theta - theta_max) - cos(theta + theta_max)
   D(phi) = 2 * sqrt(cos^2(theta_max) - cos^2(theta)) / sin(theta)
These regions can be used to find nearest neighbours in a way similar to
the one described above, this time using angular pixels.


8 Pixelization

Angular pixels are used by CUTE to calculate the radial and angular
correlation functions. CUTE creates its own pixels by dividing the
sphere in "squares" in spherical coordinates, i.e.: equal intervals
in cos(theta) and phi. The size of these pixels is determined by the
variable n_pix_sph. The sphere is divided into n_pix_sph bins in
cos(theta) and 2 x n_pix_sph bins in phi. For instance,
n_pix_sph=2048 yields pixels with an angular resolution of 0.07 deg.


9 The CUDA implementation.

The current version of CUTE comes with a GPU version (CU_CUTE)
written in CUDA-C. There exist several differences with respect to the
OpenMP version both in compilation and usage. We enumerate them here:

- The radial correlation function is not supported in CU_CUTE. This
  is not a priority, since it is numerically a much faster problem.
- In order to use CU_CUTE we must (obviously) have a CUDA-enabled
  device in our system. The current implementation assumes compute 
  capability 2.0 or higher for this device (atomic operations must be
  possible both in global and shared memory). Furthermore, the NVIDIA
  compiler driver NVCC must be installed.
- The CUDA libraries are assumed to be in /usr/local/cuda. If this is
  not the case, change the variable CUDADIR in the Makefile accordingly.
- CUDA kernels are launched with a number of thread blocks given by
  the variable n_blocks, which is taken as a command-line option (see below).
  The optimal value for this number is both data and hardware dependent, so
  some degree of trial-and-error may be necessary to find the best choice.
  For example, for an NVIDIA Tesla C2070 GPU and a catalog with 300K objects
  ~128 blocks seem to be optimal.
- For performance reasons the default compilation options for CU_CUTE
  include less precise divisions and square roots. If a good numerical
  precision is required on the GPU this can be disabled by commenting out
  the variable OPT_PRECISION in the Makefile.
- IMPORTANT: The number of histogram bins, which could be changed in CUTE
  through the compile-time options (NB_THETA,NB_R,NB_RL,NB_RT,NB_R3D,NB_CTH)
  are now hard coded (since they define the number of threads per block in GPU
  kernel launches) to 256. This shouldn't be a problem, since you can get the
  result for smaller values by re-binning these. Two options are given for 2D
  histograms: 128x128 and 64x64 bins, which must be likewise selected through
  the compile-time option NB_H2D. If so desired, these values can be changed
  by modifying the header file define.h in the folder ./src/. This may
  however affect performance or even prevent the code from working.
- The size of the object catalog must be large enough to fit inside the
  device's global memory.

Once the input files are provided and the compile-time options set,
compilation should be straightforward by just typing

  $ make CU_CUTE

This will generate the executable CU_CUTE which can then be run by
typing

  $ ./CU_CUTE <options file> <number of CUDA blocks>


10 Test suite

The current release (>0.2) comes with a folder called "test" containing
a set of files. These are:
 - shell.dat  : a test catalog with objects randomly distributed in the
                region 0 < DEC < 30 and 0 < RA < 45 with redshifts 
		around z=0.7 with a gaussian window function of
                sigma_z=0.
 - random.dat : a corresponding random catalog
 - mask.dat   : the mask file corresponding to the catalog
 - dndz.dat   : the redshift distribution for the mock catalog
 - param.txt  : an instance of the parameter file needed to run
                CUTE on the test catalog
This test suite is meant to clarify the expected format of the different
input and output files and also as a check that CUTE is running smoothly.
The parameter file is set up so that CUTE can be run on the test catalog
from the root file by just typing:

  $ ./CUTE test/param.txt

or

  $ ./CU_CUTE test/param.txt n_blocks

(n_blocks=256 should work well on a catalog of this size)


11 License.

CUTE and CU_CUTE are distributed under the GPL license (see COPYING in
the root directory). We kindly ask you to cite the program's website 
http://members.ift.uam-csic.es/dmonge/CUTE.html and accompanying
preprint arXiv:1210.1833 when using it for published results.


12 Contact.

Regarding bugs, suggestions, questions or petitions, feel free to contact
the author:
    David Alonso: david dot alonso at uam dot es 

*Edit*

Edited to include Hawkins-Maddox weighting to correct small-scall
correlation function. Implemented into angular, monopole and xi(sigma,pi).
Routine at src/woftheta.c
