GEMA
====

*GEMA* is a modular and extensible software for [random
walk](http://en.wikipedia.org/wiki/Random_walk) and other [Monte
Carlo](http://en.wikipedia.org/wiki/Monte_Carlo_method) Simulations
written in the [Clojure](http://clojure.org) programming language.

The software is built around the following main elements:

-   An engine which takes care of initializing the modules according to
    the command line arguments and running the calculations
    concurrently.
-   Modules called *models* which implement different types of
    simulations by imposing different rules for particle coordinate
    initialization and the displacements for each step.
-   Modules called *trackers* which can calculate certain values at each
    step of the simulation.
-   Modules called *reducers* which can calculate certain values using
    all the stored coordinates and tracker values, once the random walk
    is finished.

The engine and the modules are independent, but the engine "knows how to
use" the different types of modules thanks to well defined interfaces. A
module consists of a namespace and a set of functions with well defined
names and certain rules for what the module functions can return. Each
type of module (model, tracker, reducer) has a well defined interface
that is described in the relevant section of the developer
documentation.

Models, trackers and reducers can be chosen at runtime among the
available modules passing appropiate command line arguments. The models
are a set of rules for the initial coordinates assigned to particles and
the change in the coordinates for each successive step of the random
walk. The trackers are modules which take care of calculating a desired
value for each position of the particle. The software will perform the
calculation independently for a number for particles and store the
results for each. The calculation consists of assigning initial
coordinates to the particles and then moving the particles a number of
times according to the rules defines by the model. At each of these
steps the values for each active trackers are recalculated. Both the
coordinates and the tracker values are calculated for each step, but
only stored for certain values of the number of steps that can be
selected by the user. Once the random walk simulation has been done,
reducers can use all the stored particle data to compute values of
interest.

Installation
------------

Jar files can be created using Leiningen (for example with the command 
`lein uberjar`).

Usage
-----

With jar file:

\$ java -jar gema-0.4.2-standalone.jar [options]

With Leiningen:

\$ lein run -m gema.core [options]

Options
-------

| Switches                       | Default      | Description                                                               |
| ------------------------------ | ------------ | ------------------------------------------------------------------------- |
| -h, --help                     |              | Show help and version info                                                |
| -v, --verbosity NUMBER         | 1            | Verbosity level (0-\>Quiet, 1-\>Normal, 2-\>Verbose)                      |
| -q, --quiet                    |              | Suppress most of output                                                   |
| -p, --particles NUMBER         | 1000         | Number of particles                                                       |
| -s, --steps NUMBERS            | "0 1000"     | Step values (example "0 100 1000")                                        |
| -a, --auto-steps               |              | Compute steps automatically using --steps as (start stop every)           |
| -m, --model NAME               | circle       | Model (See available models below)                                        |
| -M, --model-pars PARLIST       | ""           | Parameters for selected model (example "r 1.0 l 4.5")                     |
| -i, --init-model NAME          |              | Positions initialization model (Any available model, defaults to --model) |
| -I, --init-model-pars PARLIST  | ""           | Parameters for positions initialization model                             |
| -t, --trackers NAMELIST        | "spatial"    | Value trackers (See available trackers below)                             |
| -T, --trackers-pars PARLIST    | ""           | Parameters for selected trackers (example "trackername.mypar 0.0")        |
| -r, --reducers NAMES           | ""           | Data reducers (See available reducers below)                              |
| -R, --reducers-pars PARLIST    | ""           | Parameters for selected reducers (example "reducername.mypar 0.0")        |
| -w, --workers NUMBER           |              | Number of concurrent workers                                              |
| -S, --random-seeds NUMBERS     | ""           | Random seeds for PRNGs                                                    |
| -d, --seed-device DEVICE       | /dev/random  | Binary file containing seeds or device to read them                       |

Available models: *acinarduct*, *cylinder*, *circle*, *free3d*,
*sphere*, *segment*, *boxnd*, *pointnd*.

Available trackers: *diffmr1*, *spatial*.

Available reducers: *pi*, *diffmr1adc*.

Examples
--------

Compute the value of the mathematical constant PI using a 2D square
model:

`java -jar gema-0.4.2-standalone.jar --model boxnd --model-pars "dim 2 size [1.0 1.0]" --trackers "" --particles 1e6 --steps "0" --reducers pi`

Run a simulation on a sphere of radius 0.5 with 1000 particles and 10000
steps and store values at steps 0 1 100 1000 10000, use the default
tracker:

`java -jar gema-0.4.2-standalone.jar --model sphere --model-pars "r 0.5" --particles 1e3 --steps "0 1 10 100 1000 10000"`

or using abbreviated options and adding the *diffmr1* tracker:

`java -jar gema-0.4.2-standalone.jar -m sphere -M "r 0.5" -p 1e4 -s "0 1 10 100 1000" -t "spatial diffmr1"`

Run a simulation with a cylinder using the *diffmr1* tracker with a
gradient perpendicular to the cylinder axis, store values every 1000
steps between 0 and 10000:

`java -jar gema-0.4.2-standalone.jar -m cylinder -M "r 0.5 l 4.0" -p 1e4 -s "0 10001 1000" -a -t "diffmr1" -T "diffmr1.grad-dir [1.0 0.0 0.0]"`

Run a 2D simulation in a circle of radius 1.0 with no trackers, store
all the positions from 0 to 10000 for the 10 particles and use PRNG seed
from file seeds.bin:

`java -jar gema-0.4.2-standalone.jar -m circle -M "r 1.0" -p 10 -s "0 10001 1" -a -t "" -d seeds.bin`

Model use and available models
------------------------------

Models are a set of Probability Density Functions and initialization and
displacement rules that are used during the simulation. Model with
dimensions 1-3 can be easily visualized and understood as geometric
structures, in fact many models are created as models of real objects,
in order to understand natural phenomena taking place within them, for
example gas diffusion.

One of the features of this computer program is the selection of any
available model at runtime. Furthermore, it is possible to use one model
for setting the initial positions of the particles and a different model
for restricting the random walk. The initial positions will be selected
according to the initialization model, but only those with coordinates
obeying also the initialization rule of the random walk model will be
allowed. It is the responsibility of the user to adjust model parameters
so that the model share some region of the simulation space, otherwise
the software will try indefinitely to find initial positions without
success (yes, an infinite loop). A possible use of this feature is
studying diffusion, inside a given model, of particles that start
diffusing from a smaller region within the model.

The following models are available with spatially uniform position
initialization inside the model and Gaussian displacements during the
random walk:

-   *acinarduct*: Closed 3D model of the lung acinar ducts composed by
    two concentrical cylinders of identical length with the space
    between the two cylinders divided in cells that commmunicate to the
    central section by apertures. Accepts the parameters: r-in, r-out,
    l, n-discs, n-wedges, aperture.
-   *cylinder*: Closed 3D model of a cylinder oriented with the axis in
    the Z direction. Accepts r and l as parameters.
-   *circle*: Closed 2D model of a circle at the origin. Accepts the r
    parameter.
-   *free3d*: Open 3D model in which particles are initialized at
    `[0.0 0.0 0.0]` and move randomly in open space.
-   *sphere*: Closed 3D model of a sphere at the origin. Accepts the r
    parameter.
-   *segment*: Closed 1D model of a the segment `[-l/2, l/2]`. Accepts
    the l parameter.
-   *boxnd*: N-dimensional box. Accepts the parameters dim and size
    (Dimension 1 yields a segment, dim. 2 a rectangle and dim. 3 a
    cube).
-   *pointnd*: A variable dimension model that initializes at a point
    and keeps that position in all steps. Accepts the pos paramater in
    the form [coord1 coord2 ... coordN].

Available trackers
------------------

### The *spatial* tracker

The *spatial* tracker stores the sum of the absolute values of the
displacement vector (path-len) and its components (coord-len) and the
sum of the coordinates for each position (coord-sum). It does not use or
require any initialization parameter. This tracker works with models of
any dimension and its binary values are stored as 64bit floats with
structure "float64 [float64] [float64]": For each particle and step
there is one number for the path-len and then the numbers for coord-len
and coord-sum, as many as model dimensions for each. For example, for a
three dimensionas simulation the float64 numbers would be [path-len
coord-len-x coord-len-y coord-len-z coord-sum-x coord-sum-y coord-sum-z]
for each particle and stored step.

### The *diffmr1* tracker

The *diffmr1* tracker computes NMR phase of each particle during an MR
Diffusion experiment for three different bipolar gradient shapes (Dirac
deltas, square and sinusoidal) oriented in a single direction. The
direction can be chosen by the user passing the parameter
`":diffmr1.grad-dir [Gx Gy Gz]"` in the `--trackers-pars` option. This
tracker can work with models of any dimension, just pass the right
vector length when choosing a gradient direction. This tracker computes
the parameters called `:nmr-phase-delta`, `:nmr-phase-square` and
`:nmr-phase-sin`, which are stored as an array of 64bit floats for each
particle and storage step.

Available reducers
------------------

### The *pi* reducer

The *pi* reducer computes the value of the mathematical constant PI
after a simulation has been run with a square of size 1.0 as model. It
gets PI from the fraction of particles with initial position in the
square of side 1.0 that also happen to have initial position within a
circle of radius 0.5 contained inside the square. Since the reducer uses
only the initial positions, the Monte Carlo simulation can be run
without random walk by setting the steps parameter to [0].

### The *diffmr1adc* reducer

The *diffmr1adc* reducer computes Apparent Diffusion Coefficients after
a random walk simulation has been run with the *diffmr1* tracker
activated. The reducer will use the stored NMR phase values to compute
the "true" ADC coming from the distribution of NMR phases (\<phi\^2\> /
2b) and the ADC coming from a linear fit of the logarithms of the
attenuated signal magnitudes from different values of the diffusion
weighting parameter *b*. The optimal *b* range is calculated
automatically and signal magnitudes are computed for a set of values of
*b* distributed linearly along the range. Both the true and magnitude
ADCs are computed for the three gradient shapes (Dirac delta, square,
sinusoidal) of the *diffmr1* tracker.

Storage of simulation results
-----------------------------

At the end of the simulation several binary files will be stored, one
with the positions (positions.bin) that contains 64bit floating point
numbers, and one for each tracker with the tracker name in the filename
and containing the tracker values. The binary type for the tracker files
depends on the tracker, please consult tracker documentation. The order
of the data in the files is such that as you move through the file the
particle index changest the slowest, then the steps index, and finally
the coordinate index. The coordinate index takes as many values as the
dimension of the model used and the particle and steps indexes take a
number of values equal to the number of particles and number of storage
steps specified by the user.

Random numbers
--------------

The software uses the Java [Colt](http://acs.lbl.gov/software/colt/)
numerical library implementation of the excelent [Mersenne twister
PRNG](http://en.wikipedia.org/wiki/Mersenne_twister). The seeds are read
by default from "/dev/random" in Linux, in other operating systems the
software will attempt to read from a file called "seeds.bin" in the
working directory. In any OS a source binary file for random seeds can
be specified with the -r command line argument. Files containing high
quality random number can be downloaded from The
[Random.org](http://www.random.org/) website: Go to the ["Raw
Bytes"](http://www.random.org/bytes/) section, under the "Numbers" menu
entry, and after entering the number of bytes and selecting the
"Download to file" option and press the "Get Bytes" button. You will
need four times as many bytes as the number of workers the program will
launch (double that if using different initialization and simulation
models); the number of workers is (n-cpu + 2) by default. The software
also uses the Incanter Clojure Library for some operations (ADC
calculation in the diffmr1adc reducer). The intention is to replace Colt
with Incanter for random number generation in the future, but the
possibility of manually seeding the PRNG is currently missing from the
Incanter statistical distributions.

Projected features
------------------

The following features and improvements are projected for future
versions of this software:

-   New models: **acinartree**.
-   Modules loadable at run time without re-compilation of the whole
    software.
-   A *diffmrrand* tracker with random gradient direction for each
    particle.
-   A *diffmrn* tracker supporting several directions in a single
    simulation (for DTI simulations).
-   Several diffusion encoding times in a single simulation for the
    *diffmrX* trackers.
-   Possibility of reading data files for positions or tracker values,
    and running only certain parts of the software (example: read
    positions, perform track and reduce, read positions and tracker
    values and apply reducers).
-   Better user and developer documentation.
-   Improved error handling.
-   Storage of reducers results.

Acknowledgments
---------------

The initial version of this project had about 2200 lines of code (comments and 
blank lines excluded). It was created in the period 2013-2014 by Jose M. 
Perez-Sanchez (JMPS) during his participation in the Cardio-Image program, 
funded by the Spanish National Center for Cardiovascular Research (CNIC) in 
Madrid. 

During his participation in the Cardio-Image program the author held also a 
position at the Spanish Network for Biomedical Research in Respiratory Diseases 
(CIBERES) and worked at the Mount Sinai Medical Center in New York City. The 
author is very thankful to these three institutions.

This software was a rewrite of a previous C++ version that had gotten too
difficult to maintain. Some of the initial models and many ideas behind
this software were conceived when the C++ version of GEMA was developed. 
Ignacio Rodriguez was one of the original authors of the C++ version of GEMA,
together with JMPS. The author is very thankful to him for his support and help 
during the creation of this software.

Copyright
---------

The copyright of this software is owned by Centro Nacional de Investigaciones
Cardiovasculares (http://www.cnic.es) and CIBERES (http://www.ciberes.org).

License
-------

This software is licensed under the Eclipse Public License version 1.0. Please,
check the *LICENSE.txt* file for a copy of the license.

