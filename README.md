Probabilistic latent semantic analysis (PLSA, multiprocessor version)
=====================================================================


Introduction
------------

Probabilistic latent semantic analysis is a method for analyzing co-occurrence data.  This source code is a parallelized version of PLSA which makes use of OpenMP for shared memory parallelization and the Message Passing Interface (MPI) for distributed memory parallelization.  This work is described by the poster paper:

    R. Wan, V. N. Anh, and H. Mamitsuka. Efficient Probabilistic Latent Semantic Analysis Through Parallelization. In Proc. 5th Asia Information Retrieval Symposium (poster session), volume 5839 of LNCS, pages 432-443, October 2009.

In that paper, there are five systems.  This source code generates four of them:  None, MPI, OpenMP, and All.  In the description below, we refer to all four systems simply as "PLSA".


Compiling
---------

The PLSA software is written in C and has been compiled using the GNU gcc compiler v4.1.2 under Linux.  OpenMP is included with the compiler and Open MPI v1.2.4 was used.  The system was tested on a 64-bit system, but does not make use of any features specific to 64-bit architectures.

The PLSA system uses Autotools in order to create a Makefile which is configured to your system's settings.  This is done using the "configure" script.  Simply perform the following:

    ./configure;make


Running PLSA
------------

Run `plsa` without any arguments to get a list of possible options:

    Probabilistic Latent Semantic Analysis (multiprocessor)
    =======================================================

    Usage:  ./plsa [options]

    Options:
    --base <file>      :  Base filename for output file.
    --cooccur <file>   :  Co-occurrence filename.
    --clusters <int>   :  Number of clusters.
    --seed <int>       :  Random seed.
                       :    (Default:  current time).
    --maxiter <int>    :  Maximum iterations.
    --text             :  Text mode (I/O is in text, not binary).
    --snapshot <int>   :  Output snapshots p(x,y) at regular intervals.
                       :    (Default:  Do not output).
    --openmp <int>     :  Number of OpenMP threads to use
                            (Default:  Maximum for PC).
    --verbose          :  Verbose mode.
    --debug            :  Debugging output.
    --rounding         :  Round using 100000000 as the multiplication factor.
    --nooutput         :  Suppress outputting p(x,y) to file.

    Compile-time settings:
           MPI:                              Enabled
        OpenMP:                              Enabled

    PLSA version:  Mar  7 2010 (15:10:18)


    
From top-to-bottom, the options mean:

* --base:      The filename, before the extension, of the output file.  The extension is fixed as ".plsa".
* --cooccur:   The input co-occurrence file, whose format is described below.
* --clusters:  The number of latent states.
* --seed:      The random seed to use.  If none is provided, the current system time is used.
* --maxiter:   The maximum number of iterations of the EM algorithm to perform.  One of two stopping criteria.
* --text:      Indicate that the input file is in text and not binary; useful for debugging.
* --snapshot:  Output snapshots of p(x,y) at certain intervals.  Useful if PLSA is taking a long time to run and intermediate results are required.
* --openmp:    The number of threads of execution to use for OpenMP.
* --verbose:   Verbose output.
* --debug:     Debugging output.  Output is generated as each value is read from the input file.  (Note that a lot of output will be generated.)
* --rounding:  Round the output values in p(x,y) using the specified rounding factor.  That is, if the factor is "1000", then three decimal places are used.  Useful for comparing methods due to the problem with floating point arithmetic (details below).
* --nooutput:  Do not produce the final output file.  Eliminates the creation of a fairly large file.

Many of these parameters have no defaults (such as `--maxiter` and  `--clusters`), so they will have to be explicitly given.
    

Data format
-----------

The input file is provided to the system using the --cooccur switch.  It can be in one of two formats:  text or binary.  Binary mode is the default and suggested for large files to save space; text mode is good for testing.

In either mode, the format of the file is:

    [rows][columns][row id+][column id+][w1 cos_count (w21 c21) ... (w2n c2n)]+

which means:

- rows --  The number of rows in the matrix (unsigned int)
- columns --  The number of columns in the matrix (unsigned int)
- row id --  A sequence of row id's.  The number of them should be exactly "rows".
- column id --  A sequence of column id's.  The number of them should be exactly "columns".
- [w1 cos_count (w21 c21) ... (w2n c2n)] -- A row in the matrix which
    describes only the non-zero positions.
    - w1 --  The ID of the row.
    - cos_count --  The number of non-zero values in the row.
    - w21 -- Column ID of the value.
    - c21 -- The value itself.

All values are unsigned integers.  The program does not make use of either "row id" or "column id" (it was included for a future feature of the program that has not yet been implemented).  Any value is fine (i.e., 0 or sequential integers).

In text format, the values are separated by a single whitespace (usually the tab character).  In binary mode,  they are unsigned integers (usually 4 bytes in size each).

Remember to specify whether it is a text or binary file by using (or not using) the `--text` switch.  The program does NOT check for binary or text mode and problems will occur if the mode does not match the file.

Please see the source in input.c for further details on the file format.


Sample run
----------

Suppose we have a 3 by 4 matrix like this:

    0     0     3     4
    0     0     0     6
    7     2     1     0

then the test file will appear like this:

    3     4
    0     1     2
    0     1     2     3
    0     2     2     3     3     4
    1     1     3     6
    2     3     0     7     1     2     2     1


The first 3 rows are self-explanatory.  The fourth row is interpretted as follows.  In row 0, there are 2 non-zero co-occurrences.  They are at position 2 (with a value of 3) and position 3 (with a value of 4).  Therefore, we have:  (0, 2) = 3 and (0, 3) = 4.

We execute `plsa` as follows:

    ./plsa --cooccur test.cooccur --maxiter 30 --clusters 4 --text --verbose --debug --base out

Without `--verbose` and `--debug`, there would be no output generated.  Output will vary from system to system, but some sample output would look like the following:


    Settings
    --------
    ==  Base filename:                                  out
    ==  Co-occurrence filename:                         test.cooccur
    ==  Probability data type:                          double
    ==  Clusters:                                       4
    ==  Random seed:                                    [from time]
    ==  Exponent difference [utils.h::addLogsFloat]:    23.02585093
    ==  Termination conditions
    ==    Maximum EM iterations:                        30
    ==    Percentage difference:                        0.001000
    ==  Text mode:                                      yes
    ==  Rounding:                                       no
    ==  Suppress output to file:                        no
    ==  MPI:                                            OK
    ==    My ID:                                        0
    ==    Number of processes:                          1
    ==    Block range:                                  0 - 3
    ==    Block size:                                   4
    ==  Open MP:                                        OK
    ==    Number of threads:                            4


    ==  Reading from co-occurrence file...
    ==  Applying seed from time:                        1267944777
    ==    Read (0, 2) --> 3
    ==    Read (0, 3) --> 4
    ==    Read (1, 3) --> 6
    ==    Read (2, 0) --> 7
    ==    Read (2, 1) --> 2
    ==    Read (2, 2) --> 1
    ==  ID 0 finished reading co-occurrence data.
    ==  Maximum number of pairs:                        12
    ==  Actual number of pairs in data file:            6
    ==  Percentage of zeroes:                           50.00 % (6)
    ==  Sum of co-occurrence counts:                    23
    ==  Begin initialization...
    ==  Initialization complete...
    ==  m = 3; n = 4
    [---]  Initial = -59.143897
    [  1]  -59.143897 --> -51.012234  [8.131663, 13.7489 %]
    [  2]  -51.012234 --> -46.918000  [4.094234, 8.0260 %]
    [  3]  -46.918000 --> -41.866035  [5.051965, 10.7676 %]
    [  4]  -41.866035 --> -38.945650  [2.920385, 6.9755 %]
    [  5]  -38.945650 --> -38.070771  [0.874878, 2.2464 %]
    [  6]  -38.070771 --> -37.746852  [0.323920, 0.8508 %]
    [  7]  -37.746852 --> -37.614179  [0.132673, 0.3515 %]
    [  8]  -37.614179 --> -37.561918  [0.052260, 0.1389 %]
    [  9]  -37.561918 --> -37.540558  [0.021360, 0.0569 %]
    [ 10]  -37.540558 --> -37.530982  [0.009577, 0.0255 %]
    [ 11]  -37.530982 --> -37.526141  [0.004841, 0.0129 %]
    [ 12]  -37.526141 --> -37.523373  [0.002768, 0.0074 %]
    [ 13]  -37.523373 --> -37.521621  [0.001752, 0.0047 %]
    [ 14]  -37.521621 --> -37.520430  [0.001191, 0.0032 %]
    [ 15]  -37.520430 --> -37.519583  [0.000847, 0.0023 %]
    [ 16]  -37.519583 --> -37.518963  [0.000620, 0.0017 %]
    [ 17]  -37.518963 --> -37.518506  [0.000456, 0.0012 %]
    [ 18]  -37.518506 --> -37.518162  [0.000344, 0.0009 %]
    ==  Non-probabilities:                              0
    ==  Sum of p(x,y):                                  1.000000
    ==  Total output files printed                      1
    ==  Program execution:                              0.000 secs


The resulting matrix of p(x,y) is stored as out.plsa and look like the following:

    -----
    3     4     0     1     2     0     1     2     3     
    -24.309359     -22.912351     -2.036586     -1.749422     
    -96.419985     -95.022977     -12.677112     -1.343747     
    -1.189584     -2.442347     -3.136454     -10.084368
    -----

Whose output format is the same as the input format, except that the integral co-occurrence counts are replaced with probabilities in log-space as floating point values.


Other issues
------------

This section lists some issues relevant to running or extending the program, in no particular order.

Applicable to both versions:

1.  All values are stored in log-space to prevent underflow of floating point values.

2.  Probabilities are multiplied by adding values in log-space.  To prevent underflow, a macro called `logSumsInline` is used (defined in `plsa-defn.h`).  It is provided a cut-off called `LN_LIMIT` which stipulates what cut-off (as a ln value) to use when adding two numbers.  If the smallest number is "too small", then no addition is performed and the larger value is taken.  As a result, the order in which numbers are added matters and if   MPI is used a small, negligible difference will appear.  In other words, if

    x = a + b + c + d

then when MPI with two processors is used, then

    x = y + w
    y = a + b
    w = c + d

If c is very small compared to (a + b) but not so when compared to d, then in the first case, it might be dropped.  Generally, the percentage difference in maximum likelihood will not change enough to matter.

3.  Also related to the previous point, p(x,y) values written to the 
output file may be different if the binary files are compared directly due 
to the number of digits of precision.  If sent to a text file with only 6 
digits of precision, there may not be any problems.  To solve this, the 
`--rounding` option was added.  It and the cut-off defined in plsa-defn.h 
will need to be tweaked if two runs are to be compared.

4.  Generating the output file could take half of the total execution time, depending on the matrix size.  To suppress printing this file, use the `--nooutput` switch.

5.  To get the time required for a single iteration of the loop (as reported in the paper cited in Section 1), use the `--maxiter 1` option.

6.  The four random seeds used for the paper cited in Section 1 were:  "20444 3612 31325 17062".  Of course, this information alone is not enough to generate the exact same results since every system's random number generator might be slightly different.

7.  The number of latent states must be larger than the number of processors under MPI.  If this is not the case, then the number of latent states is increased automatically.


Applicable to this version only:

1.  The maximum number of latent states is defined by the length of the message tag (for MPI).  By default, the maximum is 999.  To change this, edit the value for MAX_CLUSTERS in plsa-defn.h.

2.  To disable either MPI or OpenMP, run "./configure" and before compiling with "make", edit the file "config.h".  Search for the values for HAVE_MPI or HAVE_OPENMP and change them to 0.  Also, if either of them is 0, that means that library (MPI or OpenMP) was not detected by configure.


About PLSA
----------

This software was implemented while I was at Kyoto University (around 2009).  My contact details:

     E-mail:  rwan.work@gmail.com 

My homepage is [here](http://www.rwanwork.info/).

The latest version of PLSA can be downloaded from [GitHub](https://github.com/rwanwork/PLSA-Base).

If you have any information about bugs, suggestions for the documentation or just have some general comments, feel free to contact me via e-mail or GitHub.


Copyright and license
---------------------

    Probablistic latent semantic analysis (PLSA, multiprocessor version)
    Copyright (C) 2009-2010 by Raymond Wan

This software is distributed under the terms of the GNU General Public License (GPL, version 3 or later) -- see the file LICENSE for details.

Permission is granted to copy, distribute and/or modify this document under the terms of the GNU Free Documentation License, Version 1.3 or any later version published by the Free Software Foundation; with no Invariant Sections, no Front-Cover Texts and no Back-Cover Texts. A copy of the license is included with the archive as LICENSE.

