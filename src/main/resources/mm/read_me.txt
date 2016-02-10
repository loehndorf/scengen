Scenario-generation code scen-gen_HKW

This code implements the moment-matching heuristic from paper 'A Heuristic for
Moment-matching Scenario Generation' by Kjetil Høyland, Michal Kaut and Stein W.
Wallace, published in Computational Optimization and Applications, 24 (2-3), pp.
169–185, 2003; doi:10.1023/A:1021853807313.
Most of the code was written by Michal Kaut, except for the cubic_solve()
function written by Diego Mathieu.
The code has been tested with gcc on both Linux and Windows.
Problems, questions, and suggestions should be send to: michal.kaut@iot.ntnu.no

LICENSE
The code is distributed under the Eclipse Public License, see
http://www.eclipse.org/legal/epl-v10.html. For information about the license,
including compatibility with other licenses, see the official FAQ at
http://www.eclipse.org/legal/eplfaq.php or its Wikipedia entry at
http://en.wikipedia.org/wiki/Eclipse_Public_License. 


ALGORITHM
The algorithm uses a moment-matching approach, i.e. the user specifies first 
four moments and a correlation matrix for the random variables, and the 
algorithm creates a scenario tree that matches the given properties.

The algorithm needs fix pre-specified probabilities. By default, it uses 
equiprobable scenarios, but it is also possible to supply the probabilities in 
an additional file.


COMPILING
The package includes a project file scengen_HKW for the open-source IDE
Code::Blocks, see http://www.codeblocks.org/. This project allows to build
the stand-alone executable and also a DLL with the main procedure, plus
a small executable that uses this DLL.

In addition, there is a very simple Makefile to build the main executable.


USAGE
The file expects two files to exist in the same directory as the executable:
tg_moms.txt  - target moments
tg_corrs.txt - target correlation matrix
Both the files must be in a special "matrix" format:
  number of rows
  number of columns
  data
Rest of the file is ignored.
For an example, see the sample files in doc/sample_inputs.7z

To find the syntax of the file, run it without any arguments.
For more information, run it with a single parameter -h.

The code can accept moments in several different formats. The default is:
 - 1 = mean
 - 2 = sample standard deviation - this is the biased version, where the sum
       of deviations is divided by N, instead of (N-1).
       This is used because it can be generalized for non-equiprobable case
       (The default means we multiply all values by probability 1/N.)
 - k>2 = E{[X - mean]^k} / stDev^k
If one prefers to specify moments in other formats, there is the "-f" option,
where the value is a sum of following bits:
 -  1 -> use population estimators (as in spreadsheets)
 -  2 -> 2nd moment is Var instead of StDev
 -  4 -> 4th moment is Kurtose - 3
 -  8 -> Higher moments are not scaled by StDev
 - 16 -> TarMom[i-1] = E{X^i} ... lower bits are ignored in this case
For example, to use spreadsheet functions (average, stdev, skew, kurt), one
needs to use "-f 1". To use variance (var) instead of stdev, one would need
to use "-f 3", since we use the first two bits and 1 + 2 = 3.


DOCUMENTATION
To generate a full documentation, run doxygen in the main folder.
This will create HTML documentation in the doc/html folder.