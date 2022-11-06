# SCENARIO GENERATION FOR STOCHASTIC PROGRAMMING #

An important aspect of stochastic optimization is finding a set of scenarios of the random variables that is representative but still enumerable.
This project aims at bringing the most important methods for scenario generation together to make them accessible as well as comparable. The project
comes with a standardized interface for different methods and distributions, as well as with a number of test problems to study the performance
of the methods for different stochastic optimization problems.

The methods themselves as well as an empirical study that compares their performance is summarized in this [working paper](http://www.optimization-online.org/DB_HTML/2016/02/5319.html). The numerical section of the paper is based on the numerical experiments that can be found in [src/main/java/scengen/paper/](src/main/java/scengen/paper/).

### Binary ###
The library is packaged as an [executable jar](http://www.loehndorf.com/nils/scengen/scengen.jar) that can be run from the command line. 

Examples:
* Generate 7 scenarios from a normal distribution with mean 5 and standard deviation 2 using the quantization grids from the [optimal quantization website](http://www.quantize.maths-fi.com): </br>
`java -jar scengen.jar dist=Normal dim=1 mean=10 cov=5 scen=7 method=QuantizationGrids`
* Get 20 scenarios from a bivariate normal distribution with mean (1,1) and covariance (1,0.5,0.5,1) using quasi-Monte Carlo: </br>
`java -jar scengen.jar dist=Normal dim=2 mean=(1,1) cov=(1,0.5,0.5,1) scen=20 method=QuasiMonteCarlo`
* Get 10 scenarios from a bivariate log-normal distribution with mean (1,1) and covariance (1,0.5,0.5,1) using moment matching: </br>
`java -jar scengen.jar dist=Lognormal dim=2 mean=(1,1) cov=(1,0.5,0.5,1) scen=20 method=MomentMatching`

Arguments:
* dist=(Normal, Lognormal, Student, Uniform)
* dim=int
* mean=double,...,double
* cov=double,...,double (flat array of the full matrix)
* scen=int

Optional Arguments:
* min=double,...,double (required for Uniform)
* max=double,...,double (required for Uniform)
* correl=double,...double (flat array of the full matrix, required for Uniform and Student)
* scale=double,...,double (required for Student)
* df=double,...,double (required for Student)
* method=(MonteCarlo, QuasiMonteCarlo, MomentMatching, QuantizationGrids, QuantizationLearning, VoronoiCellSampling)
* seed=int

The required executables for the following two methods must be available in a subfolder bin/
* [Moment matching](http://work.michalkaut.net/downloads.html) by Michal Kaut
* [Scenred2](https://www.gams.com/help/index.jsp?topic=%2Fgams.doc%2Ftools%2Fscenred2%2Findex.html) by Holger Heitsch / Werner Römisch

### Source ###
1. Clone the repository and import it into an IDE, e.g., Eclipse.
2. Add lib/ to build path.
3. Add src/main/resources/vc/ to class path.
