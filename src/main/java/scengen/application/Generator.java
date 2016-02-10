package main.java.scengen.application;


import java.util.HashMap;
import java.util.Map;

import main.java.scengen.distribution.MultivariateDistribution;
import main.java.scengen.distribution.MultivariateLognormal;
import main.java.scengen.distribution.MultivariateNormal;
import main.java.scengen.distribution.MultivariateStudent;
import main.java.scengen.distribution.MultivariateUniform;
import main.java.scengen.methods.METHOD;
import main.java.scengen.methods.MomentMatching;
import main.java.scengen.methods.MonteCarlo;
import main.java.scengen.methods.QuantizationGrids;
import main.java.scengen.methods.QuantizationLearning;
import main.java.scengen.methods.QuasiMonteCarlo;
import main.java.scengen.methods.ReductionMethod;
import main.java.scengen.methods.Scenred2;
import main.java.scengen.methods.VoronoiCellSampling;
import main.java.scengen.tools.Xorshift;

/**
 * <p>The Generator class implements a simple interface to generate a reduced set of scenarios from different multivariate probability distributions.</p>
 * @author Nils Loehndorf
 *
 */
public class Generator {
	
	static long _seed = 191;
	
	/**
	 * <p>Generate a reduced set of scenarios from a normal, log-normal, or uniform, probability distribution with given parameters.</br>
	 * Use the following arguments to run the program:</p>
	 * <ul><li>dist=(Normal, Lognormal, Student, Uniform)</li>
	 * <li>dim=int</li>
	 * <li>mean=double,...,double</li>
	 * <li>cov=double,...double</li>
	 * <li>scen=int</li>
	 * <li>seed=int</li>
	 * <li>alt uniform: min=double,...,double max=double,...,double correl=double,...double</li>
	 * <li>alt student: mean=double,...,double scale=double,...,double df=double,...,double correl=double,...double</li>
	 * <li>optional: method=(MonteCarlo, QuasiMonteCarlo, MomentMatching, QuantizationGrids, QuantizationLearning, VoronoiCellSampling)</li>
	 * <p>Example usage to generate 20 scenarios from a 2-dimensional normal distribution with mean (5,10) and 
	 * covariance ((50,25),(25,100)) using quasi Quasi-Monte Carlo as scenario reduction method:</br>
	 * 'java -jar scengen.jar dist=Normal dim=2 mean=5,10 cov=50,25,25,100 scen=20 method=QuasiMonteCarlo'</p>
	 * @param args 
	 */
	public static void main (String... args) {
		_seed = new Xorshift().nextLong();
		Map<String,String> table = makeTable(args);
		if (!table.containsKey("dist"))
			throw new IllegalArgumentException("Use keyword 'dist=Lognormal,Normal,Student,Uniform' to specify distribution.");
		if (!table.containsKey("dim"))
			throw new IllegalArgumentException("Use keyword 'dim=int' to define distribution dimension.");
		if (!table.containsKey("scen"))
			throw new IllegalArgumentException("Use keyword 'scen=int' to specify the number of scenarios.");
		if (!((table.containsKey("mean") && table.containsKey("cov")) || (table.get("dist").matches("Uniform")  && table.containsKey("min") && table.containsKey("max") && table.containsKey("correl")) || (table.get("dist").matches("Student")  && table.containsKey("mean") && table.containsKey("scale") && table.containsKey("correl"))))
			throw new IllegalArgumentException("Use keywords 'mean=double,...double cov=double,...double' to specify distribution parameters.");
		String distName = table.get("dist");
		int dim = Integer.parseInt(table.get("dim"));
		long seed = new Xorshift().nextLong();
		if (table.containsKey("seed"))
			seed = Long.parseLong(table.get("seed"));
		if (table.containsKey("seed"))
			seed = Integer.parseInt(table.get("seed"));
		MultivariateNormal dist = null;
		switch (distName) {
		case "Normal": {
				String[] smean = table.get("mean").split(",");
				if (smean.length!=dim)
					throw new ArrayIndexOutOfBoundsException("Size of a parameter vector does not match distribution dimension.");
				double[] mean = new double[dim];
				for (int i=0; i<dim; i++)
					mean[i] = Double.parseDouble(smean[i]);
				String[] scov = table.get("cov").split(",");
				if (scov.length!=dim*dim)
					throw new ArrayIndexOutOfBoundsException("Size of a parameter vector does not match distribution dimension.");
				double[][]cov = new double[dim][dim];
				for (int i=0; i<dim; i++)
					for (int j=0; j<dim; j++)
						cov[i][j] = Double.parseDouble(scov[j+i*dim]);
				dist = new MultivariateNormal(mean,cov,new Xorshift(seed));}
				break;
		case "Lognormal": {
			String[] smean = table.get("mean").split(",");
			if (smean.length!=dim)
				throw new ArrayIndexOutOfBoundsException("Size of a parameter vector does not match distribution dimension.");
			double[] mean = new double[dim];
			for (int i=0; i<dim; i++)
				mean[i] = Double.parseDouble(smean[i]);
			String[] scov = table.get("cov").split(",");
			if (scov.length!=dim*dim)
				throw new ArrayIndexOutOfBoundsException("Size of a parameter vector does not match distribution dimension.");
			double[][]cov = new double[dim][dim];
			for (int i=0; i<dim; i++)
				for (int j=0; j<dim; j++)
					cov[i][j] = Double.parseDouble(scov[j+i*dim]);
			dist = new MultivariateLognormal(mean,cov,new Xorshift(seed));}
			break;
		case "Student": {
			if (table.containsKey("mean") && table.containsKey("cov")) {
				String[] smean = table.get("mean").split(",");
				if (smean.length!=dim)
					throw new ArrayIndexOutOfBoundsException("Size of a parameter vector does not match distribution dimension.");
				double[] mean = new double[dim];
				for (int i=0; i<dim; i++)
					mean[i] = Double.parseDouble(smean[i]);
				int[] df = new int[dim];
				for (int i=0; i<dim; i++)
					df[i] = 5;
				String[] scov = table.get("cov").split(",");
				if (scov.length!=dim*dim)
					throw new ArrayIndexOutOfBoundsException("Size of a parameter vector does not match distribution dimension.");
				double[][] cov = new double[dim][dim];
				for (int i=0; i<dim; i++)
					for (int j=0; j<dim; j++)
						cov[i][j] = Double.parseDouble(scov[j+i*dim]);
				double[] scale = new double[dim];
				for (int i=0; i<dim; i++) 
					scale[i] = Math.sqrt(cov[i][i]*3./5.);
				dist = new MultivariateStudent(mean,scale,getCorrelationMatrix(cov),df,new Xorshift(seed));
			}
			else {
				String[] smean = table.get("mean").split(",");
				if (smean.length!=dim)
					throw new ArrayIndexOutOfBoundsException("Size of a parameter vector does not match distribution dimension.");
				double[] mean = new double[dim];
				for (int i=0; i<dim; i++)
					mean[i] = Double.parseDouble(smean[i]);
				String[] sdf = table.get("df").split(",");
				if (smean.length!=dim)
					throw new ArrayIndexOutOfBoundsException("Size of a parameter vector does not match distribution dimension.");
				int[] df = new int[dim];
				for (int i=0; i<dim; i++)
					df[i] = Integer.parseInt(sdf[i]);
				String[] scov = table.get("correl").split(",");
				if (scov.length!=dim*dim)
					throw new ArrayIndexOutOfBoundsException("Size of a parameter vector does not match distribution dimension.");
				double[][]cor = new double[dim][dim];
				for (int i=0; i<dim; i++)
					for (int j=0; j<dim; j++)
						cor[i][j] = Double.parseDouble(scov[j+i*dim]);
				String[] sscale = table.get("scale").split(",");
				if (sscale.length!=dim)
					throw new ArrayIndexOutOfBoundsException("Size of a parameter vector does not match distribution dimension.");
				double[] scale = new double[dim];
				for (int i=0; i<dim; i++)
					scale[i] = Double.parseDouble(sscale[i]);
				dist = new MultivariateStudent(mean,scale,cor,df,new Xorshift(seed));
			}}
			break;
		case "Uniform": {
			if (table.containsKey("mean") && table.containsKey("cov")) {
				String[] smean = table.get("mean").split(",");
				if (smean.length!=dim)
					throw new ArrayIndexOutOfBoundsException("Size of a parameter vector does not match distribution dimension.");
				double[] mean = new double[dim];
				for (int i=0; i<dim; i++)
					mean[i] = Double.parseDouble(smean[i]);
				String[] scov = table.get("cov").split(",");
				if (scov.length!=dim*dim)
					throw new ArrayIndexOutOfBoundsException("Size of a parameter vector does not match distribution dimension.");
				double[][]cov = new double[dim][dim];
				for (int i=0; i<dim; i++)
					for (int j=0; j<dim; j++)
						cov[i][j] = Double.parseDouble(scov[j+i*dim]);
				dist = new MultivariateUniform(mean,cov,new Xorshift(seed));
			}
			else {
				String[] smin = table.get("min").split(",");
				if (smin.length!=dim)
					throw new ArrayIndexOutOfBoundsException("Size of a parameter vector does not match distribution dimension.");
				double[] min = new double[dim];
				for (int i=0; i<dim; i++)
					min[i] = Double.parseDouble(smin[i]);
				String[] smax = table.get("max").split(",");
				if (smax.length!=dim)
					throw new ArrayIndexOutOfBoundsException("Size of a parameter vector does not match distribution dimension.");
				double[] max = new double[dim];
				for (int i=0; i<dim; i++)
					max[i] = Double.parseDouble(smax[i]);
				String[] scorrel = table.get("correl").split(",");
				if (scorrel.length!=dim*dim)
					throw new ArrayIndexOutOfBoundsException("Size of a parameter vector does not match distribution dimension.");
				double[][]correl = new double[dim][dim];
				for (int i=0; i<dim; i++) {
					for (int j=0; j<dim; j++) {
						correl[i][j] = Double.parseDouble(scorrel[j+i*dim]);
						if (correl[i][j]>1 || correl[i][j]<-1)
							throw new IllegalArgumentException("Correlation matrix for uniform distribution invalid.");
					}
					if (correl[i][i]!=1)
						throw new IllegalArgumentException("Not all diagonal entries of correlation matrix are equal to 1.");
				}
				dist = new MultivariateUniform(min,max,correl,new Xorshift(seed));
			}}
			break;
		default: 
			throw new IllegalArgumentException(String.format("Distribution %s unknown. Use 'Normal', 'Lognormal', or 'Uniform'.",args[0]));
		}
		int numScen = Integer.parseInt(table.get("scen"));
		//optional parameters
		METHOD method = null;
		if (table.containsKey("method"))
			method = METHOD.valueOf(table.get("method"));
		else {
			if (dim<=2) 
				method = METHOD.QuantizationLearning;
			else 
				method = METHOD.VoronoiCellSampling;
		}
		Map<double[],Double> scen = runScenarioGenerator(numScen,method,dist,seed);
		printf("Multivariate distribution:\t%s \n",dist.getType());
		printf("Random vector dimenensionality:\t%s \n",dim);
		printf("Number of scenarios:\t%s \n",numScen);
		printf("Scenario generation method:\t%s \n",method);
		printf("Generated scenarios by row:\n");
		printf("probability");
		for (int i=0; i<dim; i++)
			printf(String.format("\tvalue(%d)",i));
		printf("\n");
		for (double[] x : scen.keySet()) {
			printf("%f",scen.get(x));
			for (int i=0; i<dim; i++)
				printf("\t%f",x[i]);
			printf("\n");
		}
	}
	
	static Map<String,String> makeTable(String... args) {
		Map<String,String> map = new HashMap<>();
		for (String s : args) {
			String[] ar = s.split("=");
			map.put(ar[0],ar[1]);
		}
		return map;
	}
	
	static void printf(String s, Object... args) {
		System.out.printf(s,args);
	}
	
	public static Map<double[],Double> runScenarioGenerator(int numScen, METHOD method, MultivariateDistribution mvDist, long seed) {
		ReductionMethod scenred = null;
		switch (method) {
		case MomentMatching: 
			scenred = new MomentMatching(mvDist.getMean(),mvDist.getCov(),mvDist.getSkewness(),mvDist.getKurtosis(),seed);
			break;
		case QuasiMonteCarlo: 
			scenred = new QuasiMonteCarlo(mvDist);
			break;
		case QuantizationGrids: 
			scenred = new QuantizationGrids(mvDist);
			break;
		case QuantizationLearning: 
			scenred = new QuantizationLearning(mvDist);
			break;
		case MonteCarlo: 
			scenred = new MonteCarlo(mvDist);
			break;
		case Scenred2:
			scenred = new Scenred2(new MonteCarlo(mvDist).getScenarios(10000),mvDist.getDim());
			break;
		case VoronoiCellSampling:
			scenred = new VoronoiCellSampling(mvDist);
			break;
		}
		return scenred.getScenarios(numScen);
	}
	
	public static double[][] getCorrelationMatrix(double[][] covariance) {
		if (covariance.length!=covariance[0].length)
			throw new IllegalArgumentException("No sqaure matrix.");
		int dim = covariance.length;
		double[][] correlation = new double[dim][dim];
		for (int i=0; i<dim; i++) {
			correlation[i][i] = 1;
			for (int j=i+1; j<dim; j++) {
				correlation[i][j] = covariance[i][j]/Math.sqrt(covariance[i][i]*covariance[j][j]);
				correlation[j][i] = correlation[i][j];
			}
		}
		return correlation;
	}
 	
	

	
}
