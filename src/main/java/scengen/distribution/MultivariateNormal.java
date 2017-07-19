package scengen.distribution;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;

import umontreal.ssj.probdist.NormalDistQuick;
import scengen.methods.MonteCarlo;
import scengen.methods.ReductionMethod;

public class MultivariateNormal implements MultivariateDistribution {
	
	double[][] _chol;
	double[][] _cov;
	double[] _var;
	double[] _mean;
	int _dim;
	Random _generator;
	
	public MultivariateNormal(double[] mean, double[][] covariance, Random generator) {
		if (covariance.length != mean.length || covariance[0].length != mean.length)
	         throw new IllegalArgumentException("Covariance matrix dimensions invalid.");
		_mean = mean;
		_chol = choleskyDecomposition(covariance);
		_cov = covariance;
		_dim = mean.length;
		_var = new double[_dim];
		for (int i=0; i<_dim; i++)
			_var[i] = _cov[i][i];
		_generator = generator;
	}
	
	public static void main (String... args) {
		int dim = 1;
		double[] mean = new double[dim];
		double[][] correl = new double[dim][dim];
		for (int i=0; i<dim; i++) {
			mean[i] = 10;
			correl[i][i] = 2.5*2.5;
		}
		ReductionMethod mm = new MonteCarlo(new MultivariateNormal(mean,correl,new Random()));
		Map<double[],Double> map = mm.getScenarios(20);
		for (double[] x : map.keySet()) {
			for (int i=0; i<x.length; i++)
				System.out.print(x[i]+"\t");
			System.out.println();
		}
	}
	
	
	public void setSeed(long n) {
		_generator.setSeed(n);
	}
	
	//cholesky-banachiewicz algorithm
	static double[][] choleskyDecomposition(double[][] a){
		int m = a.length;
		double[][] l = new double[m][m]; //automatically initialzed to 0's
		for(int i = 0; i< m; i++){
			for(int k = 0; k < (i+1); k++){
				double sum = 0;
				for(int j = 0; j < k; j++){
					sum += l[i][j] * l[k][j];
				}
				l[i][k] = (i == k) ? Math.sqrt(a[i][i] - sum) :
					(1.0 / l[k][k] * (a[i][k] - sum));
			}
		}
		return l;
	}
	
	public double[] getRealization() {
		double[] z = new double[_dim];
		for (int i=0; i<_dim; i++)
			z[i] = _generator.nextGaussian();
		return transform(z);
	}
	
	public double[] transform(double[] z) {
		double[] x = new double[_chol.length];
		for (int i=0; i<_dim; i++) {
			double sum = 0.0;
			for (int j=0; j<i+1; j++)
				sum += z[j]*_chol[i][j];
			x[i] = _mean[i] + sum;
		}
		return x;
	}

	@Override
	public int getDim() {
		return _dim;
	}

	@Override
	public List<double[]> getSample(int size) {
		List<double[]> list = new LinkedList<>();
		for (int i=0; i<size; i++)
			list.add(getRealization());
		return list;
	}

	@Override
	public void setGenerator(Random generator) {
		_generator = generator;
	}

	@Override
	public Random getGenerator() {
		return _generator;
	}

	@Override
	public double[] getMean() {
		return _mean;
	}

	@Override
	public double[][] getCov() {
		return _cov;
	}
	
	public double[][] getCholesky() {
		return _chol;
	}

	@Override
	public double[] getSkewness() {
		return new double[_dim];
	}

	@Override
	public double[] getKurtosis() {
		return new double[_dim];
	}

	@Override
	public DISTRIBUTION getType() {
		return DISTRIBUTION.Normal;
	}

	@Override
	public double[] getRealization(double[] u) {
		double[] z = new double[_dim];
		for (int i=0; i<_dim; i++)
			z[i] = NormalDistQuick.inverseF01(u[i]);
		return transform(z);
	}

	@Override
	public double[] getVariance() {
		return _var;
	}
	
	public boolean isIndependent() {
		for (int i=0; i<_dim; i++) {
			for (int j=i+1; j<_dim; j++) {
				if (_cov[i][j]!=0)
					return false;
			}
		}
		return true;
	}

}
