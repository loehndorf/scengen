package scengen.distribution;

import java.util.Arrays;
import java.util.Random;

import umontreal.ssj.probdist.NormalDistQuick;

public class MultivariateUniform extends MultivariateNormal {
	
	double[] _min;
	double[] _max;
	double[] _correl2;
	
	public MultivariateUniform(double[] min, double[] max, double[][] correlation, Random generator) {
		super(new double[min.length],adjustCorrelation(correlation),generator);
		_min = min;
		_max = max;
		_mean = new double[_dim];
		_cov = new double[_dim][_dim];
		for (int i=0; i<_dim; i++) {
			_mean[i] = min[i] + (max[i]-min[i])/2;
			_cov[i][i] = 1./12*Math.pow(max[i]-min[i],2);
		}
		for (int i=0; i<_dim; i++) {
			for (int j=i+1; j<_dim; j++)  {
				_cov[i][j] = correlation[i][j]*Math.sqrt(_cov[i][i]*_cov[j][j]);
				_cov[j][i] = _cov[i][j];
			}
			_var[i] = _cov[i][i];
		}
	}
	
	static double[][] adjustCorrelation(double[][] correlation) {
		int dim = correlation.length;
		double[][] correl2 = new double[dim][dim];
		for (int i=0; i<dim; i++)  {
			correl2[i][i] = 1;
			for (int j=i+1; j<dim; j++) {
				correl2[i][j] = 2*Math.sin(Math.PI/6.0*correlation[i][j]);
				correl2[j][i] = correl2[i][j];
			}
		}
		return correl2;
	}
	
	public MultivariateUniform(double[] mean, double[][] covariance, Random generator) {
		super(new double[mean.length],covariance,generator);
		double[][] correl = new double[_dim][_dim];
		_min = new double[mean.length];
		_max = new double[mean.length];
		_mean = new double[_dim];
		for (int i=0; i<mean.length; i++) {
			_min[i] = mean[i] - Math.sqrt(3*covariance[i][i]);
			_max[i] = mean[i] + Math.sqrt(3*covariance[i][i]);
			_mean[i] = _min[i] + (_max[i]-_min[i])/2;
			for (int j=0; j<_dim; j++) 
				correl[i][j] = _cov[i][j]/Math.sqrt(_cov[i][i]*_cov[j][j]);
		}
		_chol = choleskyDecomposition(adjustCorrelation(correl));
	}
	
	@Override
	public double[] transform(double[] z) {
		double[] x = new double[_dim];
		for (int i=0; i<_dim; i++) {
			x[i] = 0.0;
			for (int j=0; j<i+1; j++)
				x[i] += z[j]*_chol[i][j];
			x[i] = NormalDistQuick.cdf01(x[i])*(_max[i]-_min[i])+_min[i];
		}
		return x;
	}
	
	@Override
	public double[] getMean() {
		return _mean;
	}

	@Override
	public double[][] getCov() {
		return _cov;
	}

	@Override
	public double[] getSkewness() {
		return new double[_dim];
	}

	@Override
	public double[] getKurtosis() {
		double[] kurt = new double[_dim];
		Arrays.fill(kurt, -6./5);
		return kurt;
	}

	@Override
	public DISTRIBUTION getType() {
		return DISTRIBUTION.Uniform;
	}

}
