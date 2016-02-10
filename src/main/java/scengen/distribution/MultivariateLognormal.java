package main.java.scengen.distribution;

import java.util.Random;

public class MultivariateLognormal extends MultivariateNormal {
	
	double[] _logmean;
	double[][] _logcov;
	double[] _logvar;
	double[] _logskew;
	double[] _logkurt;
	
	public MultivariateLognormal(double[] mean, double[][] covariance, Random generator) {
		super(mean,covariance,generator);
		_logmean = new double[_dim];
		_logcov = new double[_dim][_dim];
		_logskew = new double[_dim];
		_logkurt = new double[_dim];
		for (int i=0; i<mean.length; i++) {
			_logmean[i] = Math.exp(mean[i]+covariance[i][i]/2);
			for (int j=0; j<mean.length; j++)
				_logcov[i][j] = (Math.exp(covariance[i][j])-1)*Math.exp(2*mean[i]+covariance[i][j]);
			_logskew[i] = (Math.exp(covariance[i][i])+2)*Math.sqrt(Math.exp(covariance[i][i])-1);
			_logkurt[i] = Math.exp(4*covariance[i][i])+2*Math.exp(3*covariance[i][i])+3*Math.exp(2*covariance[i][i])-6;
		}
		_logvar = new double[_dim];
		for (int i=0; i<_dim; i++)
			_logvar[i] = _logcov[i][i];
	}
	
	@Override
	public double[] transform(double[] z) {
		double[] x = super.transform(z);
		for (int i=0; i<_dim; i++)
			x[i] = Math.exp(x[i]);
		return x;
	}
	
	@Override
	public double[] getMean() {
		return _logmean;
	}

	@Override
	public double[][] getCov() {
		return _logcov;
	}
	
	@Override
	public double[] getVariance() {
		return _logvar;
	}

	@Override
	public double[] getSkewness() {
		return _logskew;
	}

	@Override
	public double[] getKurtosis() {
		return _logkurt;
	}

	@Override
	public DISTRIBUTION getType() {
		return DISTRIBUTION.Lognormal;
	}
	
	public double[] getNormalMean() {
		return super._mean;
	}
	
	public double[][] getNormalCov() {
		return super._cov;
	}

}
