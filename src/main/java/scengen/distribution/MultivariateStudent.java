package main.java.scengen.distribution;

import java.util.Random;

import umontreal.ssj.probdist.StudentDistQuick;



public class MultivariateStudent extends MultivariateNormal {
	
	double[][] _cov;
	double[][] _correlation;
	double[] _kurt;
	double[] _scale;
	int[] _df;
	
	public MultivariateStudent(double[] mean, double[] scale, double[][] correlation, int[] df, Random generator) {
		super(mean,correlation,generator);
		_df = df;
		_cov = new double[_dim][_dim];
		_kurt = new double[_dim];
		_scale = scale;
		_correlation = correlation;
		_var = new double[_dim];
		for (int i=0; i<_dim; i++) {
			for (int j=0; j<_dim; j++)
				_cov[i][j] = scale[i]*correlation[i][j]*Math.sqrt(df[i]/(df[i]-2)*df[j]/(df[j]-2));
			_kurt[i] = 6/(df[i]-4);
			_var[i] = _cov[i][i];
		}
	}
	
//	public static void main (String... args) {
//		int dim = 2;
//		int numScen = 100;
//		double[] mean = new double[dim];
//		double[] scale = new double[dim];
//		double[] df = new double[dim];
//		double[][] cov = new double[dim][dim];
//		for (int i=0; i<dim; i++) {
//			df[i] = 5;
//			mean[i] = 0;
//			scale[i] = 1;
//			cov[i][i] = 1;
////			for (int j=0; j<dim; j++)
////				if (i!=j) cov[i][j] = -0.5;
//		}
//		MultivariateDistribution mvdist = new MultivariateStudent(mean,scale,cov,df,new Random());
//		Matrix.print(mvdist.getCov());
//		MultivariateStatistics stat = new MultivariateStatistics(dim);
//		for (int i=0;i<1000000; i++)
//			stat.add(mvdist.getRealization());
//		Matrix.print(stat.getCovariance());
//		System.out.println(Arrays.toString(mvdist.getKurtosis()));
//	}

	
	@Override
	public double[] getRealization() {
		double[] z = new double[_dim];
		for (int i=0; i<_dim; i++) {
			double u,v,w;
			do {
				u = 2*_generator.nextDouble()-1;
				v = 2*_generator.nextDouble()-1;
				w = u*u+v*v;
			} while (w>1);
			double c = u*u/w;
			double r = _df[i]*(Math.pow(w,-2./_df[i])-1);
			z[i] = Math.sqrt(r*c)*(_generator.nextBoolean()?1:-1);
		}
		return transform(z);
	}
	
	@Override
	public double[] getRealization(double[] u) {
		double[] z = new double[_dim];
		for (int i=0; i<_dim; i++)
			z[i] = StudentDistQuick.inverseF(_df[i], u[i]);
		return transform(z);
	}
	
	@Override
	public double[] transform(double[] z) {
		double[] x = new double[_chol.length];
		for (int i=0; i<_dim; i++) {
			double sum = 0.0;
			for (int j=0; j<i+1; j++)
				sum += z[j]*_chol[i][j];
			x[i] = _mean[i] + sum*_scale[i];
		}
		return x;
	}
	
	public double[] getScale() {
		return _scale;
	}
	
	public int[] getDf() {
		return _df;
	}
	
	public double[][] getCorrelation() {
		return _correlation;
	}

	@Override
	public double[][] getCov() {
		return _cov;
	}

	@Override
	public double[] getKurtosis() {
		return _kurt;
	}

	@Override
	public DISTRIBUTION getType() {
		return DISTRIBUTION.Student;
	}

}
