package scengen.measures;

import java.util.Map;

import scengen.distribution.MultivariateDistribution;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

public class CovarianceDiscrepancy {
	
	public static double getDiscrepancy(Map<double[],Double> map, MultivariateDistribution dist) {
		return getLogEuclideanRiemannianMetric(getCovariance(dist.getDim(),map),dist.getCov());
	}
	
	public static double[][] getCovariance(int dim, Map<double[],Double> map) {
		double[][] scenCov = new double[dim][dim];
		double[] avg = new double[dim];
		for (double[] scen : map.keySet())
			for (int i=0; i<dim; i++) 
				avg[i] += map.get(scen)*scen[i];
		for (double[] scen : map.keySet()) 
			for (int i=0; i<dim; i++) 
				for (int j=0; j<dim; j++) 
					scenCov[i][j] += map.get(scen)*(scen[i]-avg[i])*(scen[j]-avg[j]);
		return scenCov;
	}
	
	/* V. Arsigny, P. Fillard, X. Pennec, and N. Ayache.
	 * Log- Euclidean metrics for fast and simple calculus on diffusion tensors. 
	 * Magnetic Resonance in Medicine, 56(2):411–421, 2006.
	 */
	public static double getLogEuclideanRiemannianMetric(double[][] cov1, double[][] cov2) {
		double[][] logcov1 = matrixLog(cov1);
		double[][] logcov2 = matrixLog(cov2);
		//Froebenius norm (simplifies to Euclidean)
		double dist = 0.;
		for (int i=0; i<cov1.length; i++) {
			for (int j=0; j<cov1.length; j++) {
				double diff = logcov1[i][j] - logcov2[i][j];
				dist += diff*diff;
			}
		}
		return Math.sqrt(dist);
	}
	
	/*
	 * Returns the log matrix of the given matrix
	 */
	static double[][] matrixLog (double[][] matrix) {
		if (matrix.length==1)
			return new double[][]{{Math.log(matrix[0][0])}};
		if (matrix.length!=matrix[0].length)
			throw new IllegalArgumentException("No sqaure matrix.");
		RealMatrix A = new Array2DRowRealMatrix(matrix);
		EigenDecomposition ev = new EigenDecomposition(A);
		RealMatrix V = ev.getV();
		RealMatrix Vinv = (new LUDecomposition(V).getSolver().getInverse());
		RealMatrix logA = Vinv.multiply(A).multiply(V);
		for (int i=0; i<matrix.length; i++)
			logA.setEntry(i, i,Math.log(logA.getEntry(i, i)));
		return V.multiply(logA).multiply(Vinv).getData();
	}

}
