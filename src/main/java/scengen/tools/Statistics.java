package main.java.scengen.tools;

import java.io.Serializable;

import org.apache.commons.math3.distribution.TDistribution;

/**
 * The {@code Statistics} object keeps track of the first two moments of a stream of data points.
 * 
 * @author Nils Loehndorf
 */
public class Statistics implements Serializable {
	
	private static final long serialVersionUID = 677360888189461348L;
	int count;
	double avgX;
	public double avgX2;
	
	/**
	 * Adds a sample point to the list.
	 * @param X add a data point
	 */
	public synchronized void add(double X){
		if (Double.isNaN(X)) return;
		double a = 1./++count;
		avgX += a*(X-avgX);
		avgX2 += a*(X*X-avgX2);
	}

	
	/**
	 * Returns the mean.
	 * @return batch mean
	 */
	public double getMean() {
		if (count==0) return Double.NaN;
		return avgX;
	}
	
	/**
	 * Returns the mean.
	 * @return batch mean
	 */
	public double getMeanOfSquares() {
		if (count==0) return Double.NaN;
		return avgX2;
	}
	
	/**
	 * Returns the variance.
	 * @return variance
	 */
	public double getVariance() {
		if (count<=1) return Double.NaN;
		return (avgX2-avgX*avgX)*count/(count-1);
	}
	
	/**
	 * Returns the standard deviation.
	 * @return variance
	 */
	public double getStandardDeviation() {
		if (count<=1) return Double.NaN;
		return Math.sqrt(getVariance());
	}
	
	/**
	 * Returns the helf width of the confidence interval (e.g, for the 95% confidence interval alpha=0.95).
	 * @return variance
	 */
	public double getHalfwidth(double alpha) {
		if (count<5) return Double.NaN;
		TDistribution t = new TDistribution(count);
		return t.inverseCumulativeProbability((1.+alpha)/2.)*getStandardError();
	}
	
	/**
	 * Returns the standard error.
	 * @return standard error
	 */
	public double getStandardError() {
		if (count<=1) return Double.NaN;
		return Math.sqrt((avgX2-avgX*avgX)/(count-1));
	}
	
	/**
	 * @return Returns the coefficient of variation (CV).
	 */
	public double getCV() {
		if (count<=1) return Double.NaN;
		return getStandardDeviation()/getMean();
	}
	
	/**
	 * Deletes the sample and resets the variables.
	 */
	public void reset() {
		avgX = 0.0; avgX2 = 0.0; count = 0;
	}
	
	public int size() {
		return count;
	}
	
}
