package scengen.distribution;

import java.util.List;
import java.util.Random;

public interface MultivariateDistribution {
	
	public double[] getRealization();
	
	public double[] getRealization(double[] u);
	
	public int getDim();
	
	public List<double[]> getSample(int size);
	
	public void setGenerator(Random generator);
	
	public Random getGenerator();
	
	public double[] getMean();
	
	public double[] getVariance();
	
	public double[][] getCov();
	
	public double[] getSkewness();
	
	public double[] getKurtosis();
	
	public DISTRIBUTION getType();

}
