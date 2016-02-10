package main.java.scengen.methods;

import java.util.Map;

public abstract class ReductionMethod {
	
	public abstract Map<double[],Double> getScenarios(int numScen);
	
	/**
	 * @param dim
	 * @param scenarios
	 * @return Returns the scenario along one particular dimension.
	 */
	public static double[] getScenarios(int dim, Map<double[],Double> scenarios) {
		double[] scen = new double[scenarios.size()];
		int i=0;
		for (double[] x : scenarios.keySet())
			scen[i++] = x[dim];
		return scen;
	}
	
	public static double[] getWeights(Map<double[],Double> scenarios) {
		double[] weights = new double[scenarios.size()];
		int i=0;
		for (double[] x : scenarios.keySet())
			weights[i++] = scenarios.get(x);
		return weights;
	}

}
