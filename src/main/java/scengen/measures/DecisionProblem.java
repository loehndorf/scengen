package main.java.scengen.measures;

public interface DecisionProblem {
	
	/**
	 * @return optimal objective value
	 */
	public double getOptimum();
	
	/**
	 * @param solution
	 * @return objective value for the given solution vector
	 */
	public double getObjective(double[] solution);
	
	/**
	 * 
	 * @param solution
	 * @return true if solution is feasible, false otherwise
	 */
	public boolean isFeasible(double[] solution);
	
	/**
	 * @return optimal solution
	 */
	public double[] getSolution();
	
}
