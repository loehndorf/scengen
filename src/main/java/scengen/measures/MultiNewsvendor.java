package scengen.measures;

import java.util.Collection;
import java.util.Map;

import scengen.distribution.MultivariateLognormal;
import scengen.distribution.MultivariateNormal;
import scengen.distribution.MultivariateStudent;
import scengen.distribution.MultivariateUniform;
import scengen.methods.ReductionMethod;

public class MultiNewsvendor implements DecisionProblem {
	
	Newsvendor[] _problems;
	int _dim;
	
	private MultiNewsvendor(Newsvendor[] problems) {
		_dim = problems.length;
		_problems = problems;
	}
	
	public static MultiNewsvendor solveNormal(MultivariateNormal mvDist, double price, double cost) {
		double[] mean = mvDist.getMean();
		double[][] cov = mvDist.getCov();
		int dim = mvDist.getDim();
		Newsvendor[] problems = new Newsvendor[dim];
		for (int i=0; i<dim; i++) 
			problems[i] = Newsvendor.solveNormalNewsvendor(price, cost, mean[i], Math.sqrt(cov[i][i]));
		return new MultiNewsvendor(problems);
	}
	
	public static MultiNewsvendor solveLognormal(MultivariateLognormal mvDist, double price, double cost) {
		double[] mean = mvDist.getMean();
		double[][] cov = mvDist.getCov();
		int dim = mvDist.getDim();
		Newsvendor[] problems = new Newsvendor[dim];
		for (int i=0; i<dim; i++) 
			problems[i] = Newsvendor.solveLognormalNewsvendor(price, cost, mean[i], Math.sqrt(cov[i][i]));
		return new MultiNewsvendor(problems);
	}
	
	public static MultiNewsvendor solveUniform(MultivariateUniform mvDist, double price, double cost) {
		double[] mean = mvDist.getMean();
		double[][] cov = mvDist.getCov();
		int dim = mvDist.getDim();
		Newsvendor[] problems = new Newsvendor[dim];
		for (int i=0; i<dim; i++) 
			problems[i] = Newsvendor.solveUniformNewsvendor(price, cost, mean[i], Math.sqrt(cov[i][i]));
		return new MultiNewsvendor(problems);
	}
	
	public static MultiNewsvendor solveStudent(MultivariateStudent mvDist, double price, double cost, Collection<double[]> sample) {
		double[] mean = mvDist.getMean();
		int dim = mvDist.getDim();
		Newsvendor[] problems = new Newsvendor[dim];
		for (int i=0; i<dim; i++) {
			double[] col = new double[sample.size()];
			int n=0;
			for (double[] d : sample)
				col[n++] = d[i];
			problems[i] = Newsvendor.solveStudentNewsvendor(price, cost, mean[i], mvDist.getScale()[i],mvDist.getDf()[i],col);
		}
		return new MultiNewsvendor(problems);
	}
	
	public static MultiNewsvendor solveEmpirical(int dim, Map<double[],Double> scens, double price, double cost) {
		Newsvendor[] problems = new Newsvendor[dim];
		for (int i=0; i<dim; i++) {
			double[] scenarios = ReductionMethod.getScenarios(i, scens);
			double[] weights = ReductionMethod.getWeights(scens);
			problems[i] = Newsvendor.solveEmpiricalNewsvendor(price,cost,scenarios,weights);
		}
		return new MultiNewsvendor(problems);
	}

	@Override
	public double[] getSolution() {
		double[] quantity = new double[_dim];
		for (int i=0; i<_dim; i++) 
			quantity[i] = _problems[i].getOptimalQuantity();
		return quantity;
	}

	@Override
	public double getObjective(double[] solution) {
		double actProfit = 0.;
		for (int i=0; i<_dim; i++) 
			actProfit += _problems[i].getProfit(solution[i]);
		return actProfit;
	}

	@Override
	public boolean isFeasible(double[] solution) {
		return true;
	}

	@Override
	public double getOptimum() {
		double sum = 0.;
		for (int i=0; i<_dim; i++) 
			sum += _problems[i].getOptimalProfit();
		return sum;
	}

}
