package scengen.measures;

import java.util.Iterator;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import com.google.ortools.linearsolver.*;
import com.google.ortools.linearsolver.MPSolver.ResultStatus;

public class RiskAverseNewsvendor implements DecisionProblem {
	
	
	static { 
	    System.loadLibrary("jniortools");
	}	
	
	Map<double[], Double> _scenarios;
	double _price;
	double _cost;
	double _alpha;
	double[] _quantity;
	double _profit;
	double _var;
	int _dim;
	
	private RiskAverseNewsvendor(Map<double[],Double> referenceScenarios, int dim, double price, double cost, double alpha) {
		_scenarios = referenceScenarios;
		_price = price;
		_alpha = alpha;
		_cost = cost;
		_dim = dim;
	}
	
	public static RiskAverseNewsvendor solve(Map<double[],Double> referenceScenarios, int dim, double price, double cost, double alpha) {
		RiskAverseNewsvendor nwr = new RiskAverseNewsvendor(referenceScenarios,dim,price,cost,alpha);
		nwr.computeCVaRBenders();
		return nwr;
	}
	
	
	public double getCVaR(double[] quantity) {
		SortedMap<Double,Double> profits = new TreeMap<>();
		for (double[] demand : _scenarios.keySet()) {
			double prob = _scenarios.get(demand);
			Double profit = 0.;
			for (int i=0; i<_dim; i++) 
				profit += Math.min(quantity[i],demand[i])*_price - quantity[i]*_cost;
			profit += Math.random()*1.e-6;
			profits.put(profit,prob);
		}
		double cvar = 0.;
		double cum = 0;
		Iterator<Double> iter = profits.keySet().iterator();
		do  {
			double profit = iter.next();
			double prob = profits.get(profit);
			cum += prob;
			cvar += profit*prob;
		} while(cum < _alpha && iter.hasNext());
		return cvar/cum;
	}
	

	
	/*
	 * cut
	 * 0: objective value
	 * 1: shadow prices of var
	 * 2..n+1: shadow prices of buy
	 */
	void getCut(int dim, double price, double cost, double alpha, double[] demand, double[] buy, double var, double prob, double[] cut) {
		double profit = 0.;
		for (int i=0; i<dim; i++) 
			profit += Math.min(buy[i],demand[i])*price - buy[i]*cost;
		if (var>profit) { //constraint is binding so that shadow prices are non-zero
			cut[0] -= prob/alpha*(var-profit);
			cut[1] -= prob/alpha; //increase var by one means deacreasing objective by -1/alpha
			for (int i=0; i<dim; i++) 
				cut[i+2] += prob/alpha*(demand[i]>buy[i] ? price-cost : -cost);
		} //otherwise constraint is non-binding and has zero dual value
	}
	
	
	void computeCVaRBenders() {
		//initialize LP solver
		double maxVar = 10000;
		MPSolver lp = new MPSolver("NewsvendorAtRisk",MPSolver.OptimizationProblemType.valueOf("CLP_LINEAR_PROGRAMMING"));
		double inf = MPSolver.infinity();
		MPVariable varVar = lp.makeNumVar(-inf, inf, "var");
		MPVariable[] buyVar = lp.makeNumVarArray(_dim, -maxVar, maxVar);
		MPVariable costToGo = lp.makeNumVar(-inf, inf, "cost_to_go");
		//add objective function
		MPObjective obj = lp.objective();
		obj.setMaximization();
		obj.setCoefficient(varVar, 1.);
		obj.setCoefficient(costToGo, 0);
		lp.makeConstraint(-maxVar,maxVar).setCoefficient(varVar, 1);
		//get a first solution
		ResultStatus status = lp.solve();
		if (status != ResultStatus.OPTIMAL)
			throw new IllegalStateException("No optimal solution found. Return code: "+status);
		double varSol = varVar.solutionValue();
		double[] buySol = new double[_dim];
		for (int i=0; i<_dim; i++)
			buySol[i] = buyVar[i].solutionValue();
		double[] cut = new double[_dim+2];
		//int dim, double price, double cost, double alpha, double[] demand, double[] buy, double var, double prob, double[] cut
		for (double[] demand :_scenarios.keySet()) 
			getCut(_dim, _price, _cost , _alpha, demand, buySol, varSol, _scenarios.get(demand), cut);
		double ub = obj.value();
		double lb = varSol+cut[0];
		obj.setCoefficient(costToGo, 1);
		//now successively add cuts until convergence is reached
		while(Math.abs(ub-lb) > 0.001) {
			double intercept = cut[0]-cut[1]*varSol;
			for (int i=0; i<_dim; i++)
				intercept -= cut[i+2]*buySol[i];
			MPConstraint cutCtr = lp.makeConstraint(-inf,intercept);
			cutCtr.setCoefficient(costToGo, 1);
			cutCtr.setCoefficient(varVar, -cut[1]);
			for (int i=0; i<_dim; i++)
				cutCtr.setCoefficient(buyVar[i], -cut[i+2]);
			//cut added now solve again
			status = lp.solve();
			if (status != ResultStatus.OPTIMAL)
				throw new IllegalStateException("No optimal solution found. Return code: "+status);
			varSol = varVar.solutionValue();
			buySol = new double[_dim];
			for (int i=0; i<_dim; i++)
				buySol[i] = buyVar[i].solutionValue();
			cut = new double[_dim+2];
			for (double[] demand :_scenarios.keySet()) 
				getCut(_dim, _price, _cost , _alpha, demand, buySol, varSol, _scenarios.get(demand), cut);
			ub = obj.value();
			lb = varSol+cut[0];
		}
		_var = varSol;
		_profit = lb;
		_quantity = buySol;
				
	}

	@Override
	public double getObjective(double[] solution) {
		return getCVaR(solution);
	}

	@Override
	public boolean isFeasible(double[] solution) {
		return true;
	}

	@Override
	public double[] getSolution() {
		return _quantity;
	}

	@Override
	public double getOptimum() {
		return _profit;
	}
	
	

}
