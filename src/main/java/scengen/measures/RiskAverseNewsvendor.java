package main.java.scengen.measures;

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
	
//	public double  getCVaR(double[] quantity) {
//		double cvar = 0;
//		for (double[] demand : _scenarios.keySet()) {
//			double prob = _scenarios.get(demand);
//			double profit = 0.;
//			for (int i=0; i<_dim; i++) 
//				profit += Math.min(quantity[i],demand[i])*_price - quantity[i]*_cost;
//			cvar += prob*Math.max(0,_var-profit);
//		}
//		return _var - cvar/_alpha;
//	}
	
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
	
//	public double[] getError(Map<double[],Double> scens) {
//		NewsvendorAtRisk scenNewsvendor = NewsvendorAtRisk.solve(scens,_dim,_price,_cost,_alpha);
//		double expProfit = scenNewsvendor.getOptimum();
//		double actProfit = getCVaR(scenNewsvendor.getOptimalQuantity());
//		return new double[]{expProfit,actProfit};
//	}
//	
//	public double getError2(Map<double[],Double> scens) {
//		NewsvendorAtRisk scenNewsvendor = NewsvendorAtRisk.solve(scens,_dim,_price,_cost,_alpha);
//		return getCVaR(scenNewsvendor.getOptimalQuantity());
//	}

//	public static void main(String[] args) {
//		int dim = 25;
//		int numScen = 100000;
//		double[] mean = new double[dim];
//		double[][] cov = new double[dim][dim];
//		for (int i=0; i<dim; i++) {
//			mean[i] = 20;
//			cov[i][i] = 25;
//		}
//		MultivariateDistribution mvdist = new MultivariateNormal(mean,cov,new Xorshift(7));
//		Map<double[],Double> scens = new MonteCarlo(mvdist,mvdist.getGenerator()).getScenarios(numScen);
//		
//
//		NewsvendorAtRisk nwr = NewsvendorAtRisk.solve(scens,dim,1,0.95,1.0);
//		Newsvendor nw = Newsvendor.solveNormalNewsvendor(1, 0.95, 20, 5);
////		NewsvendorAtRisk ref = NewsvendorAtRisk.solve(new MonteCarlo(mvdist,mvdist.getGenerator()).getScenarios(10000),mvdist,10000,191,0.05))
//		System.out.println(nw.getOptimalProfit()+" "+nwr.getOptimum()/dim);
//		System.out.println(nwr.getCVaR(nwr.getOptimalQuantity())/dim);
////		System.out.println(Arrays.toString(NewsvendorAtRisk.getError(scens,mvdist,10000,191,0.05)));
//		
//	}
	
//	void computeCVaRCplex() {
//		int numScen = _scenarios.size();;
//		double[] prob = new double[numScen];
//		double[][] demand = new double[numScen][];
//		int k=0;
//		for (double[] scen : _scenarios.keySet()) {
//			demand[k] = scen;
//			prob[k] = _scenarios.get(scen);
//			k++;
//		}
//		int dim = demand[0].length;
//		
//		try {
//		
//			IloCplex lp = new IloCplex();
//			lp.setParam(IloCplex.IntParam.RootAlgorithm,IloCplex.Algorithm.Barrier);
//			lp.setParam(IloCplex.IntParam.BarCrossAlg,-1);
//			lp.setParam(IloCplex.IntParam.SimDisplay, 0);
//			lp.setParam(IloCplex.IntParam.BarDisplay, 0);
//	//		lp.setOut(null); 
//			//create decision variables
//			LinkedList<IloNumVar> vars = new LinkedList<>();
//			IloObjective obj = lp.addMaximize();
//			IloLPMatrix  mat = lp.addLPMatrix();
//			IloNumVar varVar = lp.numVar(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, "var");
//			lp.setLinearCoef(obj,1,varVar);
//			vars.add(varVar);
//			int j=0;
//			int var = j++;
//			int[] buy = new int[dim];
//			IloNumVar[] buyVar = new IloNumVar[_dim];
//			for (int d=0; d<dim; d++) {
//				buyVar[d] = lp.numVar(0, Double.POSITIVE_INFINITY, "buy_"+d);
//				vars.add(buyVar[d]);
//				buy[d] = j++;
//			}
//			int[] tailLoss = new int[numScen];
//			for (int n=0; n<numScen; n++) {
//				IloNumVar y = lp.numVar(0, Double.POSITIVE_INFINITY,"tailLoss_"+n);
//				lp.setLinearCoef(obj,-prob[n]/_alpha,y);
//				tailLoss[n] = j++;
//				vars.add(y);
//			}
//			int[][] sell = new int[numScen][dim];
//			for (int n=0; n<numScen; n++) 
//				for (int d=0; d<dim; d++) {
//					sell[n][d] = j++;
//					vars.add(lp.numVar(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, "sell_"+n+"_"+d));
//				}
//			mat.addCols(vars.toArray(new IloNumVar[0]));
//			
//			double[] lbs = new double[numScen*2*dim+numScen];
//			double[] ubs = new double[numScen*2*dim+numScen];
//			int[][] x = new int[numScen*2*dim+numScen][];
//			double[][] a = new double[numScen*2*dim+numScen][];
//			int i=0;
//			for (int s=0; s<numScen; s++) {
//				x[i] = new int[2+2*dim];
//				a[i] = new double[2+2*dim];
//				int m=0;
//				a[i][m]=1; x[i][m] = tailLoss[s]; m++;
//				a[i][m]=-1; x[i][m] = var; m++;
//				for (int d=0; d<dim; d++) {
//					a[i][m] = _price; x[i][m] = sell[s][d]; m++;
//				}
//				for (int d=0; d<dim; d++) {
//					 a[i][m] = -_cost; x[i][m] = buy[d]; m++;
//				}
//				ubs[i] = Double.POSITIVE_INFINITY;
//				i++;
//				for (int d=0; d<dim; d++) {
//					x[i] = new int[]{sell[s][d], buy[d]};
//					a[i] = new double[]{1,-1};
//					lbs[i] = Double.NEGATIVE_INFINITY;
//					i++;
//					x[i] = new int[]{sell[s][d]};
//					a[i] = new double[]{1};
//					lbs[i] = Double.NEGATIVE_INFINITY;
//					ubs[i] = demand[s][d];
//					i++;
//				}
//			}
//			mat.addRows(lbs, ubs, x, a);
//			boolean status = lp.solve(); 
//			if (!status) {
////				lp.exportModel("newsvendor.lp");
//				throw new IllegalStateException("No optimal solution found.");
//			}
//			_profit = lp.getObjValue();
//			_quantity = lp.getValues(buyVar);
//			_var = lp.getValue(varVar);
//		}  catch (IloException e) {e.printStackTrace();}
//	}
	
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
	
//	/*
//	 * cut
//	 * 0: objective value
//	 * 1: shadow prices of var
//	 * 2..n+1: shadow prices of buy
//	 */
//	public void getCut2(int dim, double price, double cost, double alpha, double[] demand, double[] buy, double var, double prob, double[] cut) {
//		MPSolver lp = new MPSolver("NewsvendorAtRisk",MPSolver.OptimizationProblemType.valueOf("CLP_LINEAR_PROGRAMMING"));
//		MPObjective obj = lp.objective();
//		obj.setMaximization();
//		MPVariable[] sell = new MPVariable[dim];
//		MPConstraint[] sellCtr = new MPConstraint[dim];
//		for (int i=0; i<dim; i++) {
//			sell[i] = lp.makeNumVar(demand[i]<0?demand[i]:0, demand[i]>0?demand[i]:0, "sell_"+i);
//			sellCtr[i] = lp.makeConstraint(Double.NEGATIVE_INFINITY,buy[i]);
//			sellCtr[i].setCoefficient(sell[i], 1.);
//		}
//		MPVariable tailLoss = lp.makeNumVar(0, Double.POSITIVE_INFINITY, "tailLoss");
//		obj.setCoefficient(tailLoss,-1./alpha);
//		double costs = 0;
//		for (int i=0; i<dim; i++) 
//			costs += buy[i]*cost;
//		MPConstraint tailLossCtr = lp.makeConstraint(var+costs, Double.POSITIVE_INFINITY);
//		tailLossCtr.setCoefficient(tailLoss, 1.);
//		for (int i=0; i<dim; i++) 
//			tailLossCtr.setCoefficient(sell[i], price);
//		ResultStatus status = lp.solve();
//		if (status != ResultStatus.OPTIMAL)
//			throw new IllegalStateException("No optimal solution found. Return code: "+status);
//		cut[0] += prob*obj.value();
//		double tailLossDual = tailLossCtr.dualValue();
//		cut[1] += prob*tailLossDual;
//		for (int i=0; i<dim; i++)
//			cut[i+2] += prob*(sellCtr[i].dualValue()+tailLossDual);
//	}
	
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
//		System.out.println(ub+"\t"+lb+"\t"+lp.numConstraints()+"\t"+Arrays.toString(cut));
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
			//int dim, double price, double cost, double alpha, double[] demand, double[] buy, double var, double prob, double[] cut
			for (double[] demand :_scenarios.keySet()) 
				getCut(_dim, _price, _cost , _alpha, demand, buySol, varSol, _scenarios.get(demand), cut);
			ub = obj.value();
			lb = varSol+cut[0];
//			System.out.println(ub+"\t"+lb+"\t"+lp.numConstraints()+"\t"+Arrays.toString(cut));
		}
//		System.out.println(lb);
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
