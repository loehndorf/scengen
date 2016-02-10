package main.java.scengen.paper;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import main.java.scengen.application.Generator;
import main.java.scengen.distribution.DISTRIBUTION;
import main.java.scengen.distribution.MultivariateDistribution;
import main.java.scengen.distribution.MultivariateLognormal;
import main.java.scengen.distribution.MultivariateNormal;
import main.java.scengen.distribution.MultivariateStudent;
import main.java.scengen.distribution.MultivariateUniform;
import main.java.scengen.measures.CovarianceDiscrepancy;
import main.java.scengen.measures.DecisionProblem;
import main.java.scengen.measures.MultiNewsvendor;
import main.java.scengen.measures.PROBLEM;
import main.java.scengen.measures.RiskAverseNewsvendor;
import main.java.scengen.measures.WassersteinDistance;
import main.java.scengen.methods.METHOD;
import main.java.scengen.methods.MonteCarlo;
import main.java.scengen.tools.Statistics;
import main.java.scengen.tools.Xorshift;

public class Experiment {
	
	static int _baseDemand = 1;
	static int _sampleSize = 100000;
	static boolean _useGrids = true;
	static int _numRepititions = 5;
	static long _seed = 191;
	static Random _rand = new Xorshift(_seed);
	static int _df = 5; //only Student distribution
	static String filename = "experiment_"+System.nanoTime()+".txt";
	static double _alpha = 0.05;
	static double[] correlation = new double[]{0.0,0.5};
	static double[] coefficientOfVariation = new double[]{0.3,0.7};
	static double[] fractiles = new double[]{0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
	static int[][] dimScens = new int[][]{{2,5},{2,50},{10,25},{10,250},{20,50},{20,500}};
	static DISTRIBUTION[] dists = new DISTRIBUTION[]{DISTRIBUTION.Normal,DISTRIBUTION.Uniform,DISTRIBUTION.Lognormal,DISTRIBUTION.Student};
	static METHOD[] methods = new METHOD[]{METHOD.MomentMatching,METHOD.QuantizationLearning,METHOD.QuasiMonteCarlo,METHOD.VoronoiCellSampling};
	static PROBLEM[] problems = new PROBLEM[]{PROBLEM.MultiNewsvendor,PROBLEM.RiskAverseNewsvendor};
		
	public static void main (String... args) throws IOException {	
		FileOutputStream stream = new FileOutputStream(filename);
		BufferedWriter br = new BufferedWriter(new OutputStreamWriter(stream));
		String s = "Distribution\tDimensions\tScenarios\tCV\trho\tfractile\tMeasure\tMethod\tOptimum\tExpected\tActual\tWasserstein\tCovDet\tExpSE\tActSE\tWasserSE\tCovDetSE\tRepititions\n";
		br.write(s);br.flush();
		System.out.print(s);
		for (DISTRIBUTION dist : dists) {
			for (int[] dimScen : dimScens) {
				int dim = dimScen[0];
				int numScen = dimScen[1];
				for (double cv : coefficientOfVariation) {
					for (double cor : correlation) {
						MultivariateDistribution mvDist = null;
						double[] means = new double[dim];
						double[][] cov = new double[dim][dim];
						for (int i=0; i<dim; i++) {
							means[i] = _baseDemand;
							cov[i][i] = _baseDemand*cv*_baseDemand*cv;
						}
						for (int i=0; i<dim; i++) {
							for (int j=i+1; j<dim; j++) {
								cov[i][j] = Math.sqrt(cov[i][i]*cov[j][j])*cor;
								cov[j][i] = cov[i][j];
							}
						}
						switch(dist) {
						case Normal : 
							mvDist = new MultivariateNormal(means,cov,_rand);	
							break;
						case Lognormal : 
							mvDist = new MultivariateLognormal(means,cov,_rand);	
							break;
						case Uniform : 
							mvDist = new MultivariateUniform(means,cov,_rand);	
							break;
						case Student : 
							int[] df = new int[dim];
							double[] scale = new double[dim];
							double[][] correl = new double[dim][dim];
							for (int i=0; i<dim; i++) {
								scale[i] = Math.sqrt(cov[i][i]);
								for (int j=0; j<dim; j++)
									if (i!=j)
										correl[i][j] = cor;
									else
										correl[i][j] = 1;
							}
							
							Arrays.fill(df,_df);
							mvDist = new MultivariateStudent(means,scale,correl,df,_rand);	
							break;
						}
						Map<double[],Double> referenceScen = new MonteCarlo(mvDist).getScenarios(_sampleSize);
						Map<METHOD,Collection<Map<double[],Double>>> scenMap = new HashMap<>();
						Map<METHOD,Statistics> wasserstein = new HashMap<>();
						Map<METHOD,Statistics> covdet = new HashMap<>();
						//generate scenarios using the given methods
						for (METHOD method : methods) {
							List<Map<double[],Double>> scens = new ArrayList<>(_numRepititions);
							Statistics wasserStat = new Statistics();
							Statistics covdetStat = new Statistics();
							for (int i=0; i<_numRepititions; i++) {
								try {
									//use Sylvain's grids if applicable
									METHOD m = method;
									if (method==METHOD.QuantizationLearning && _useGrids && mvDist.getDim()<=10 && numScen<=1000 && mvDist.getType()==DISTRIBUTION.Normal && ((MultivariateNormal)mvDist).isIndependent())
										m = METHOD.QuantizationGrids;
									Map<double[],Double> scen = Generator.runScenarioGenerator(numScen,m,mvDist,_rand.nextLong());
									if (scen!=null) scens.add(scen);
									covdetStat.add(CovarianceDiscrepancy.getDiscrepancy(scen, mvDist));
									wasserStat.add(WassersteinDistance.getLowerBound(mvDist.getDim(),scen,referenceScen.keySet()));
								} catch (Exception e) {
									e.printStackTrace();
								}
							}
							wasserstein.put(method, wasserStat);
							covdet.put(method,covdetStat);
							scenMap.put(method, scens);
						}
						for (double fractile : fractiles) {
							double price = 1;
							double cost = 1-fractile;
							//solve the decision problems
							Map<PROBLEM,DecisionProblem> problemMap = new HashMap<>();
							for (PROBLEM problem : problems) {
								if (problem==PROBLEM.MultiNewsvendor) {
									MultiNewsvendor mnw = null;
									switch(dist) {
										case Normal :  mnw = MultiNewsvendor.solveNormal((MultivariateNormal)mvDist, price, cost); break;
										case Lognormal :  mnw = MultiNewsvendor.solveLognormal((MultivariateLognormal)mvDist, price, cost); break;
										case Uniform :  mnw = MultiNewsvendor.solveUniform((MultivariateUniform)mvDist, price, cost); break;
										case Student :  mnw = MultiNewsvendor.solveStudent((MultivariateStudent)mvDist, price, cost, referenceScen.keySet()); break;
									}
									problemMap.put(problem,mnw);
								}
								if (problem==PROBLEM.RiskAverseNewsvendor) {
									problemMap.put(problem,RiskAverseNewsvendor.solve(referenceScen, dim, price, cost, _alpha));
								} 
							}
							//evaluate each method on the given problem	
							for (PROBLEM problem : problems) {
								for (METHOD method : methods) {
									s = dist.toString()+"\t"+dim+"\t"+numScen+"\t"+cv+"\t"+cor+"\t"+fractile+"\t"
											+problem.toString()+"\t"+method.toString()+"\t";
									DecisionProblem probTrue = problemMap.get(problem);
									double trueOpt = probTrue.getOptimum();
									Statistics approxOpt = new Statistics();
									Statistics actualOpt = new Statistics();
									Statistics wasserStat= wasserstein.get(method);
									Statistics covdetStat = covdet.get(method);
									for (Map<double[],Double> scen : scenMap.get(method)) {
										DecisionProblem probApprox = null;
										if (problem==PROBLEM.MultiNewsvendor) 
											probApprox = MultiNewsvendor.solveEmpirical(dim,scen, price, cost);
										if (problem==PROBLEM.RiskAverseNewsvendor)
											probApprox = RiskAverseNewsvendor.solve(scen, dim, price, cost, _alpha);
										approxOpt.add(probApprox.getOptimum());
										actualOpt.add(probTrue.getObjective(probApprox.getSolution()));
									}
									s += trueOpt+"\t"+approxOpt.getMean()+"\t"+actualOpt.getMean()+"\t"+wasserStat.getMean()+"\t"+covdetStat.getMean();
									s += "\t"+approxOpt.getStandardError()+"\t"+actualOpt.getStandardError()+"\t"+wasserStat.getStandardError()+"\t"+covdetStat.getStandardError();
									s += "\t"+covdetStat.size()+"\n";
									br.write(s);br.flush();
									System.out.print(s);
								}
							}
						}
					}
				}
			}
		}
		br.close();
	}

}
