package main.java.scengen.paper;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
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
import main.java.scengen.methods.METHOD;
import main.java.scengen.tools.Statistics;
import main.java.scengen.tools.Xorshift;

public class MomentMatchingNaN {
	
	static int _baseDemand = 1;
	static int _sampleSize = 1000;
	static int _numRepititions = 1;
	static long _seed = 191;
	static Random _rand = new Xorshift(_seed);
	static int _df = 5; //only Student distribution
	static String filename = "experiment_"+System.nanoTime()+".txt";
	static double _alpha = 0.05;
	static double[] correlation = new double[]{0.0,0.5};
	static double[] coefficientOfVariation = new double[]{0.3,0.7};
	static int[][] dimScens = new int[][]{{2,5},{2,50},{10,25},{10,250},{25,50},{25,500}};
	static DISTRIBUTION[] dists = new DISTRIBUTION[]{DISTRIBUTION.Lognormal,DISTRIBUTION.Student};
	static METHOD[] methods = new METHOD[]{METHOD.MomentMatching};
		
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
//						Map<double[],Double> referenceScen = new MonteCarlo(mvDist,mvDist.getGenerator()).getScenarios(_sampleSize);
						Map<METHOD,Statistics> covdet = new HashMap<>();
						//generate scenarios using the given methods
						for (METHOD method : methods) {
							List<Map<double[],Double>> scens = new ArrayList<>(_numRepititions);
							Statistics covdetStat = new Statistics();
							for (int i=0; i<_numRepititions; i++) {
								try {
									//use Sylvain's grids if applicable
									METHOD m = method;
									Map<double[],Double> scen = Generator.runScenarioGenerator(numScen,m,mvDist,_rand.nextLong());
									if (scen!=null) scens.add(scen);
									covdetStat.add(CovarianceDiscrepancy.getDiscrepancy(scen, mvDist));
								} catch (Exception e) {
									e.printStackTrace();
								}
							}
							covdet.put(method,covdetStat);
						}
							//solve the decision problems
							//evaluate each method on the given problem	
						for (METHOD method : methods) {
							s = dist.toString()+"\t"+dim+"\t"+numScen+"\t"+cv+"\t"+cor+"\t"+method.toString()+"\t";
							Statistics covdetStat = covdet.get(method);
							s += covdetStat.getMean();
							s += "\t"+covdetStat.getStandardError();
							s += "\t"+covdetStat.size()+"\n";
							br.write(s);br.flush();
							System.out.print(s);
						}
					}
				}
			}
		}
		br.close();
	}

}
