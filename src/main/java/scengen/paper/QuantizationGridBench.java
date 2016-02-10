package main.java.scengen.paper;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.Map;
import java.util.Random;

import main.java.scengen.application.Generator;
import main.java.scengen.distribution.MultivariateDistribution;
import main.java.scengen.distribution.MultivariateNormal;
import main.java.scengen.measures.WassersteinDistance;
import main.java.scengen.methods.METHOD;
import main.java.scengen.methods.MonteCarlo;
import main.java.scengen.tools.Statistics;
import main.java.scengen.tools.Xorshift;

public class QuantizationGridBench {
	
	static int _sampleSize = 100000;
	static int _numRepititions = 10;
	static long _seed = 191;
	static Random _rand = new Xorshift(_seed);
	static String filename = "comparison_"+System.nanoTime()+".txt";
	static int[][] dimScens = new int[][]{{2,4},{2,40},{10,20},{10,200}};
		
	public static void main (String... args) throws IOException {	
		FileOutputStream stream = new FileOutputStream(filename);
		BufferedWriter br = new BufferedWriter(new OutputStreamWriter(stream));
		String s = "Dimensions\tScenarios\tMethod\tWasserstein\tSE\n";
		br.write(s);br.flush();
		System.out.print(s);
		for (int[] dimScen : dimScens) {
			int dim = dimScen[0];
			int numScen = dimScen[1];
			double[] means = new double[dim];
			double[][] cov = new double[dim][dim];
			for (int i=0; i<dim; i++)
				cov[i][i] = 1;
			MultivariateDistribution mvDist =  new MultivariateNormal(means,cov,_rand);	
			for (METHOD method : new METHOD[]{METHOD.QuantizationGrids,METHOD.QuantizationLearning}) {
				Statistics wasserStat = new Statistics();
				for (int i=0; i<_numRepititions; i++) {
					Map<double[],Double> referenceScen = new MonteCarlo(mvDist).getScenarios(_sampleSize);
					Map<double[],Double> scen = Generator.runScenarioGenerator(numScen,method,mvDist,_rand.nextLong());
					wasserStat.add(WassersteinDistance.getLowerBound(mvDist.getDim(),scen,referenceScen.keySet()));
				}
				s = dim+"\t"+numScen+"\t"+method.toString()+"\t";
				s += wasserStat.getMean()+"\t"+wasserStat.getStandardError()+"\n";
				br.write(s);br.flush();
				System.out.print(s);
			}						
		}
		br.close();
	}

}
