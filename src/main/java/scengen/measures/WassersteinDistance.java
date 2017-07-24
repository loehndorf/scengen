package scengen.measures;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;

import scengen.tools.Metrics;

public class WassersteinDistance  {
	
	static String fileToString(String filename) {
		StringBuilder sb = new StringBuilder();
		try {
			FileInputStream in = new FileInputStream(filename);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			String s;
			while ((s = br.readLine())!=null) {
				sb.append(s).append("\n");
			}
			br.close();
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return sb.toString();
	}
	
//	public static double getExactDistance(Map<double[],Double> scenarios, Collection<double[]> sample) {
//		Env env = null;
//		Model model = null;
////		double time = 0;
//		try {
//			int sampleSize = sample.size();
//			env = new Env();
////			env.PutLicenseString("<"+fileToString("lib/sulum/sulumlp.lic")+">");
//			model = new Model(env);
//			model.SetIntParam(SlmParamInt.SlmPrmIntPresolve,SlmOnOff.SlmOff.get());
////			model.SetIntParam(SlmParamInt.SlmPrmIntSimMaxIter,1000);
////			model.SetIntParam(SlmParamInt.SlmPrmIntUpdateSolQuality,SlmOnOff.SlmOn.get());
//			model.SetIntParam(SlmParamInt.SlmPrmIntSimUseNetwork,SlmOnOff.SlmOn.get());
//			model.SetIntParam(SlmParamInt.SlmPrmIntLogLevel,0);
////			model.SetIntParam(SlmParamInt.SlmPrmIntDebug,SlmOnOff.SlmOn.get());
//			int numScen = scenarios.size();
//			int cons = sampleSize+numScen;
//			int vars = sampleSize*numScen;
//			//constraints
//			SlmBoundKey[] cbnd = new SlmBoundKey[cons];
//			int[] abegin = new int[cons+1];
//			int[] aidxj = new int[vars*2];
//			double[] avalj = new double[vars*2];
//			Arrays.fill(avalj, 1.0);
//			double[] clo = new double[cons];
//			double[] cup = clo;
//			for (int i=0; i<cons; i++) {
//				cbnd[i] = SlmBoundKey.SlmBndFx;
//				if (i<numScen){
//					abegin[i] = (i>0?abegin[i-1]+sampleSize:0);
//				}
//				else {
//					clo[i] = -1./sampleSize;
//					abegin[i] = abegin[i-1]+(i==numScen?sampleSize:numScen);
//				}
//			}
//			abegin[cons] = abegin[cons-1]+numScen;
//			//variables
//			double[] obj = new double[vars];
//			double[] vlo = new double[vars];
//			double[] vup = new double[vars];
//			SlmBoundKey[] vbnd = new SlmBoundKey[vars];
//			//construct matrix
//			for (int i=0; i<numScen; i++) 
//				for (int j=0; j<sampleSize; j++) {
//					aidxj[i*sampleSize+j] = i*sampleSize+j;
//					avalj[i*sampleSize+j] = 1;
//				}
//			for (int i=0; i<sampleSize; i++) 
//				for (int j=0; j<numScen; j++) {
//					aidxj[vars+i*numScen+j] = j*sampleSize+i;
//					avalj[vars+i*numScen+j] = -1;
//				}
//			int i=0; double probSum = 0.;
//			//create variables and objective
//			for (double[] d1 : scenarios.keySet()) {
//				double prob = scenarios.get(d1);
//				clo[i] = prob;
//				if (i==numScen-1) clo[i] = 1-probSum;
//				probSum += clo[i];
//				int j=0;
//				for (double[] d2 : sample) {
//					int idx = i*sampleSize+j;
//					//variable bounds
//					vbnd[idx] = SlmBoundKey.SlmBndRa;
//					vup[idx] = 1./sampleSize;
//					//objective
//					obj[idx] = Metrics.squaredDistance(d1,d2);
//					j++;
//				}
//				i++;
//			}
////			
////			model.SetIntParam(SlmParamInt.SlmPrmIntObjSense,SlmObjSense.SlmObjSenseMin.get());
//			model.SetAllData(cons,vars,0.0,cbnd,clo,cup,vbnd,obj,vlo,vup,0,abegin,aidxj,avalj);
////			model.SetIntParam(SlmParamInt.SlmPrmIntObjSense,SlmObjSense.SlmObjSenseMin.get());
////			System.err.println("Model time: "+(System.nanoTime()-time)/1000000000.0);
////			time = System.nanoTime();
////			model.WriteProblem("wasserstein_log_dim25_scen500_cv1.0_rho0.0.mps.gz");
////			System.err.println("Save time: "+(System.nanoTime()-time)/1000000000.0);
////			System.out.println(model.WriteStringProblem(SlmCompressionType.SlmCompressionTypeNone, SlmEncodingType.SlmEncodingTypeNone, SlmProblemFormat.SlmProblemFormatLp));
//			model.Optimize();
////			model.PrintSolQuality(SlmStream.SlmStrInfo);
////			model.WriteSolution("wasserstein_neg5_simplex.sol");
//			return Math.sqrt(model.GetDbInfo(SlmInfoDb.SlmInfoDbPrimObj));
//		} catch (Exception e) {
//			e.printStackTrace();
//			try {
////				System.err.println("Solver time: "+(System.nanoTime()-time)/1000000000.0);
//				return Math.sqrt(model.GetDbInfo(SlmInfoDb.SlmInfoDbPrimObj));
//			} catch (Exception e1) {
//				e1.printStackTrace();
//			}
//			
//		}
//		finally {
//			env.Dispose();
//			model.Dispose();
//		}
//		return Double.NaN;
//	}
	
	
	public static double getLowerBound(int dim, Map<double[],Double> scen, Collection<double[]> sample) {
		int sampleSize = sample.size();
		scengen.tools.MetricTree<Cluster> tree = new scengen.tools.MetricTree<>(dim);
		for (double[] x : scen.keySet())
			tree.add(x, new Cluster(x,scen.get(x)*sampleSize));
		double distance = 0.;
		for (double[] x : sample) {
			Cluster nearest = tree.get(x);
			distance += Metrics.squaredDistance(x, nearest.center);
		}
		return Math.sqrt(distance/sampleSize);
	}

	public static double getUpperBound(int dim, Map<double[],Double> scen, Collection<double[]> sample) {
		int sampleSize = sample.size();
		scengen.tools.MetricTree<Cluster> tree = new scengen.tools.MetricTree<>(dim);
		for (double[] x : scen.keySet())
			tree.add(x, new Cluster(x,scen.get(x)*sampleSize));
		double distance = 0.0;
		for (double[] x : sample) {
			Cluster next = tree.get(x);
			double dec1 = 1;
			double dec2 = next.decrement(dec1);
			while (dec1!=dec2) {
				distance += dec2*Metrics.squaredDistance(next.center,x);
				scengen.tools.MetricTree<Cluster> tree2 = new scengen.tools.MetricTree<>(dim);
				for (Cluster c : tree.getAll())
					if (c != next)
						tree2.add(c.center, new Cluster(c.center,c.weight));
				tree = tree2;
				if (tree.size()>0) {
					next = tree.get(x);
				}
				else {
					tree.add(next.center,next);
					next.weight = sampleSize;
				}
				dec1 -= dec2;
				dec2 = next.decrement(dec1);
			}
			distance += dec2*Metrics.squaredDistance(next.center,x);
		}
		return Math.sqrt(distance/sampleSize);
	}
	
	static class Cluster {
		double[] center;
		double weight;
		
		Cluster(double[] x, double w) {
			center = x;
			weight = w;
		}
		
		double decrement(double w) {
			double dec = Math.min(weight,w);
			weight -= w;
			return dec;
		}
		
	}
	
}
