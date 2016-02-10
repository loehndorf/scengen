package main.java.scengen.methods;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.LinkedHashMap;
import java.util.Map;

import main.java.scengen.tools.Xorshift;

public class MomentMatching extends ReductionMethod {
	
	static String _sep = System.getProperty("file.separator");
	static String _mmpath = System.getProperty("user.dir")+"/bin/";
	static String _userdir = System.getProperty("user.dir")+_sep;
	long _time = System.nanoTime();
	int _dim;
	double[] _mean;
	double[] _stdev;
	double[] _skewness;
	double[] _kurtosis;
	double[][] _correl;
	boolean _debugmode = false;
	long _seed;
	
	public MomentMatching(double[] mean, double[][] covariance, long seed) {
		_dim = mean.length;
		if (covariance.length != mean.length || covariance[0].length != mean.length)
	         throw new IllegalArgumentException("Covariance matrix dimensions invalid.");
		_mean = mean;
		_skewness = new double[_dim];
		_kurtosis = new double[_dim];
		_correl = correlation(mean,covariance);
		_stdev = stdev(covariance);
		_seed = seed;
	}
	
	public MomentMatching(double[] mean, double[][] covariance, double[] skewness, double[] kurtosis, long seed) {
		_dim = mean.length;
		if (covariance.length != mean.length || covariance[0].length != mean.length)
	         throw new IllegalArgumentException("Covariance matrix dimensions invalid.");
		_mean = mean;
		_skewness = skewness;
		_kurtosis = kurtosis;
		_correl = correlation(mean,covariance);
		_stdev = stdev(covariance);
		_seed = seed;
	}
	
	static double[][] correlation (double[] mean, double[][] covariance) {
		double[][] correl = new double[covariance.length][covariance.length];
		for (int i=0; i<covariance.length; i++) {
			correl[i][i] = 1.;
			for (int j=i+1; j<covariance.length; j++) {
				correl[i][j] = covariance[i][j]/(Math.sqrt(covariance[i][i])*Math.sqrt(covariance[j][j]));
				correl[j][i] = correl[i][j];
			}
		}	
		return correl;
	}
	
	static double[] stdev (double[][] covariance) {
		double[] stdev = new double[covariance.length];
		for (int i=0; i<covariance.length; i++)
			stdev[i] = Math.sqrt(covariance[i][i]);
		return stdev;
	}
	
	@Override
	public synchronized Map<double[], Double> getScenarios(int numScen) {
		Map<double[],Double> map = null;
		StringBuilder sb = new StringBuilder();
		int seed = new Xorshift(_seed).nextInt();
		try {
			writeMoments(_dim);
			writeCorrelation(_dim);
			String filename = null;
			if(System.getProperty("os.name").startsWith("Mac")){
				filename = "scengen_HKW.mac";
			}
			else if(System.getProperty("os.name").startsWith("Linux")){
				filename = "scengen_HKW.linux";
			}
			else if(System.getProperty("os.name").startsWith("Windows")){
				filename = "scen-gen_HKW.exe";
			}
			else if(System.getProperty("os.name").startsWith("Microsoft")){
				filename = "scen-gen_HKW.exe";
			}
			else {
				throw new IllegalStateException("No compatible binary found for your OS.");
			}
			Process p = Runtime.getRuntime().exec(_mmpath+filename+" "+numScen+" -f 1 -r "+seed);
			BufferedReader reader=new BufferedReader(new InputStreamReader(p.getInputStream())); 
			String line = reader.readLine();
			while(line!=null) { 
				sb.append(line);
				line=reader.readLine(); 
			}
			map = readOutput(numScen,_dim);
		} catch (IOException e) { 
			System.out.println(sb.toString()); 
			e.printStackTrace(); }
		finally {
			new File(_userdir+"out_scen.txt").delete();
			new File(_userdir+"tg_moms.txt").delete();
			new File(_userdir+"tg_corrs.txt").delete();
		}
		return map;
	}
	
	void writeMoments(int dim) throws IOException {
		FileOutputStream stream = new FileOutputStream(_userdir+"tg_moms.txt");
		BufferedWriter br = new BufferedWriter(new OutputStreamWriter(stream));
		br.write("4\n");
		br.write(dim+"\n\n");
		StringBuilder sb = new StringBuilder();
		for (int i=0; i<dim; i++)
			sb.append(_mean[i]+"\t");
		sb.append("\n");
		for (int i=0; i<dim; i++)
			sb.append(_stdev[i]+"\t");
		sb.append("\n");
		for (int i=0; i<dim; i++)
			sb.append(_skewness[i]+"\t");
		sb.append("\n");
		for (int i=0; i<dim; i++)
			sb.append(_kurtosis[i]+"\t");
		br.write(sb.toString());
		br.close();
	}
	
	void writeCorrelation(int dim) throws IOException {
		FileOutputStream stream = new FileOutputStream(_userdir+"tg_corrs.txt");
		BufferedWriter br = new BufferedWriter(new OutputStreamWriter(stream));
		br.write(dim+"\n");
		br.write(dim+"\n");
		StringBuilder sb = new StringBuilder();
		for (int i=0; i<dim; i++) {
			sb.append("\n");
			for (int j=0; j<dim; j++)
				sb.append(_correl[i][j]+"\t");
		}
		br.write(sb.toString());
		br.close();
	}
	
	Map<double[], Double> readOutput(int numScen, int dim) throws IOException  {
		Map<double[],Double> map = new LinkedHashMap<>();
		FileInputStream stream = new FileInputStream(_userdir+"out_scen.txt");
		BufferedReader br = new BufferedReader(new InputStreamReader(stream));
		for (int i=0; i<4; i++) br.readLine();
		String s;
		double[] values = null;
		while ((s = br.readLine())!=null) {
			String[] sar = s.split("\\s+");
//			System.out.println(Arrays.toString(sar));
			values = new double[dim];
			for (int i=0; i<dim; i++)
				values[i] = Double.parseDouble(sar[i+1]);
			map.put(values,Double.parseDouble(sar[0]));
		}
		br.close();
		return map;
	}
	
	

}
