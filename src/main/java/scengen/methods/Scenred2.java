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

public class Scenred2 extends ReductionMethod {
	
	static String _sep = System.getProperty("file.separator");
	static String _scenredpath = System.getProperty("user.dir")+"/src/main/resources/scenred2/";
	static String _userdir = System.getProperty("user.dir")+_sep;
	
	static int _timeLimit = 600; //seconds
	static int _metricType = 3; //1:transport, 2:fortet-mourier ,3:wasserstein
	static int _reductionMethod = 1; //1:forward, 2:backward
	static int _pNorm = 2; //0: max, 1:sum, 2:euclidean
	static int _order = 2; //metric order 
	static int _dim;
	static int _size;
	Map<double[],Double> _input;
	
	public Scenred2(Map<double[],Double> input, int dim) {
		_dim = dim;
		_input = input;
	}
	
	@Override
	public synchronized Map<double[], Double> getScenarios(int numScen) {
		Map<double[],Double> map = null;
		StringBuilder sb = new StringBuilder();
		try {
			writeCommand();
			writeParam(numScen);
			writeSample();
			Process p = null;
			if(System.getProperty("os.name").startsWith("Windows"))
				p = Runtime.getRuntime().exec(_scenredpath+"scenred2.exe command.txt -nogams");
			else {
				throw new IllegalStateException("No compatible binary found for your OS.");
			}
			BufferedReader reader=new BufferedReader(new InputStreamReader(p.getInputStream())); 
			String line = reader.readLine();
			while(line!=null) { 
				sb.append(line);
				line=reader.readLine(); 
			}
			map = readOutput(numScen,_dim);
		} catch (IOException e) { 
			System.out.println(sb.toString()); 
			e.printStackTrace(); 
		}
		finally {
			new File(_userdir+"scenred2_in.txt").delete();
			new File(_userdir+"scenred2_out.txt").delete();
			new File(_userdir+"scenred2_param.txt").delete();
		}
		return map;
	}
	
	void writeSample() throws IOException {
		FileOutputStream stream = new FileOutputStream(_userdir+"scenred2_in.txt");
		BufferedWriter br = new BufferedWriter(new OutputStreamWriter(stream));
		br.write(String.format("TYPE FAN\n"
				+ "TIME 2\n"
				+ "SCEN %d\n"
				+ "RANDOM %d\n"
				+ "DATA\n",_input.size(),_dim));
		StringBuilder sb = new StringBuilder();
		for (double[] d : _input.keySet()) {
			sb.append(_input.get(d)+"\n");
			for (int j=0; j<_dim; j++)
				sb.append(1.0+"\t");
			sb.append("\n");
			for (int j=0; j<_dim; j++)
				sb.append(d[j]+"\t");
			sb.append("\n\n");
		}
		sb.append("END");
		br.write(sb.toString());
		br.close();
	}
	
	void writeCommand() throws IOException {
		FileOutputStream stream = new FileOutputStream(_userdir+"command.txt");
		BufferedWriter br = new BufferedWriter(new OutputStreamWriter(stream));
		br.write(String.format("report_level 1\n"
				+ "runtime_limit %d\n"
				+ "read_scen %sscenred2_in.txt\n"
				+ "scen_red %sscenred2_param.txt\n"
				+ "out_scen %sscenred2_out.txt",_timeLimit,_userdir,_userdir,_userdir));
		br.close();
	}
	
	void writeParam(int numScen) throws IOException {
		FileOutputStream stream = new FileOutputStream(_userdir+"scenred2_param.txt");
		BufferedWriter br = new BufferedWriter(new OutputStreamWriter(stream));
		br.write(String.format("red_num_leaves %d\n"
				+ "metric_type %d\n"
				+ "p_norm %d\n"
				+ "reduction_method %d\n"
				+ "order %d"
				+ "\n",numScen,_metricType,_pNorm,_reductionMethod,_order));
		br.close();
	}
	
	Map<double[], Double> readOutput(int numScen, int dim) throws IOException  {
		Map<double[],Double> map = new LinkedHashMap<>();
		FileInputStream stream = new FileInputStream(_userdir+"scenred2_out.txt");
		BufferedReader br = new BufferedReader(new InputStreamReader(stream));
//		for (int i=0; i<9; i++) br.readLine();
		String s;
		double[] values = null;
		while (!(s = br.readLine()).matches("DATA"));
		while ((s = br.readLine())!=null) {
			double prob = Double.parseDouble(s);
			s = br.readLine();
			s = br.readLine();
			String[] sar = s.split("\\s+");
			values = new double[dim];
			for (int i=0; i<dim; i++)
				values[i] = Double.parseDouble(sar[i]);
			map.put(values,prob);
			s = br.readLine();
			if (s.matches("END"))
				break;
		}
		br.close();
		return map;
	}
	
	

}
