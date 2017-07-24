package scengen.methods;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.LinkedHashMap;
import java.util.Map;

import scengen.distribution.DISTRIBUTION;
import scengen.distribution.MultivariateDistribution;
import scengen.distribution.MultivariateNormal;

public class QuantizationGrids extends ReductionMethod {

	MultivariateNormal _normal;
	int _dim;
	
	public QuantizationGrids (MultivariateDistribution mvDist) {
		if (mvDist.getDim()>10 || mvDist.getType()!=DISTRIBUTION.Normal)
			throw new IllegalArgumentException("VectorQuantization1 accepts only multivariate normal distributions with at most 10 dimensions.");
		_dim = mvDist.getDim();
		_normal = (MultivariateNormal) mvDist;
	}
	
	public Map<double[],Double> getScenarios(int numScen) {
		Map<double[],Double> map = new LinkedHashMap<>();
		try {
			String filename = String.format("%d_%d_%s",numScen,_dim,_dim>1?"nopti":"dualopti");
			InputStream stream = ClassLoader.getSystemClassLoader().getResourceAsStream(filename);
			BufferedReader br = new BufferedReader(new InputStreamReader(stream));
			String s;
			double[] values = null;
			int shift = _dim>1?1:0;
			while ((s = br.readLine())!=null) {
				String[] sar = s.split("\\s+");
				values = new double[_dim];
				for (int i=0; i<_dim; i++)
					values[i] = Double.parseDouble(sar[i+1+shift]);
				map.put(values,Double.parseDouble(sar[shift]));
			}
			//remove last line
			map.remove(values);
			br.close();
		} catch (IOException e) {e.printStackTrace();}
		return fitMultivariate(map);
	}
	
	public Map<double[],Double> fitMultivariate(Map<double[],Double> scen) {
		double[] mean = _normal.getMean();
		double[][] chol = _normal.getCholesky();
		Map<double[],Double> map = new LinkedHashMap<>();
		for (double[] d: scen.keySet()) {
			double[] x = new double[_dim];
			for (int i=0; i<_dim; i++) {
				double sum = 0.0;
				for (int j=0; j<i+1; j++)
					sum += d[j]*chol[i][j];
				x[i] = mean[i] + sum;
			}
			map.put(x,scen.get(d));
		}
		return map;
	}
	
	
	
	
}
