package main.java.scengen.methods;

import java.util.LinkedHashMap;
import java.util.Map;

import main.java.scengen.distribution.MultivariateDistribution;

public class MonteCarlo extends ReductionMethod {
	
	MultivariateDistribution _mvDist;
	
	public MonteCarlo (MultivariateDistribution mvDist) {
		_mvDist = mvDist;
		
	}
	
	@Override
	public Map<double[], Double> getScenarios(int numScen) {
		Map<double[],Double> resultMap = new LinkedHashMap<>();
		for (int i=0; i<numScen; i++)
			resultMap.put(_mvDist.getRealization(),1./numScen);
		return resultMap;
	}
}
