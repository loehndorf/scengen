package main.java.scengen.methods;

import java.util.LinkedHashMap;
import java.util.Map;

import umontreal.ssj.hups.PointSetIterator;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import main.java.scengen.distribution.MultivariateDistribution;

public class QuasiMonteCarlo extends ReductionMethod {
	
	MultivariateDistribution _mvDist;
	static RandomStream rand = new MRG32k3a();
	
	public QuasiMonteCarlo (MultivariateDistribution mvDist) {
		_mvDist = mvDist;		
	}
	
	@Override
	public Map<double[], Double> getScenarios(int numScen) {
		int dim = _mvDist.getDim();
		SobolSequence sobol = new SobolSequence(numScen+1,dim);
		sobol.leftMatrixScramble(rand);
		PointSetIterator iter = sobol.iterator();
		iter.resetToNextPoint();
		Map<double[],Double> resultMap = new LinkedHashMap<>();
		for (int i=0; i<numScen; i++) {
			double[] u = new double[dim];
			iter.nextPoint(u, dim);
			resultMap.put(_mvDist.getRealization(u),1./numScen);
		}
		return resultMap;
	}
}
