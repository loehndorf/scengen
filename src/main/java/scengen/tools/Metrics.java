package scengen.tools;

public final class Metrics {
	
	public static final double FortetMourier(double[] d1, double[] d2) {
        double x1 = 0;
        double x2 = 0;
    	double d = 0;
        for (int i = 0; i < d1.length; i++) {
            double diff = (d1[i] - d2[i]);
            if (!Double.isNaN(diff)) {
                d += diff * diff;
            }
            x1 += d1[i]*d1[i];
            x2 += d2[i]*d2[i];
        }
        return Math.sqrt(d)*Math.max(1,Math.max(Math.sqrt(x1),Math.sqrt(x2)));
    }
	
	public static final double WeightedEuclidean(double[] weights, double[] d1, double[] d2) {
		if (d1==d2) return 0.;
		double sumSq = 0.;
		for (int i=0; i<d1.length; i++) {
			double d = d1[i]-d2[i];
			sumSq += d*d;
		}
		return Math.sqrt(sumSq);
	}
	
	public static final double Euclidean(double[] d1, double[] d2) {
		if (d1==d2) return 0.;
		double sumSq = 0.;
		for (int i=0; i<d1.length; i++) {
			double d = d1[i]-d2[i];
			sumSq += d*d;
		}
		return Math.sqrt(sumSq);
	}
	
	public static final double squaredDistance(double[] d1, double[] d2) {
		if (d1==d2) return 0.;
		double sumSq = 0.;
		for (int i=0; i<d1.length; i++) {
			double d = d1[i]-d2[i];
			sumSq += d*d;
		}
		return sumSq;
	}
	
	public static final double manhattanDistance(double[] d1, double[] d2) {
		if (d1==d2) return 0.;
		double sumSq = 0.;
		for (int i=0; i<d1.length; i++) {
			double d = Math.abs(d1[i]-d2[i]);
			sumSq += d;
		}
		return sumSq;
	}

	public static final double squaredDistance(double[] d1, double[] d2, double[] weights) {
		if (d1==d2) return 0.;
		double sumSq = 0.;
		for (int i=0; i<d1.length; i++) {
			double d = d1[i]-d2[i];
			sumSq += d*d*weights[i];
		}
		return sumSq;
	}
	
	public static final double infDistance(double[] d1, double[] d2, double[] weights) {
		if (d1==d2) return 0.;
		double maxDist = 0.;
		for (int i=0; i<d1.length; i++) {
			double d = Math.abs(d1[i]-d2[i])*weights[i];
			if (d>maxDist)
				maxDist = d;
		}
		return maxDist;
	}

}
