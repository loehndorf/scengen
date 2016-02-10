package main.java.scengen.tools.kdtree;

/**
 *
 */
class CityBlock implements DistanceFunction {
    @Override
    public double distance(double[] p1, double[] p2) {
        double d = 0;

        for (int i = 0; i < p1.length; i++) {
            d += Math.abs(p1[i] - p2[i]);
        }

        return d;
    }

    @Override
    public double distanceToRect(double[] point, double[] min, double[] max) {
        double d = 0;

        for (int i = 0; i < point.length; i++) {
            double diff = 0;
            if (point[i] > max[i]) {
                diff = Math.abs(point[i] - max[i]);
            }
            else if (point[i] < min[i]) {
                diff = Math.abs(point[i] - min[i]);
            }
            d += diff;
        }

        return d;
    }
}