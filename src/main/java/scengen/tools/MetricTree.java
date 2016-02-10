package main.java.scengen.tools;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.Map;

import main.java.scengen.tools.kdtree.KdTree;
import main.java.scengen.tools.kdtree.SquaredEuclidean;


/**
 * Wrapper class for a metric tree to facilitate nearest neighbor queries over a fixed set of points. 
 * @author Nils Loehndorf
 *
 * @param <E>
 */
public class MetricTree<E> implements Serializable {

	private static final long serialVersionUID = 5100441495706307337L;

	KdTree<E> kdtree;
	SquaredEuclidean dist;
	Map<double[],E> map = new LinkedHashMap<>();
	int _dimension;
	double[] _weights;
	
	public MetricTree(int dimension) {
		_dimension = dimension;
		_weights = new double[dimension];
		Arrays.fill(_weights,1);
		kdtree = new KdTree<E>(dimension);
		dist = new SquaredEuclidean();
	}
	
	/**
	 * Add the given element to the tree at the given location.
	 * @param point
	 * @param e
	 */
	public void add(double[] point, E e) {
		map.put(point,e);
		kdtree.addPoint(point, e);
	}
	
	/**
	 * Retrieve the element closest to the given location from the tree.
	 * @param point
	 * @return Return the nearest element.
	 */
	public E get(double[] point) {
		return kdtree.findNearestNeighbors(point, 1, dist).getMax();

	}
	
	/**
	 * Checks whether an element can be found that has the given point as its key. Compares by references not by value.
	 * @param point
	 * @return
	 */
	public boolean contains(double[] point) {
		if (map.containsKey(point))
			return true;
		return false;
	}
	
	/**
	 * Retrieves all elements that are stored with this tree.
	 * @return
	 */
	public Collection<E> getAll() {
		return map.values();
	}
	
	/**
	 * @return Returns the size of the tree.
	 */
	public int size() {
		return map.size();
	}
	

}
