package main.java.scengen.tools.kdtree;

import java.util.HashMap;

/**
 * A SparseArray is a multidimensional (sparse) array based on HashMap. As with arrays, the user has to provide the size in each dimension. Unlike
 * with arrays, heap space is not reserved in advance. The data structure relies on varargs. Like an arrays, it throws an ArrayIndexOutOfBoundException
 * if the index exceeds the bounds of the array. To specify an index, MultiArray makes use of varagrs. For example, if we want to insert an object 'obj' into the array
 * at position (1,2,3,4), we simply call put(obj,1,2,3,4). To retrieve the same object, we call get(1,2,3,4). 
 * 
 * @author Nils
 *
 * @param <V>
 */
public class SparseArray<V> extends HashMap<Integer,V> {
	
	private static final long serialVersionUID = 6642522335359990799L;
	int _dim;
	int[] _cumSize;
	int[] _size;
	
	/**
	 * Create a sparse array with a given dimensionality.
	 * @param dimensions maximum index in each dimension
	 */
	public SparseArray(Integer... dimensions) {
		_dim = dimensions.length;
		_size = new int[_dim];
		_cumSize = new int[_dim];
		for (int i=_dim-1; i>=0; i--) {
			_size[i] = dimensions[i];
			_cumSize[i] = i==_dim-1?1:dimensions[i+1]*_cumSize[i+1];
		}
	}
	
	int key (Integer... index) {
		int key = 0;
		for (int i=0; i<_dim; i++) {
			if (index[i] >= _size[i])
				throw new IndexOutOfBoundsException();
			key += index[i]*_cumSize[i];
		}
		return key;
	}
	
	@Override @Deprecated 
	public V put(Integer i, V value) {
		System.err.println("ERROR: Use put(V value, Integer... index) to insert values into this map.");
		throw new UnsupportedOperationException();
	}
	
	public V put(V value, Integer... index) {
		return super.put(key(index),value);
	}
	
	/**
	 * Adds all entries of m into this array. 
	 * @param m another {@code SparseArrays}
	 * @throws Throws an {@link IndexOutOfBoundsException} if the dimensionality of m does not match that of this array.
	 */
	public void putAll(SparseArray<V> m) {
		if (equalDimensions(m))
			for (Integer i : m.keySet())
				super.put(i, m.get(i));
		else
			throw new IndexOutOfBoundsException();
	}
	
	public V get(Integer... index) {
		return super.get(key(index));
	}
	
	public boolean contains(Integer... index) {
		return super.containsKey(key(index));
	}
	
	public V remove(Integer... index) {
		return super.remove(key(index));
	}
	
	public boolean equalDimensions(SparseArray<V> m) {
		if (_dim!=m._dim)
			return false;
		for (int i=0; i<_dim; i++)
			if(_size[i]!=m._size[i])
				return false;
		return true;
	}
	
	public static boolean test() {
		SparseArray<Integer> map = new SparseArray<>(10,24,7);
		int i=9;
		int j=23;
		int k=6;
		int key = 24*7*i+7*j+k;
		map.put(key,i,j,k);
		if (map.get(i,j,k)==key)
			return true;
		return false;
	}

}
