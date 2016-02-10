package main.java.scengen.tools;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
//import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

public class Matrix extends Array2DRowRealMatrix {

	private static final long serialVersionUID = 4469929056984141554L;

	public Matrix(double[][] rawData) throws DimensionMismatchException,
			NotStrictlyPositiveException {
		super(rawData);
	}
	
	public Matrix(int rows, int columns) throws DimensionMismatchException,
	NotStrictlyPositiveException {
		super(rows,columns);
	}
	
	/**
	 * 
	 * @return a new matrix of the element-wise logarithm of this matrix
	 */
	public Matrix log() {
		Matrix out = new Matrix(getRowDimension(), getColumnDimension());
		for (int i=0; i<getRowDimension(); i++) 
			for (int j=0; j<getColumnDimension(); j++) {
				double d = getEntry(i,j);
				d = Math.log(d);
				out.setEntry(i, j, d);
			}
		return out;	
	}
	
	/**
	 * @return a new matrix of the element-wise exponentiation of this matrix
	 */
	public Matrix exp() {
		Matrix out = new Matrix(getRowDimension(), getColumnDimension());
		for (int i=0; i<getRowDimension(); i++) 
			for (int j=0; j<getColumnDimension(); j++) {
				double d = getEntry(i,j);
				d = Math.exp(d);
				out.setEntry(i, j, d);
			}
		return out;	
	}
	
	/**
	 * @return a new matrix of the element-wise inverse hyperbolic sine of this matrix
	 */
	public Matrix arcsinh() {
		Matrix out = new Matrix(getRowDimension(), getColumnDimension());
		for (int i=0; i<getRowDimension(); i++) 
			for (int j=0; j<getColumnDimension(); j++) {
				double d = getEntry(i,j);
				d = Math.log(d+Math.sqrt(d*d + 1.));
				out.setEntry(i, j, d);
			}
		return out;	
	}
	
	/**
	 * @return a new matrix of the element-wise hyperbolic sine of this matrix
	 */
	public Matrix sinh() {
		Matrix out = new Matrix(getRowDimension(), getColumnDimension());
		for (int i=0; i<getRowDimension(); i++) 
			for (int j=0; j<getColumnDimension(); j++) {
				double d = getEntry(i,j);
				d = Math.sinh(d);
				out.setEntry(i, j, d);
			}
		return out;	
	}
	
	/**
	 * @param p
	 * @return the p-th root matrix of this matrix
	 */
	public Matrix root (double p) {
		if (!isSquare())
			throw new IllegalArgumentException("No sqaure matrix.");
		if (p==1) return this;
		if (p<1) throw new IllegalArgumentException(String.format("Illegal fractional root: %f", p));
		EigenDecomposition ev = new EigenDecomposition(this);
		RealMatrix V = ev.getV();
		RealMatrix D = ev.getD();
		for (int i=0; i<D.getRowDimension(); i++) {
			if (D.getEntry(i,i)<0)
				throw new IllegalStateException("Root of the matrix is not defined.");
			D.setEntry(i,i,Math.pow(D.getEntry(i,i),1./p));
		}
		RealMatrix B = V.multiply(D).multiply(new LUDecomposition(V).getSolver().getInverse());
		return new Matrix(B.getData());
	}
	
	/**
	 * 
	 * @param other another Matrix
	 * @return a new matrix with {@code Double.NaN} replaced by the valze of the other dataframe.
	 */
	public Matrix replaceNaN(RealMatrix other) {
		Matrix out = new Matrix(getRowDimension(), getColumnDimension());
		for (int i=0; i<getRowDimension(); i++) 
			for (int j=0; j<getColumnDimension(); j++) {
				double d = getEntry(i,j);
				if (Double.isNaN(d))
					d = other.getEntry(i,j);
				out.setEntry(i, j, d);
			}
		return out;	
	}
	
	/**
	 * Print tab-separated matrix into the console.
	 */
	public void print() {
		for (int i=0; i<getRowDimension(); i++) {
			System.out.print(getEntry(i,0));
			for (int j=1; j<getColumnDimension(); j++)
				System.out.print("\t"+getEntry(i,j));
			System.out.println();
		}
	}
	
	/**
	 * Return matrix as a text string, with entries separated by commas and rows in curly braces.
	 */
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("{");
		for (int i=0; i<getRowDimension(); i++) {
			sb.append("{").append(getEntry(i,0));
			for (int j=1; j<getColumnDimension(); j++)
				sb.append(",").append(getEntry(i,j));
			if (i<getRowDimension()-1)
				sb.append("}\n,");
		}
		sb.append("}}");
		return sb.toString();
	}

}
