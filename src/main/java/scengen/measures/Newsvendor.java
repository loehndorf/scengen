package scengen.measures;

import java.util.Iterator;
import java.util.SortedMap;
import java.util.TreeMap;

import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.distribution.LogNormalDistribution ;
import org.apache.commons.math3.distribution.NormalDistribution;



public abstract class Newsvendor {
	
	double _quantity;
	double _profit;
	double _fractile;
	double _safetyfactor;
	double _price;
	double _cost;
	
	private Newsvendor(double price, double cost) {
		if (price<cost) throw new IllegalArgumentException("Negative contribution margin.");
		_fractile = (price-cost)/price;
		_price = price;
		_cost = cost;
	};
	
	
	public double getFractile() {
		return _fractile;
	}
	
	public static Newsvendor solveNormalNewsvendor (final double price, final double cost, final double mu, final double sigma) {
		final NormalDistribution dist = new NormalDistribution();
		return new Newsvendor(price,cost) {{
			_safetyfactor = dist.inverseCumulativeProbability((price-cost)/price);
			_quantity = mu+sigma*_safetyfactor;
			_profit = (price-cost)*mu - price*sigma*dist.density(_safetyfactor);
		}
		@Override
		public double getProfit(double quantity) {
			double z = (quantity-mu)/sigma;
			double lostSales = sigma*(dist.density(z)-z*(1-dist.cumulativeProbability(z)));
			return _price*mu -_cost*quantity - _price*lostSales;
		}
		};
	}
	
	
	public static Newsvendor solveStudentNewsvendor (final double price, final double cost, final double mu, final double scale, final double deg, final double[] sample) {
		final TDistribution dist = new TDistribution(deg);
		return new Newsvendor(price,cost) {{
			_safetyfactor = dist.inverseCumulativeProbability((price-cost)/price);
			_quantity = mu+scale*Math.sqrt((deg/(deg-2)))*_safetyfactor;
			_profit = getProfit(_quantity);
		}
		@Override
		public double getProfit(double quantity) {
			double profit = -quantity*cost;
			for (int i=0; i<sample.length; i++) {
				double demand = sample[i];
				profit += Math.min(quantity,demand)*price/sample.length;
			}
			return profit;
		}};
	}
	
	public static Newsvendor solveUniformNewsvendor (final double price, final double cost, final double mu, final double sigma) {
		return new Newsvendor(price,cost) {
			double _min = mu-sigma*Math.sqrt(3);
			double _max = mu+sigma*Math.sqrt(3);
			{
			_quantity = _min+(_max-_min)*(price-cost)/price;
			_profit = getProfit(_quantity);
		}
		@Override
		public double getProfit(double quantity) {
			return -_cost*quantity + 0.5*price*(quantity*quantity-_min*_min)/(_max-_min)+price*quantity*(1-(quantity-_min)/(_max-_min));
		}
		};
	}
	
	public static Newsvendor solveLognormalNewsvendor (final double price, final double cost, final double mu, final double sigma) {
		final NormalDistribution dist1 = new NormalDistribution();
		final double cv = sigma/mu;
		final double nu = Math.log(mu)-Math.log(Math.sqrt(1+cv*cv));
		final double tau = Math.sqrt(Math.log(1+cv*cv));
		final LogNormalDistribution dist2 = new LogNormalDistribution(nu,tau);
		return new Newsvendor(price,cost) {{
			_safetyfactor = dist1.inverseCumulativeProbability((price-cost)/price);
			_quantity = Math.exp(nu+tau*_safetyfactor);
			_profit = (price-cost)*mu - price*mu*dist1.cumulativeProbability(tau-_safetyfactor)+cost*mu;
		}
		@Override
		public double getProfit(double quantity) {
			double lostSales = quantity*(1-dist2.cumulativeProbability(quantity))-Math.exp(nu+tau*tau/2)*dist1.cumulativeProbability((nu+tau*tau-Math.log(quantity))/tau);
			return _price*mu -_cost*quantity + _price*lostSales;
		}
		};
	}
	
	public static Newsvendor solveEmpiricalNewsvendor(final double price, final double cost, final double[] scenarios, final double[] weights) {
		return new Newsvendor(price,cost) {
			SortedMap<Double,Double> _sortedScen;
			{
			if (scenarios.length<1) throw new IllegalArgumentException("Number of scenarios less than one.");
			if (scenarios.length!=weights.length) throw new IllegalArgumentException("Number of scenarios does not match number of weights.");
			_sortedScen = new TreeMap<>();
			for (int i=0; i<scenarios.length; i++)
				_sortedScen.put(scenarios[i], weights[i]);
			double cum = 0.;
			Iterator<Double> iter = _sortedScen.keySet().iterator();
			while(iter.hasNext() && cum < _fractile) {
				Double demand = iter.next();
				cum += _sortedScen.get(demand);
				_quantity = demand;
			}
			_profit = getProfit(_quantity);
		}
		@Override
		public double getProfit(double quantity) {
			double profit = -quantity*cost;
			for (Double demand : _sortedScen.keySet())
				profit += _sortedScen.get(demand)*Math.min(quantity,demand)*price;
			return profit;
		}
		};
	}
	
	public double getOptimalQuantity() {
		return _quantity;
	}
	
	public double getOptimalProfit() {
		return _profit;
	}
	
	public abstract double getProfit(double quantity);
	
}
