package test.java.scengen.application;

import org.junit.Test;


public class Generator {

	@Test
	public void testDistNormal() {
		main.java.scengen.application.Generator.main("dist=Normal",
				"dim=2",
				"scen=10",
				"mean=0,0",
				"cov=1,0,0,1",
				"method=MonteCarlo");
	}
	
	@Test
	public void testDistLognormal() {
		main.java.scengen.application.Generator.main("dist=Lognormal",
				"dim=2",
				"scen=10",
				"mean=0,0",
				"cov=1,0,0,1",
				"method=MonteCarlo");
	}
	
	@Test
	public void testDistUniform1() {
		main.java.scengen.application.Generator.main("dist=Uniform",
				"dim=2",
				"scen=10",
				"mean=0,0",
				"cov=1,0,0,1",
				"method=MonteCarlo");
	}
	
	@Test
	public void testDistUniform2() {
		main.java.scengen.application.Generator.main("dist=Uniform",
				"dim=2",
				"scen=10",
				"min=-2,-2",
				"max=2,2",
				"correl=1,0,0,1",
				"method=MonteCarlo");
	}
	
	@Test
	public void testDistStudent() {
		main.java.scengen.application.Generator.main("dist=Student",
				"dim=2",
				"scen=10",
				"mean=0,0",
				"cov=1,0,0,1",
				"df=5,5",
				"method=QuantizationLearning");
	}
	
	@Test
	public void testMethodQMC() {
		main.java.scengen.application.Generator.main("dist=Normal",
				"dim=2",
				"scen=10",
				"mean=0,0",
				"cov=1,0,0,1",
				"method=QuasiMonteCarlo");
	}
	
	@Test
	public void testMethodMomentMatching() {
		main.java.scengen.application.Generator.main("dist=Normal",
				"dim=2",
				"scen=10",
				"mean=0,0",
				"cov=1,0,0,1",
				"method=MomentMatching");
	}
	
	@Test
	public void testMethodQuantization1() {
		main.java.scengen.application.Generator.main("dist=Normal",
				"dim=2",
				"scen=10",
				"mean=0,0",
				"cov=1,0,0,1",
				"method=QuantizationLearning");
	}
	
	@Test
	public void testMethodQuantization2() {
		main.java.scengen.application.Generator.main("dist=Normal",
				"dim=2",
				"scen=10",
				"mean=0,0",
				"cov=1,0,0,1",
				"method=QuantizationGrids");
	}
	
	@Test
	public void testMethodVoronoiCellSampling() {
		main.java.scengen.application.Generator.main("dist=Normal",
				"dim=2",
				"scen=10",
				"mean=0,0",
				"cov=1,0,0,1",
				"method=VoronoiCellSampling");
	}

}
