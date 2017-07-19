package scengen.application;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

public class GeneratorTest {
    @Before
    public void setUp() throws Exception {
    }

    @After
    public void tearDown() throws Exception {
    }

    @Test
    public void testDistNormal() {
        scengen.application.Generator.main("dist=Normal",
                "dim=2",
                "scen=10",
                "mean=0,0",
                "cov=1,0,0,1",
                "method=MonteCarlo");
    }

    @Test
    public void testDistLognormal() {
        scengen.application.Generator.main("dist=Lognormal",
                "dim=2",
                "scen=10",
                "mean=0,0",
                "cov=1,0,0,1",
                "method=MonteCarlo");
    }

    @Test
    public void testDistUniform1() {
        scengen.application.Generator.main("dist=Uniform",
                "dim=2",
                "scen=10",
                "mean=0,0",
                "cov=1,0,0,1",
                "method=MonteCarlo");
    }

    @Test
    public void testDistUniform2() {
        scengen.application.Generator.main("dist=Uniform",
                "dim=2",
                "scen=10",
                "min=-2,-2",
                "max=2,2",
                "correl=1,0,0,1",
                "method=MonteCarlo");
    }

    @Test
    public void testDistStudent() {
        scengen.application.Generator.main("dist=Student",
                "dim=2",
                "scen=10",
                "mean=0,0",
                "cov=1,0,0,1",
                "df=5,5",
                "method=QuantizationLearning");
    }

    @Test
    public void testMethodQMC() {
        scengen.application.Generator.main("dist=Normal",
                "dim=2",
                "scen=10",
                "mean=0,0",
                "cov=1,0,0,1",
                "method=QuasiMonteCarlo");
    }

//    @Test
//    public void testMethodMomentMatching() {
//        scengen.application.Generator.main("dist=Normal",
//                "dim=2",
//                "scen=10",
//                "mean=0,0",
//                "cov=1,0,0,1",
//                "method=MomentMatching");
//    }

    @Test
    public void testMethodQuantization1() {
        scengen.application.Generator.main("dist=Normal",
                "dim=2",
                "scen=10",
                "mean=0,0",
                "cov=1,0,0,1",
                "method=QuantizationLearning");
    }

//    @Test
//    public void testMethodQuantization2() {
//        scengen.application.Generator.main("dist=Normal",
//                "dim=2",
//                "scen=10",
//                "mean=0,0",
//                "cov=1,0,0,1",
//                "method=QuantizationGrids");
//    }

    @Test
    public void testMethodVoronoiCellSampling() {
        scengen.application.Generator.main("dist=Normal",
                "dim=2",
                "scen=10",
                "mean=0,0",
                "cov=1,0,0,1",
                "method=VoronoiCellSampling");
    }

}