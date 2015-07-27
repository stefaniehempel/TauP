package edu.sc.seis.TauP;

import static org.junit.Assert.*;

import org.junit.Test;


public class ShootTest {

    @Test
    public void testShootExistingRayParam() throws Exception {
        TauModelLoader.clearCache();
        double depth = 119;
        SeismicPhase phase = new SeismicPhase("P", "iasp91", depth);
        for (int i = 0; i < phase.getRayParams().length; i++) {
            Arrival maxRPArrival = phase.shootRay(phase.getRayParams()[i]);
            assertEquals(i+"th ray param dist", phase.getDist()[i], maxRPArrival.getDist(), 0.0001);
            assertEquals(i+"th ray param time", phase.getTime()[i], maxRPArrival.getTime(), 0.0001);
        }
    }
    
    @Test
    public void testShootMiddleRayParam() throws Exception {
        TauModelLoader.clearCache();
        double depth = 119;
        SeismicPhase phase = new SeismicPhase("P", "iasp91", depth);
        for (int i = 0; i < phase.getRayParams().length-1; i++) {
            double rp = (phase.getRayParams()[i]+phase.getRayParams()[i+1])/2;
            double timeTol = Math.abs(phase.getTime()[i]-phase.getTime()[i+1]);
            Arrival maxRPArrival = phase.shootRay(rp);
            assertEquals(i+"th ray param dist", phase.getDist()[i], maxRPArrival.getDist(), 0.1);
            assertEquals(i+"th ray param time, neighbors: "+phase.getTime()[i]+" "+phase.getTime()[i+1], phase.getTime()[i], maxRPArrival.getTime(), timeTol);
            assertEquals(i+"th+1 ray param dist", phase.getDist()[i+1], maxRPArrival.getDist(), 0.1);
            assertEquals(i+"th+1 ray param time", phase.getTime()[i+1], maxRPArrival.getTime(), timeTol);
        }
    }
}