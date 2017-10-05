package edu.sc.seis.TauP;

public class Amplitude {

	private SeismicPhase thisPhase;
	
	private double amplitudeFactor;
	
	private double RTproduct;
	
	private boolean criticalOnly = false;
	
	private boolean polarization = false;
	
	private boolean horizontal;
	
	public Amplitude(SeismicPhase thisPhase) {
		this.thisPhase = thisPhase;
	}
	
	public void setCritical(boolean criticalOnly) {
		this.criticalOnly = criticalOnly;
	}
	
	public void setPolarizationControl(boolean pol) {
		setPolarizationControl(pol, true);
	}
	
	public void setPolarizationControl(boolean pol, boolean horizontal) {
		this.polarization = true;
		this.horizontal = horizontal;
	}
	
	public void setPolarization(boolean horizontal) {
		this.horizontal = horizontal;
	}
	
	// need to compute geometrical spreading
	private void computeAmplitudeFactor(Arrival interpArrival) throws NoSuchLayerException, NoSuchMatPropException {
		double geomFact = 0.;
		double takeoffVelocity, takeoffAngle, incidentAngle, incidentVelocity;
		char lastLeg;
		if (thisPhase.getName().endsWith("kmps")) {
	        geomFact = 0.;
	        amplitudeFactor = 0.;
	    } else {
	    	// values to be obtained from SeismicPhase
	    	VelocityModel vMod = thisPhase.getTauModel().getVelocityModel();
	    	double sourceDepth = thisPhase.sourceDepth;
	    	double receiverDepth = thisPhase.receiverDepth;
	    	String name = thisPhase.name;
	    	double radiusOfEarth = thisPhase.getTauModel().getRadiusOfEarth();
	    	
	    	// keeping radians for a while longer
	    	takeoffVelocity = (thisPhase.getDownGoing()[0] ? 1 : -1 ) *
					vMod.evaluateBelow(sourceDepth, name.charAt(0));
            takeoffAngle = thisPhase.calcTakeoffAngle(interpArrival.getRayParam())/180.*Math.PI;
	        lastLeg = thisPhase.getLegs().get(thisPhase.getLegs().size()-2).charAt(0);
	        incidentVelocity = vMod.evaluateBelow(receiverDepth, lastLeg);
	        incidentAngle = thisPhase.calcIncidentAngle(interpArrival.getRayParam())/180*Math.PI;
	        //Math.asin(incidentVelocity*interpArrival.getRayParam()/radiusOfEarth);
	        
	        // geometrical spreading and amplitude factors
	        if (interpArrival.getDddp() != 0.) {
	        	geomFact = Math.sin(takeoffAngle);
	        	if (geomFact != 0.) {
	        		geomFact = Math.cos(incidentAngle)*Math.cos(takeoffAngle)*
	                		Math.sin(interpArrival.getDist()) * interpArrival.getDddp() /
	                		interpArrival.getRayParam();
	                geomFact = (radiusOfEarth - receiverDepth)*(radiusOfEarth-sourceDepth)
	                		/takeoffVelocity*1e3 //considering R in km and velocity in km/s > *1e3*1e3/1e3 to get to SI-units
	                		*Math.sqrt(Math.abs(geomFact)); //unit: m*s
	                
	                amplitudeFactor = Math.pow(takeoffVelocity,5)*incidentVelocity*
	                		vMod.evaluateBelow(sourceDepth, 'r') * //density at source, in g/cm3
	        				vMod.evaluateBelow(receiverDepth, 'r')*1.e24;
	                amplitudeFactor = 1/(4.*Math.PI*Math.signum(amplitudeFactor)*Math.sqrt(Math.abs(amplitudeFactor)));
	                amplitudeFactor /= geomFact;//*Math.cos(incidentAngle); deprecated definition of amplFact
	                		//amplfact unit kg m/s^2
	        	}
	        }
	    }
	}
	
	public void compute(Arrival interpArrival) throws NoSuchLayerException, NoSuchMatPropException {
		computeAmplitudeFactor(interpArrival);
		// computeRTfactor(interpArrival);
	}
	
	public double getAmplFact() {
		return amplitudeFactor;
	}
	
	public double getRT() {
		return RTproduct;
	}
	
}
