package edu.sc.seis.TauP;

import java.util.ArrayList;
import java.util.List;

public class Amplitude {

	private SeismicPhase thisPhase;
	
	private double amplitudeFactor;
	
	private double RTproduct;
	
	private boolean criticalOnly = TauP_Time.amplCriticalRTOnly; //default false
	
	private boolean polarization = TauP_Time.amplPolarizationSwitch; //default false
	
	private boolean horizontal = TauP_Time.amplHorizontal; //default false
	
	private boolean freeSurfaceSwitch = TauP_Time.freeSurfaceSwitch; //default true
	
	private boolean verticalComp = TauP_Time.verticalComp; //default true
	
	private boolean compControl = TauP_Time.compControl; //default false
	
	public Amplitude(SeismicPhase thisPhase) {
		this.thisPhase = thisPhase;
	}
	
	/*
	 * computing the product of reflection/ transmission coefficients and the free surface correction
	 */
	private void computeRTfactor(double arrivalRayParam) throws NoSuchLayerException, NoSuchMatPropException {
    	//computing RT for leaving the branch/ leg 
    	double inVp, inVs, inDensity, outVp, outVs, outDensity;
    	double interfaceR;
    	double cTopDepth, interfaceDepth, cBotDepth;
    	
    	// get info from current phase
    	TauModel tMod = thisPhase.getTauModel();
    	VelocityModel vMod=tMod.getVelocityModel();
    	List<Integer> branchSeq = thisPhase.branchSeq;
    	List<Boolean> waveType = thisPhase.waveType;
    	List<Integer> legAction = thisPhase.legAction;
    	List<Boolean> downGoing = thisPhase.downGoing;
    	double planetRadius=vMod.getRadiusOfEarth();
    	
    	// automatically assign vertical if not ScS, if polarization control not chosen by user
    	if ( !(thisPhase.name.contains("K") &&
    			( thisPhase.name.contains("ScS") || thisPhase.name.contains("SS") ||
					thisPhase.name.contains("sS") || thisPhase.name.contains("Ss") ) ) ) {
    			if (!polarization) horizontal=true;
    			if (!compControl) verticalComp=false;
    	}
    	
    	// use new RTcoeff for now
    	RTcoeff coeff;
    	boolean waveIn, waveOut = true, legActionChange, topFluidLayer, botFluidLayer;
    	int RTaction;
    	RTproduct=1.; // initialize output
    	
		for (int j=0; j < branchSeq.size()-1; j++) {
    		cTopDepth=tMod.getTauBranch(branchSeq.get(j), true).getTopDepth();
    		cBotDepth=tMod.getTauBranch(branchSeq.get(j), true).getBotDepth();
    		
    		if (downGoing.get(j)) { //RT at bottom of branch 
    			interfaceDepth = cBotDepth;
    			RTaction = 4; //transmission downwards by default
    		} else { //RT at top of branch
    			interfaceDepth = cTopDepth;
    			RTaction = 3; //transmission upwards by default
    		}
    		
			interfaceR=planetRadius-interfaceDepth;
			
			if (interfaceR > 0.) {
				
				// if change in legAction to last leg-action then deal with it, else transmission
	    		// always dealing with action at the end of the branch
	    		waveIn=waveType.get(j);
	    		waveOut=waveType.get(j+1);
	    		legActionChange= ( Math.abs(legAction.get(j)-legAction.get(j+1)) > 0 );
			
				if (interfaceDepth > 0. || downGoing.get(j)) {
	        		inVp=vMod.evaluateAbove(interfaceDepth, 'P');
	        		inVs=vMod.evaluateAbove(interfaceDepth, 'S');
	        		inDensity=vMod.evaluateAbove(interfaceDepth, 'r');
				}
				else {
					inVp = 0; inVs = 0; inDensity = 0;
				}
				
				outVp=vMod.evaluateBelow(interfaceDepth, 'P');
				outVs=vMod.evaluateBelow(interfaceDepth, 'S');
				outDensity=vMod.evaluateBelow(interfaceDepth, 'r');
				
	    		if (legActionChange) {
	    			RTaction = legAction.get(j);
	    		}
	    		
	    		topFluidLayer = ( tMod.getSlownessModel().depthInFluid(interfaceDepth*.99) ) ;
	    		botFluidLayer = ( tMod.getSlownessModel().depthInFluid(interfaceDepth) );
	    		
	    		// assigning K!
	    		if (botFluidLayer && !downGoing.get(j) || //upgoing wave in outer core
	    				topFluidLayer && downGoing.get(j)) { //downgoing wave in outer core
	    			waveIn = true;
	    		}
	    		if (downGoing.get(j)) {
	    			if ((topFluidLayer && RTaction<3) || //reflection at ICB from above
	    					botFluidLayer && RTaction>2 ) { //transmission into outer core from mantle
	    				waveOut = true;
	    			}
	    		}
	    		else {
	    			if ((botFluidLayer && RTaction<3) || //reflection at CMB from below
	    					topFluidLayer && RTaction>2 ) { //transmission into outer core from inner core
	    				waveOut = true;
	    			}
	    		}
	    			
	    		coeff =  new RTcoeff(inVp, inVs, inDensity,
	                    outVp, outVs, outDensity);
	
				if ( !( downGoing.get(j) ) ) coeff=coeff.flip();
				//if ( coeff.getTopVs()==0 ) coeff=coeff.flip(); // need this for surface reflection for depth phases only?
				if ((criticalOnly) ? (legActionChange || topFluidLayer || botFluidLayer) : true ) RTproduct *= getIndividualRT(coeff, arrivalRayParam,interfaceR,waveIn,waveOut,RTaction);
			}
		}
		if (freeSurfaceSwitch) RTproduct *= getFreeSurfaceCorrection(arrivalRayParam,waveOut);
	}
	
	/*
	 * get reflection coefficient at specific interface
	 */
	private double getIndividualRT(RTcoeff coeff, double arrivalRayParam, double interfaceR, boolean waveIn, boolean waveOut, Integer action) {
		double RT=1;
		
		if (action>0 && action<3) { //reflection
			if ( waveIn && waveOut ) RT = coeff.getPtoPRefl(arrivalRayParam/interfaceR);
			if ( waveIn && !waveOut ) RT = coeff.getPtoSVRefl(arrivalRayParam/interfaceR);
			
			if ( !waveIn && !waveOut ) {
				if (!horizontal) RT = coeff.getSVtoSVRefl(arrivalRayParam/interfaceR);
				else RT = coeff.getSHtoSHRefl(arrivalRayParam/interfaceR);
			}
			if ( !waveIn && waveOut ) RT = coeff.getSVtoPRefl(arrivalRayParam/interfaceR);
		}
		else if (action>2) { // transmission
			if (waveIn && waveOut) RT = coeff.getPtoPTrans(arrivalRayParam/interfaceR);
			if (waveIn && !waveOut) RT = coeff.getPtoSVTrans(arrivalRayParam/interfaceR);
			if (!waveIn && !waveOut)  {
				if (!horizontal) RT = coeff.getSVtoSVTrans(arrivalRayParam/interfaceR);
				else RT = coeff.getSHtoSHTrans(arrivalRayParam/interfaceR);
			}
			if (!waveIn && waveOut) RT = coeff.getSVtoPTrans(arrivalRayParam/interfaceR);
		}
		
		return RT;
	}
	
	/*
	 * compute free surface correction
	 */
	private double getFreeSurfaceCorrection(double arrivalRayParam, boolean isPWave) 
    		throws NoSuchLayerException, NoSuchMatPropException {
    	double FSC=0.; //only returned if asked for vertical FSC for SH arrival
    	// need velocity model
    	VelocityModel vMod = thisPhase.getTauModel().getVelocityModel();
    	double planetRadius = vMod.getRadiusOfEarth();
		double vp = vMod.evaluateBelow(0., 'P');
		double vs = vMod.evaluateBelow(0., 'S');
		double density=vMod.evaluateBelow(0., 'r');
    	
    	RTcoeff coeff = new RTcoeff(vp,vs,density,0.,0.,0.); 
    	
    	if (isPWave) {
    		if (verticalComp) FSC=coeff.getPtoZFSC(arrivalRayParam/planetRadius);
    		else FSC=coeff.getPtoXFSC(arrivalRayParam/planetRadius);
    	}
    	else {
    		if (horizontal && !verticalComp) FSC=2;
    		else { // for SV wave
    			if (verticalComp) FSC=coeff.getSVtoZFSC(arrivalRayParam/planetRadius);
        		else FSC=coeff.getSVtoXFSC(arrivalRayParam/planetRadius);
    		}
    	}
    	return FSC;
    }
	
	/*
	 * compute geometrical spreading and amplitude factor acc. Aki and Richards
	 */
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
	    	
	    	// source side
	    	takeoffVelocity = (thisPhase.getDownGoing()[0] ? 1 : -1 ) *
					vMod.evaluateBelow(sourceDepth, name.charAt(0));
            takeoffAngle = thisPhase.calcTakeoffAngle(interpArrival.getRayParam())/180.*Math.PI;
	        
            // receiver side
            lastLeg = thisPhase.getLegs().get(thisPhase.getLegs().size()-2).charAt(0);
	        incidentVelocity = vMod.evaluateBelow(receiverDepth, lastLeg);
	        incidentAngle = thisPhase.calcIncidentAngle(interpArrival.getRayParam())/180*Math.PI;
	        
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
	                amplitudeFactor /= geomFact; //amplfact unit kg m/s^2
	        	}
	        }
	    }
	}
	
	public void compute(Arrival interpArrival) throws NoSuchLayerException, NoSuchMatPropException {
		computeAmplitudeFactor(interpArrival);
		computeRTfactor(interpArrival.getRayParam());
	}
	
	public double getAmplFact() {
		return amplitudeFactor;
	}
	
	public double getRT() {
		return RTproduct;
	}
	
	public String getPol() {
		return ((horizontal) ? "H" : "V");
	}
	
	public String getComp() {
		return ((verticalComp) ? "Z" : "X");
	}
}
