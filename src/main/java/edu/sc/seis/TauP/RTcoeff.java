/*
 * <pre> The TauP Toolkit: Flexible Seismic Travel-Time and Raypath Utilities.
 * Copyright (C) 1998-2000 University of South Carolina This program is free
 * software; you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version. This program
 * is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU General Public License for more details. You
 * should have received a copy of the GNU General Public License along with this
 * program; if not, write to the Free Software Foundation, Inc., 59 Temple Place -
 * Suite 330, Boston, MA 02111-1307, USA. The current version can be found at <A
 * HREF="www.seis.sc.edu">http://www.seis.sc.edu </A> Bug reports and comments
 * should be directed to H. Philip Crotwell, crotwell@seis.sc.edu or Tom Owens,
 * owens@seis.sc.edu </pre>
 */
/**
 * RTcoeff.java Reflection and transmission coefficients for body
 * waves. Methods for calculating coefficients for each of the possible
 * interactions are provided. Calculations are done using the
 * com.visualnumerics.javagrande.Complex class from VisualNumerics with two added methods. It is
 * further assume that the incoming ray is coming from the "top" for solid-solid
 * interactions and from the bottom for free surface interactions. If the ray is
 * actually coming from the bottom, flip the velocities. The convention for
 * free surface and solid solid is a little strange, but can be thought of as
 * the top velocities correspond to the layer that they ray starts in.
 * 
 * @see "Cerveny 2001 page 480+"
 * @author Stefanie Hempel
 * @version 0 Fri Feb 19 11:00:00 GMT 2016
 */
package edu.sc.seis.TauP;

import java.io.Serializable;

public class RTcoeff implements Serializable {

    // IMPORTANT!!!!
    // Where ever "CX" appears in this class, it is used as a shorthand for
    // the Complex class, so CX.times() is the same as Complex.times, but
    // the code is, IMHO, less cluttered (Philip Crotwell).
    /** just to avoid having Complex all over the place. */
    private static final Complex CX = new Complex();
    
    private static boolean matlabTests = false;

    protected double topVp;

    protected double topVs;

    protected double topDensity;

    protected double botVp;

    protected double botVs;

    protected double botDensity;

    // "flat earth" ray parameter
    protected double rp;

    // temp variables to make calculations less ugly
    protected double q, X, Y, Z;

    protected Complex P1, P2, P3, P4;

    /** used only in free surface calculations */
    protected Complex DS;
    protected Complex D;

    /**
     * delta for SH-SH equations
     */
    protected Complex D1;

    // we need the squared terms so often that it is worthwhile to store them
    protected double sqBotVs; // botVs squared

    protected double sqTopVs; // topVs squared

    protected double sqBotVp; // botVp squared

    protected double sqTopVp; // topVp squared

    protected double sqRP; // rp squared

    // remember last calculated ray param and wave type to avoid repeating
    protected double lastRayParam = -1.0;

    protected boolean lastInIsPWave = true;

    // CM
    protected boolean firstTime = true;

    // CM
    public RTcoeff(double topVp,
                                double topVs,
                                double topDensity,
                                double botVp,
                                double botVs,
                                double botDensity) {
        this.topVp = topVp;
        this.topVs = topVs;
        this.topDensity = topDensity;
        this.botVp = botVp;
        this.botVs = botVs;
        this.botDensity = botDensity;
    }

    /**
     * Flips the sense of the layers, useful when you have a ray going through
     * the same layer in the opposite direction.
     */
    public RTcoeff flip() {
        return new RTcoeff(botVp,
                                        botVs,
                                        botDensity,
                                        topVp,
                                        topVs,
                                        topDensity);
    }

    protected void calcTempVars(double rayParam, boolean inIsPWave) {
        if(rayParam < 0) {
            throw new IllegalArgumentException("rayParam cannot be negative");
        }
        this.rp = rayParam; // ray parameter
        Complex tmp;
        
        // CM
        // if (rayParam != lastRayParam && inIsPWave == lastInIsPWave ) {
        // if ( (rayParam != lastRayParam || inIsPWave != lastInIsPWave ) ||
        // firstTime ) {
        if(rayParam != lastRayParam || inIsPWave != lastInIsPWave) {
            lastRayParam = -1.0; // in case of failure in method
            // CM
            firstTime = false;
            sqBotVs = botVs * botVs; // botVs squared
            sqTopVs = topVs * topVs; // topVs squared
            sqBotVp = botVp * botVp; // botVp squared
            sqTopVp = topVp * topVp; // topVp squared
            sqRP = rp * rp; // rp squared
            
            q = 2 * (botDensity * sqBotVs - topDensity * sqTopVs);
            X = botDensity - q * sqRP;
            Y = topDensity + q * sqRP;
            Z = botDensity - topDensity - q * sqRP;
            
            // these should be positive, even if imaginary
            P1 = CX.sqrt(new Complex ( 1 - sqTopVp * sqRP)).posImag();
            P2 = CX.sqrt(new Complex ( 1 - sqTopVs * sqRP)).posImag();
            P3 = CX.sqrt(new Complex ( 1 - sqBotVp * sqRP)).posImag();
            P4 = CX.sqrt(new Complex ( 1 - sqBotVs * sqRP)).posImag();
            
            // first summand to D
            Complex ATerm = CX.times( CX.times( P1, P2 ), CX.times( P3, P4 ) ).times( q*q*sqRP  );
            // second summand to D
            Complex BTerm = CX.plus( CX.times( P1, P4 ).times(topVs * botVp),
            		CX.times( P2, P3 ).times(topVp * botVs)).times(topDensity * botDensity);
            // third summand to D
            Complex CTerm = CX.times( P3, P4 ).times(topVp * topVs ).times( Math.pow(Y, 2) );
            // forth summand
            Complex DTerm = CX.times( P1, P2 ).times(botVp * botVs ).times( Math.pow(X, 2) );
            // fifth summand
            D = CX.plus(CX.plus(ATerm,BTerm),CX.plus(CTerm,DTerm)).plus(topVp * botVp * topVs * botVs * sqRP * Math.pow(Z, 2));
            
            // D1
            D1 = CX.plus( CX.times(topDensity * topVs, P2),
            		CX.times(botDensity * botVs, P4) );
            
            //D for Free Surface Correction
            DS = CX.plus(Math.pow(1 - 2*Math.pow(topVs,2)*sqRP,2),
            		CX.times(4*sqRP*Math.pow(topVs,3)/topVp,CX.times(P1,P2)));

            lastRayParam = rayParam;
            lastInIsPWave = inIsPWave;
        }
    }

    // extra handling of free surface not necessary in this formulation

    // Solid-Solid interface
    /**
     * Calculates incident P wave to reflected P wave Complex coefficient.
     * [q² p² P1 P2 P3 P4 - ρ1 ρ2 (α1 β2 P2 P3 − β1 α2 P1 P4 )
	 *	− α1 β1 P3 P4 Y² + α2 β2 P1 P2 X² − α1 α2 β1 β2 p² Z²] / D
     */
    public Complex getComplexPtoPRefl(double rayParam) {
        calcTempVars(rayParam, true);
//        Complex ATerm = CX.times(CX.times(P1, P2),
//                                          CX.times(P3, P4)).
//        		times(new Complex ( q*q * sqRP ) );
//        Complex BTerm = CX.plus(CX.times(P2, P3).times(new Complex ( topVp * botVs )),
//        		CX.times(P1, P4).times(new Complex ( topVs * botVp )).times(new Complex ( topDensity * botDensity ) ) );
//        Complex CTerm = CX.times(P3, P4).times(new Complex ( topVp * topVs * Y * Y));
//        Complex DTerm = CX.times(P1, P2).times(new Complex ( botVp * botVs * X * X));
//        Complex ETerm = new Complex ( topVp * botVp * topVs *botVs * sqRP * Z * Z);
//
//        Complex numerator = CX.plus(ATerm,BTerm).plus(DTerm).plus(CTerm).plus(ETerm);
        
        // first summand
        Complex ATerm = CX.times( CX.times( P1, P2), CX.times( P3, P4)).times( q*q*sqRP  );
        // second summand
        Complex BTerm = CX.minus( CX.times( P1, P4 ).times(topVs * botVp),
        		CX.times( P2, P3 ).times(topVp * botVs)).times(topDensity * botDensity);
        // third summand
        Complex CTerm = CX.times( P3, P4 ).times(topVp * topVs ).times( Math.pow(Y, 2) );
        // forth summand
        Complex DTerm = CX.times( P1, P2 ).times(botVp * botVs ).times( Math.pow(X, 2) );
        // combine w fifth summand
        Complex numerator = CX.plus(CX.plus(ATerm,BTerm),CX.minus(DTerm,CTerm)).minus(topVp * botVp * topVs * botVs * sqRP * Math.pow(Z, 2));
        
        return CX.over(numerator, D);
    }

    /**
     * Calculates incident P wave to reflected P wave coefficient.
     */
    public double getPtoPRefl(double rayParam) {
    	Complex val=getComplexPtoPRefl(rayParam);
        return Math.signum(CX.real(val))*CX.abs(val);
        //return CX.real(getComplexPtoPRefl(rayParam));
    }
    public double getPtoPReflPhase(double rayParam) {
        // return CX.imag(getComplexPtoPRefl(rayParam));
    	return CX.argument(getComplexPtoPRefl(rayParam));
    }

    /**
     * Calculates incident P wave to reflected SV wave Complex coefficient.
     * 2 e α1 p P1 (q P3 P4 Y + α2 β2 X Z ) / D
     */
    public Complex getComplexPtoSVRefl(double rayParam) {
        calcTempVars(rayParam, true);
        double epsilon = 1;
        Complex numerator = CX.plus( CX.times(P3, P4).times(new Complex ( q * Y ) ) , new Complex ( botVp * botVs * X * Z)). 
         times (CX.times(P1, new Complex (2 * epsilon * rp * topVp) ) );
        return CX.over(numerator, D);
    }

    /**
     * Calculates incident P wave to reflected SV wave coefficient.
     */
    public double getPtoSVRefl(double rayParam) {
    	return CX.abs(getComplexPtoSVRefl(rayParam));
        //return CX.real(getComplexPtoSVRefl(rayParam));
    }
    public double getPtoSVReflPhase(double rayParam) {
        // return CX.imag(getComplexPtoSVRefl(rayParam));
        return CX.argument(getComplexPtoSVRefl(rayParam));
    }

    /**
     * Calculates incident P wave to transmitted P wave Complex coefficient.
     * 2α1 ρ1 P1 (β2 P2 X + β1 P4 Y ) / D
     */
    public Complex getComplexPtoPTrans(double rayParam) {
        calcTempVars(rayParam, true);
        Complex numerator = CX.times( CX.times (2* topVp * topDensity, P1 ),
        		CX.plus( CX.times( botVs, P2 ).times ( X ),
        		CX.times(topVs, P4).times( Y ) ) );
        return CX.over(numerator, D);
    }

    /**
     * Calculates incident P wave to transmitted P wave coefficient.
     */
    public double getPtoPTrans(double rayParam) {
        // return CX.real(getComplexPtoPTrans(rayParam));
    	Complex val=getComplexPtoPTrans(rayParam);
        return Math.signum(CX.real(val))*CX.abs(val);
    }
    public double getPtoPTransPhase(double rayParam) {
        // return CX.imag(getComplexPtoPTrans(rayParam));
        return CX.argument(getComplexPtoPTrans(rayParam));
    }
    
    public Complex getComplexPtoXFreeSurfaceConversion(double rayParam) {
        calcTempVars(rayParam, true);
        Complex numerator =CX.times(4 * topVs * rayParam, CX.times(P1, P2 ));
        return CX.over(numerator, DS);
    }
    public double getPtoXFSC(double rayParam) {
        // return CX.real(getComplexPtoPTrans(rayParam));
    	Complex val=getComplexPtoXFreeSurfaceConversion(rayParam);
        return Math.signum(CX.real(val))*CX.abs(val);
    }
    public double getPtoXFSCPhase(double rayParam) {
        // return CX.imag(getComplexPtoPTrans(rayParam));
        return CX.argument(getComplexPtoXFreeSurfaceConversion(rayParam));
    }
    
    public Complex getComplexPtoZFreeSurfaceConversion(double rayParam) {
        calcTempVars(rayParam, true);
        double epsilon = 1;
        Complex numerator =CX.times(2 * epsilon, P1).times(1-2*Math.pow(topVs,2)*sqRP);
        return CX.over(numerator, DS);
    }
    public double getPtoZFSC(double rayParam) {
        // return CX.real(getComplexPtoPTrans(rayParam));
    	Complex val=getComplexPtoZFreeSurfaceConversion(rayParam);
        return Math.signum(CX.real(val))*CX.abs(val);
    }
    public double getPtoZFSCPhase(double rayParam) {
        // return CX.imag(getComplexPtoPTrans(rayParam));
        return CX.argument(getComplexPtoZFreeSurfaceConversion(rayParam));
    }

    /**
     * Calculates incident P wave to transmitted SV wave Complex coefficient.
     * −2 epsilon α1 ρ1 p P1 (q P2 P3 − β1 α2 Z ) / D
     */
    public Complex getComplexPtoSVTrans(double rayParam) {
        calcTempVars( rayParam , true );
        double epsilon = 1;
        Complex numerator = CX.times( CX.times(-2 * epsilon * topVp * topDensity * rp, P1 ), 
        		CX.minus( CX.times( q , P2 ).times( P3 ),
        		 topVs * botVp * Z ));
        return CX.over( numerator , D );
    }

    /**
     * Calculates incident P wave to transmitted SV wave coefficient.
     */
    public double getPtoSVTrans(double rayParam) {
        // return CX.real(getComplexPtoSVTrans(rayParam));
    	Complex val=getComplexPtoSVTrans(rayParam);
        return Math.signum(CX.real(val))*CX.abs(val);
    }
    public double getPtoSVTransPhase(double rayParam) {
    	// return CX.imag(getComplexPtoSVTrans(rayParam));
        return CX.argument( getComplexPtoSVTrans( rayParam ) );
    }
    
    /**
     * Calculates incident SV wave to reflected P wave Complex coefficient.
     * −2 epsilon β1 p P2 (q P3 P4 Y + α2 β2 X Z ) / D
     */
    public Complex getComplexSVtoPRefl(double rayParam) {
        calcTempVars(rayParam, false);
        double epsilon = 1;
        Complex numerator = CX.plus( CX.times(P3, P4).times(new Complex ( q * Y ) ) , new Complex ( botVp * botVs * X * Z)). 
         times (CX.times(P2, new Complex (- 2 * epsilon * rp * topVs) ) );
        return CX.over(numerator, D);
    }

    /**
     * Calculates incident SV wave to reflected P wave coefficient.
     */
    public double getSVtoPRefl(double rayParam) {
        // return CX.real(getComplexSVtoPRefl(rayParam));
    	Complex val=getComplexSVtoPRefl(rayParam);
        return Math.signum(CX.real(val))*CX.abs(val);
    }
    public double getSVtoPReflPhase(double rayParam) {
        // return CX.imag(getComplexSVtoPRefl(rayParam));
        return CX.argument(getComplexSVtoPRefl(rayParam));
    }

    /**
     * Calculates incident SV wave to reflected SV wave Complex coefficient.
     * <P>= [q² p² P1 P2 P3 P4 + ρ1 ρ2 (α1 β2 P2 P3 − β1 α2 P1 P4 )
     *	− α1 β1 P3 P4 Y² + α2 β2 P1 P2 X² − α1 α2 β1 β2 p² Z²] / D
     */
    public Complex getComplexSVtoSVRefl(double rayParam) {
    	calcTempVars(rayParam, true);
    	
    	// first summand
        Complex ATerm = CX.times( CX.times( P1, P2), CX.times( P3, P4)).timesReverse( q*q*sqRP  );
        // second summand
        Complex BTerm = CX.plus( CX.times( P1, P4 ).timesReverse(-1 * topVs * botVp),
        		CX.times( P2, P3 ).timesReverse(topVp * botVs)).timesReverse(topDensity * botDensity);
//    	Complex BTerm = CX.minus( CX.times( P1, P4 ).times(topVs * botVp),
//        		CX.times( P2, P3 ).times(topVp * botVs)).times(topDensity * botDensity);
        // third summand
        Complex CTerm = CX.times( P3, P4 ).timesReverse(-1 * topVp * topVs ).times( Math.pow(Y, 2) );
        // forth summand
        Complex DTerm = CX.times( P1, P2 ).timesReverse(botVp * botVs ).times( Math.pow(X, 2) );
        // combine w fifth summand
        Complex numerator = CX.plus(CX.plus(ATerm,BTerm),CX.plus(DTerm,CTerm)).plus(-1 * topVp * botVp * topVs * botVs * sqRP * Math.pow(Z, 2));
            
    	return CX.over(numerator, D);
	}

    /**
     * Calculates incident SV wave to reflected SV wave coefficient.
     */
    public double getSVtoSVRefl(double rayParam) {
        // return CX.real(getComplexSVtoSVRefl(rayParam));
    	Complex val=getComplexSVtoSVRefl(rayParam);
        return Math.signum(CX.real(val))*CX.abs(val);
    }
    public double getSVtoSVReflPhase(double rayParam) {
        // return CX.imag(getComplexSVtoSVRefl(rayParam));
        return CX.argument(getComplexSVtoSVRefl(rayParam));
    }
    
    /**
     * Calculates incident SV wave to transmitted P wave Complex coefficient.
     * 2 epsilon β1 ρ1 p P2 (q P1 P4 − α1 β2 Z ) / D
     */
    public Complex getComplexSVtoPTrans(double rayParam) {
        calcTempVars(rayParam, false);
        double epsilon = 1;
        Complex numerator = CX.minus( CX.times( P1 , P4 ).times( q ),
        		new Complex( topVp * botVs * Z ) ).
        		times( CX.times (P2, new Complex( 2 * epsilon * topVs * topDensity * rp ) ) );
        return CX.over(numerator, D);
    }

    /**
     * Calculates incident SV wave to transmitted P wave coefficient.
     */
    public double getSVtoPTrans(double rayParam) {
        // return CX.real(getComplexSVtoPTrans(rayParam));
    	Complex val=getComplexSVtoPTrans(rayParam);
        return Math.signum(CX.real(val))*CX.abs(val);
    }
    public double getSVtoPTransPhase(double rayParam) {
        // return CX.imag(getComplexSVtoPTrans(rayParam));
        return CX.argument(getComplexSVtoPTrans(rayParam));
    }

    /**
     * Calculates incident SV wave to transmitted SV wave Complex coefficient.
     * 2β1 ρ1 P2 (α1 P3 Y + α2 P1 X ) / D
     */
    public Complex getComplexSVtoSVTrans(double rayParam) {
        calcTempVars(rayParam, false);
        Complex numerator = CX.plus( CX.times(P3, new Complex ( topVp * Y) ) ,
        		CX.times(P1, new Complex( botVp * X ) ) ).
        		times( CX.times (P2, new Complex( 2 * topVs * topDensity ) ) );
        return CX.over(numerator, D);
    }

    /**
     * Calculates incident SV wave to transmitted SV wave coefficient.
     */
    public double getSVtoSVTrans(double rayParam) {
        // return CX.real(getComplexSVtoSVTrans(rayParam));
    	Complex val=getComplexSVtoSVTrans(rayParam);
        return Math.signum(CX.real(val))*CX.abs(val);
    }
    public double getSVtoSVTransPhase(double rayParam) {
        // return CX.imag(getComplexSVtoSVTrans(rayParam));
        return CX.argument(getComplexSVtoSVTrans(rayParam));
    }
    
    //Free Surface Correction 
    public Complex getComplexSVtoXFreeSurfaceConversion(double rayParam) {
        calcTempVars(rayParam, false);
        double epsilon = 1;
        Complex numerator =CX.times(2 * epsilon, P2).times(1-2*Math.pow(topVs,2)*sqRP);
        return CX.over(numerator, DS);
    }
    public double getSVtoXFSC(double rayParam) {
        // return CX.real(getComplexPtoPTrans(rayParam));
    	Complex val=getComplexSVtoXFreeSurfaceConversion(rayParam);
        return Math.signum(CX.real(val))*CX.abs(val);
    }
    public double getSVtoXFSCPhase(double rayParam) {
        // return CX.imag(getComplexPtoPTrans(rayParam));
        return CX.argument(getComplexSVtoXFreeSurfaceConversion(rayParam));
    }
    
    public Complex getComplexSVtoZFreeSurfaceConversion(double rayParam) {
        calcTempVars(rayParam, false);
        Complex numerator =CX.times(-4 * Math.pow(topVs,2) * rayParam, CX.times(P1, P2 )).over(topVp);
        return CX.over(numerator, DS);
    }
    public double getSVtoZFSC(double rayParam) {
        // return CX.real(getComplexPtoPTrans(rayParam));
    	Complex val=getComplexSVtoZFreeSurfaceConversion(rayParam);
        return Math.signum(CX.real(val))*CX.abs(val);
    }
    public double getSVtoZFSCPhase(double rayParam) {
        // return CX.imag(getComplexPtoPTrans(rayParam));
        return CX.argument(getComplexSVtoZFreeSurfaceConversion(rayParam));
    }

    // SH waves
    /**
     * Calculates incident SH wave to reflected SH wave Complex coefficient.
     *  ̄(ρ1 β1 P2 − ρ2 β2 P4 ) / D1
     */
    public Complex getComplexSHtoSHRefl(double rayParam) {
        calcTempVars(rayParam, false);
        Complex numerator = CX.minus( CX.times( new Complex ( topDensity * topVs ) , P2 ),
        		CX.times( new Complex ( botDensity * botVs ) , P4 ) );
        return CX.over(numerator, D1);
    }

    /**
     * Calculates incident SH wave to reflected SH wave coefficient.
     */
    public double getSHtoSHRefl(double rayParam) {
        // return CX.real(getComplexSHtoSHRefl(rayParam));
    	Complex val=getComplexSHtoSHRefl(rayParam);
        return Math.signum(CX.real(val))*CX.abs(val);
        }
    
    public double getSHtoSHReflPhase(double rayParam) {
        // return CX.imag(getComplexSHtoSHRefl(rayParam));
        return CX.argument(getComplexSHtoSHRefl(rayParam));
    }

    /**
     * Calculates incident SH wave to transmitted SH wave Complex coefficient.
     *  ̄2ρ1 β1 P2 / D1
     */
    public Complex getComplexSHtoSHTrans(double rayParam) {
        calcTempVars(rayParam, false);
        Complex numerator = CX.times( P2, new Complex ( 2 * topVs * topDensity ) );
        return CX.over(numerator, D1);
    }

    /**
     * Calculates incident SH wave to transmitted SH wave coefficient.
     */
    public double getSHtoSHTrans(double rayParam) {
        // return CX.real(getComplexSHtoSHTrans(rayParam));
    	Complex val=getComplexSHtoSHTrans(rayParam);
        return Math.signum(CX.real(val))*CX.abs(val);
    }
    public double getSHtoSHTransPhase(double rayParam) {
        // return CX.imag(getComplexSHtoSHTrans(rayParam));
        return CX.argument(getComplexSHtoSHTrans(rayParam));
    }

    public static void main(String[] args) {
    	// check free surface correction
    	double vp1=6.400;
    	double vs1=3.698;
    	double rho1=2.980;
    	
    	// check vs. Cerveny's plots
        double topVp = 6.4;
        double topVs = 3.698;
        double topDensity = 2.98;
        double botVp = 8.0;
        double botVs = 4.618;
        double botDensity = 3.3;
    	// Mars CMB
//	    double topVp = 9.569;
//	    double topVs = 5.053;
//	    double topDensity = 4.065;
//	    double botVp = 4.825;
//	    double botVs = 0.000;
//	    double botDensity = 5.937;
	    // Earth ICB
//	    double topVp = 10.35568;
//	    double topVs = 0.000;
//	    double topDensity = 12.16634;
//	    double botVp = 11.02827;
//	    double botVs = 3.50432;
//	    double botDensity = 12.76360;
	    // Earth CMB
//	    double topVp = 13.71660;
//	    double topVs = 7.26466;
//	    double topDensity = 5.56645;
//	    double botVp = 8.06482;
//	    double botVs = 0.00000;
//	    double botDensity = 9.90349;
	    
        double depth;
        double radiusOfEarth;
        double DtoR = Math.PI / 180.0;
        double RtoD = 180.0 / Math.PI;
        double[] RPP = new double[91];
        double[] RPS = new double[91];
        double[] RSP = new double[91];
        double[] RSS = new double[91];
        double[] TPP = new double[91];
        double[] TPS = new double[91];
        double[] TSP = new double[91];
        double[] TSS = new double[91];
        double[] RPPphase = new double[91];
        double[] RPSphase = new double[91];
        double[] RSPphase = new double[91];
        double[] RSSphase = new double[91];
        double[] TPPphase = new double[91];
        double[] TPSphase = new double[91];
        double[] TSPphase = new double[91];
        double[] TSSphase = new double[91];
        double[] RSH = new double[91];
        double[] TSH = new double[91];
        double[] RSHphase = new double[91];
        double[] TSHphase = new double[91];
        RTcoeff coeff = new RTcoeff(topVp,
                                                              topVs,
                                                              topDensity,
                                                              botVp,
                                                              botVs,
                                                              botDensity);
        // Free Surface Correction
        double[] FPX = new double[91];
        double[] FPZ = new double[91];
        double[] FPXphase = new double[91];
        double[] FPZphase = new double[91];
        double[] FSVX = new double[91];
        double[] FSVZ = new double[91];
        double[] FSVXphase = new double[91];
        double[] FSVZphase = new double[91];
        RTcoeff fcoeff = new RTcoeff(vp1,vs1,rho1,0.,0.,0.);
        
        //coeff = coeff.flip();
        double rayParam;
        System.out.println("current test version 1.0");
        for(int i = 0; i <= 90; i++) {
        	try {
        	rayParam = Math.sin(DtoR * i) / topVp;
            RPP[i] = coeff.getPtoPRefl(rayParam);
            RPS[i] = coeff.getPtoSVRefl(rayParam);
            RPPphase[i] = coeff.getPtoPReflPhase(rayParam)/Math.PI;
            RPSphase[i] = coeff.getPtoSVReflPhase(rayParam)/Math.PI;
            TPP[i] = coeff.getPtoPTrans(rayParam);
            TPS[i] = coeff.getPtoSVTrans(rayParam);
            TPPphase[i] = coeff.getPtoPTransPhase(rayParam)/Math.PI;
            TPSphase[i] = coeff.getPtoSVTransPhase(rayParam)/Math.PI;
            
            // Free Surface Correction
            FPX[i] = fcoeff.getPtoXFSC(rayParam);
            FPXphase[i] = fcoeff.getPtoXFSCPhase(rayParam);
            FPZ[i] = fcoeff.getPtoZFSC(rayParam)/Math.PI;
            FPZphase[i] = fcoeff.getPtoZFSCPhase(rayParam)/Math.PI;
            
            
            
            rayParam = Math.sin(DtoR * i) / topVs;
            RSP[i] = coeff.getSVtoPRefl(rayParam);
            RSS[i] = coeff.getSVtoSVRefl(rayParam);
            RSPphase[i] = coeff.getSVtoPReflPhase(rayParam)/Math.PI;
            RSSphase[i] = coeff.getSVtoSVReflPhase(rayParam)/Math.PI;
            RSH[i] = coeff.getSHtoSHRefl(rayParam);
            RSHphase[i] = coeff.getSHtoSHReflPhase(rayParam)/Math.PI;
            
            TSP[i] = coeff.getSVtoPTrans(rayParam);
            TSS[i] = coeff.getSVtoSVTrans(rayParam);
            TSPphase[i] = coeff.getSVtoPTransPhase(rayParam)/Math.PI;
            TSSphase[i] = coeff.getSVtoSVTransPhase(rayParam)/Math.PI;
            TSH[i] = coeff.getSHtoSHTrans(rayParam);
            TSHphase[i] = coeff.getSHtoSHTransPhase(rayParam)/Math.PI;
            
            // Free Surface Correction
            FSVX[i] = fcoeff.getSVtoXFSC(rayParam);
            FSVXphase[i] = fcoeff.getSVtoXFSCPhase(rayParam);
            FSVZ[i] = fcoeff.getSVtoZFSC(rayParam)/Math.PI;
            FSVZphase[i] = fcoeff.getSVtoZFSCPhase(rayParam)/Math.PI;
        	}
        	catch(RuntimeException ex) {
        		throw new RuntimeException( ex.getMessage() );
        	}
        }
        try {
        	if (matlabTests) {
        		// Cerveny reference
        		java.io.Writer out = new java.io.BufferedWriter(new java.io.FileWriter("freesurfacecorrection.dat"));
        		for (int i=0; i<=90; i++) {	 
                	 //rayParam = Math.sin(DtoR * i) / topVp;
                	 //out.write(i+" "+RPP[i]+"\n");
                	 out.write(i +" "+FPX[i]+" "+FPXphase[i]+" "+FPZ[i]+" "+FPZphase[i]+" "
                			 		 +FSVX[i]+" "+FSVXphase[i]+" "+FSVZ[i]+" "+FSVZphase[i]+"\n");
                }
                out.close();
                
//        		java.io.Writer out = new java.io.BufferedWriter(new java.io.FileWriter("refltrans2.dat"));
//        		for (int i=0; i<=90; i++) {	 
//                	 //rayParam = Math.sin(DtoR * i) / topVp;
//                	 //out.write(i+" "+RPP[i]+"\n");
//                	 out.write(i +" "+RPP[i]+" "+RPPphase[i]+" "+RPS[i]+" "+RPSphase[i]+" "
//                			 		 +TPP[i]+" "+TPPphase[i]+" "+TPS[i]+" "+TPSphase[i]+" "
//                			 		 +RSS[i]+" "+RSSphase[i]+" "+RSP[i]+" "+RSPphase[i]+" "
//                			 		 +TSS[i]+" "+TSSphase[i]+" "+TSP[i]+" "+TSPphase[i]+" "
//                			 		 +RSH[i]+" "+RSHphase[i]+" "+TSH[i]+" "+TSHphase[i]+"\n");
//                }
//                out.close();
                
//                //for Mars CMB tests
//                java.io.Writer out = new java.io.BufferedWriter(new java.io.FileWriter("refltrans2PICB.dat"));
//        		for (int i=0; i<=90; i++) {	 
//                	 rayParam = Math.sin(DtoR * i) / topVp;
//                	 //out.write(i+" "+RPP[i]+"\n");
//                	 out.write(i +" "+rayParam+" "+RPP[i]+" "+RPPphase[i]+" "
//                			 		 +TPP[i]+" "+TPPphase[i]+" "
//                			 		 +TPS[i]+" "+TPSphase[i]+"\n");
//                }
//                out.close();
//                out = new java.io.BufferedWriter(new java.io.FileWriter("refltrans2SICB.dat"));
//        		for (int i=0; i<=90; i++) {	 
//                	 rayParam = Math.sin(DtoR * i) / topVs;
//                	 //out.write(i+" "+RPP[i]+"\n");
//                	 out.write(i +" "+rayParam+" "+RSS[i]+" "+RSSphase[i]+" "
//                			 		 +TSS[i]+" "+TSSphase[i]+" "
//        			 				 +TSP[i]+" "+TSPphase[i]+" "
//                			 		 +RSH[i]+" "+RSHphase[i]+" "+TSH[i]+" "+TSHphase[i]+"\n");
//                }
//                out.close();
        	}
        	else {
	            java.io.Writer out = new java.io.BufferedWriter(new java.io.FileWriter("refltrans2.gmt"));
	            out.write("#!/bin/sh\n\n");
	            out.write("/bin/rm -f refltrans2.ps\n\n");
	            //out.write("psbasemap -K -P -JX12 -R0/90/-1/2  -B10/1 > refltrans2.ps\n");
	            out.write("gmt psbasemap -K -P -JX12 -R0/95/-1/2  -Bp10/.5,s1/.1 > refltrans2.ps\n");
	             out.write("gmt psxy -K -O -JX -R -W2,black >> refltrans2.ps <<END\n");
	             for (int i=0; i<=90; i++) {	 
	            	 rayParam = Math.sin(DtoR * i) / topVp;
	            	 //out.write(i+" "+RPP[i]+"\n");
	            	 out.write(i +" "+RPP[i]+"\n");
	             }
	             out.write("END\n");
	             
	             out.write("gmt psxy -K -O -JX -R -W2,blue >> refltrans2.ps <<END\n");
	             for (int i=0; i<=90; i++) {
	            	 rayParam = Math.sin(DtoR * i) / topVp;
	            	 //out.write(i+" "+RPS[i]+"\n");
	            	 out.write(i+" "+RPS[i]+"\n");
	             }
	             out.write("END\n");
	             
	             out.write("gmt psxy -K -O -JX -R -W2,black,. >> refltrans2.ps <<END\n");
	             for (int i=0; i<=90; i++) {
	            	 rayParam = Math.sin(DtoR * i) / topVp;
	            	 //out.write(i+" "+TPP[i]+"\n");
	            	 out.write(i+" "+TPP[i]+"\n");
	             }
	             out.write("END\n");
	             
	             out.write("gmt psxy -K -O -JX -R -W2,blue,. >> refltrans2.ps <<END\n");
	             for (int i=0; i<=90; i++) {
	            	 rayParam = Math.sin(DtoR * i) / topVp;
	            	 //out.write(i+" "+TPS[i]+"\n");
	            	 out.write(i+" "+TPS[i]+"\n");
	             }
	             out.write("END\n");
	             
	             out.write("gmt psxy -K -O -JX -R -W2,red >> refltrans2.ps <<END\n");
	             out.write("gmt psxy -K -O -JX -R -W2,red >> refltrans2.ps <<END\n");
	             out.write("gmt psxy -K -O -JX -R -W2,red >> refltrans2.ps <<END\n");
	             out.write("gmt psxy -K -O -JX -R -W2,red >> refltrans2.ps <<END\n");
	             out.write("gmt psxy -K -O -JX -R -W2,red >> refltrans2.ps <<END\n");
	             out.write("gmt psxy -K -O -JX -R -W2,red >> refltrans2.ps <<END\n");
	             for(int i = 0; i <= 90; i++) {
	            	 rayParam = Math.sin(DtoR * i) / topVs;
	            	 //out.write(i + " " + RSP[i] + "\n");
	            	 out.write(i + " " + RSP[i] + "\n");
	          	 }
	            out.write("END\n");
	            out.write("gmt psxy -K -O -JX -R -W2,magenta >> refltrans2.ps <<END\n");
	            for(int i = 0; i <= 90; i++) {
	            	rayParam = Math.sin(DtoR * i) / topVs;
	                //out.write(i + " " + RSS[i] + "\n");
	            	out.write(i + " " + RSS[i] + "\n");
	            }
	            out.write("END\n");
	            out.write("END\n");
	            out.write("END\n");
	            out.write("END\n");
	            out.write("END\n");
	            out.write("END\n");
	            out.write("END\n");
	            out.write("END\n");
	            out.write("gmt psxy -K -O -JX -R -W2,red,. >> refltrans2.ps <<END\n");
	            for(int i = 0; i <= 90; i++) {
	            	rayParam = Math.sin(DtoR * i) / topVs;
	                //out.write(i + " " + TSP[i] + "\n");
	            	out.write(i + " " + TSP[i] + "\n");
	            }
	            out.write("END\n");
	            out.write("gmt p:qsxy -O -JX -R -W2,magenta,. >> refltrans2.ps <<END\n");
	            for(int i = 0; i <= 90; i++) {
	            	rayParam = Math.sin(DtoR * i) / topVs;
	                //out.write(i + " " + TSS[i] + "\n");
	                out.write(i + " " + TSS[i] + "\n");
	            }
	            out.write("END\n");
	            
	            out.close();
        	}
        } catch(java.io.IOException e) {
            System.err.println(e);
        }
    }
    
    public double getMinVs() {
    	double sign = Math.signum(botVs-topVs); //negative if bot is lower
    	double minimum = Math.min(botVs, topVs);
    	if (minimum==0) { return sign; }
    	else {
    		return 1;
    	}
    }
    
    public double getTopVs() {
    	return topVs;
    }
    
    public double getBotVs() {
    	return botVs;
    }
    
    public RTcoeff autoFlip() {
    	if (topVs == 0) {
    		return this.flip();
    	}
    	else {
    		return this;
    	}
    }
} // RTcoeff