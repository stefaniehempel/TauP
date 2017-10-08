/*
 * The TauP Toolkit: Flexible Seismic Travel-Time and Raypath Utilities.
 * Copyright (C) 1998-2000 University of South Carolina
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place - Suite 330, Boston, MA 02111-1307, USA.
 * 
 * The current version can be found at <A
 * HREF="www.seis.sc.edu">http://www.seis.sc.edu</A>
 * 
 * Bug reports and comments should be directed to H. Philip Crotwell,
 * crotwell@seis.sc.edu or Tom Owens, owens@seis.sc.edu
 * 
 */
package edu.sc.seis.TauP;

import java.time.Duration;

/**
 * convenience class for storing the parameters associated with a phase arrival.
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 * @author H. Philip Crotwell
 * 
 * Modified to account for dddp,  amplitude computation.
 * S. Hempel/ ISAE Toulouse Sep 2017
 * 
 */
public class Arrival {


    public Arrival(SeismicPhase phase,
                   double time,
                   double dist,
                   double rayParam,
                   int rayParamIndex,
                   String name,
                   String puristName,
                   double sourceDepth) {
        this(phase,
             time,
             dist,
             rayParam,
             rayParamIndex,
             name,
             puristName,
             sourceDepth,
             phase.calcTakeoffAngle(rayParam),
             phase.calcIncidentAngle(rayParam));
    }

    public Arrival(SeismicPhase phase,
                   double time,
                   double dist,
                   double rayParam,
                   int rayParamIndex,
                   String name,
                   String puristName,
                   double sourceDepth,
                   double takeoffAngle,
                   double incidentAngle) {
        if (Double.isNaN(time)) {
            throw new IllegalArgumentException("Time cannot be NaN");
        }
        if (rayParamIndex < 0) {
            throw new IllegalArgumentException("rayParamIndex cannot be negative: "+rayParamIndex);
        }
        this.phase = phase;
        this.time = time;
        this.dist = dist;
        this.rayParam = rayParam;
        this.rayParamIndex = rayParamIndex;
        this.name = name;
        this.puristName = puristName;
        this.sourceDepth = sourceDepth;
        this.takeoffAngle = takeoffAngle;
        this.incidentAngle = incidentAngle;
    }
    
    public Arrival(SeismicPhase phase,
            double time,
            double dist,
            double rayParam,
            int rayParamIndex,
            String name,
            String puristName,
            double sourceDepth,
            double takeoffAngle,
            double incidentAngle,
            double dddp,
            double amplFact,
            double RTFact,
            double tstar,
            double dtel,
            double dtdh,
            double corr,
            double corrError) {
    	if (Double.isNaN(time)) {
            throw new IllegalArgumentException("Time cannot be NaN");
        }
        if (rayParamIndex < 0) {
            throw new IllegalArgumentException("rayParamIndex cannot be negative: "+rayParamIndex);
        }
		this.phase = phase;
		this.time = time;
		this.dist = dist;
		this.rayParam = rayParam;
		this.rayParamIndex = rayParamIndex;
		this.name = name;
		this.puristName = puristName;
		this.sourceDepth = sourceDepth;
		this.takeoffAngle = takeoffAngle;
		this.incidentAngle = incidentAngle;
		this.amplFact = amplFact;
		this.RTFact = RTFact;
		this.dtdh = dtdh;
		this.dddp = dddp;
		this.dtel = dtel;
		this.tstar = tstar;
		this.corr = corr;
		this.corrError = corrError;
    }
    
    /** phase that generated this arrival. */
    private SeismicPhase phase;

    /** travel time in seconds */
    private double time;

    /** angular distance (great circle) in radians */
    private double dist;

    /** ray parameter in seconds per radians. */
    private double rayParam;

    private int rayParamIndex;

    /** phase name */
    private String name;

    /** phase name changed for true depths */
    private String puristName;

    /** source depth in kilometers */
    private double sourceDepth;

    /** pierce and path points */
    private TimeDist[] pierce, path;

    private double incidentAngle;
    
    private double takeoffAngle;
    
    //SH
    private double RTFact;
    
    private String polarization = "V";
    
    private String component = "Z";
    
    private double amplFact;
    
    private double dtdh;
    
    private double dddp;
    
    private double dtel;
    
    private double corr;
    
    private double corrError;
    
    private double[] tau = new double[3];
    
    private double tstar;
    
    // get set methods
    /** @return the phase used to calculate this arrival. */
    public SeismicPhase getPhase() {
        return phase;
    }

    /** @return travel time in seconds */
    public double getTime() {
        return time;
    }
    
    /**@return travel time as a Duration */
    public Duration getDuration() {
        return Duration.ofNanos(Math.round(getTime()*1000000000));
    }

    /** returns travel distance in radians */
    public double getDist() {
        return dist;
    }
    
    //SH
    public double getAmplFact() {
		return this.amplFact;
    }
    
    public double getRTFact() {
		return this.RTFact;
    }
    
    public String getPol() {
		return this.polarization;
    }
    
    public String getComp() {
		return this.component;
    }
    
    public double getDtdh() {
    	return this.dtdh;
    }
    
    public double getDddp() {
    	return this.dddp;
    }
    
    public double[] getTaus() {
    	return this.tau;
    }
    
    public double getDtel() {
    	return this.dtel;
    }
    
    public double getTstar() {
    	return this.tstar;
    }
    
    public double getCorr() {
	return this.corr;
    }
    
    public double getCorrError() {
    	return this.corrError;
    }
    
    public void setDtdh(double dtdh) {
    	this.dtdh = dtdh;
    }

    public void setDddp(double dddp) {
    	this.dddp = dddp;
    }
    
    public void setCorr(double corr) {
    	this.corr = corr;
    }
    
    public void setCorrError(double corrError) {
    	this.corrError = corrError;
    }
    
    public void setTaus(double[] tau) {
    	for (int i=0; i<3; i++) {
    		this.tau[i] = tau[i];
    	}
    }
    public void setTaus(double tau0, double tau1, double tau2) {
		this.tau[0] = tau0;
		this.tau[1] = tau1;
		this.tau[2] = tau2;
    }
    
    public void setDtel(double dtel) {
    	this.dtel = dtel; 
    }
    
    public void setTstar(double tstar) {
    	this.tstar = tstar; 
    }
    
    public void setAmplFact(double amplFact) {
		this.amplFact = amplFact;
    }
    
    public void setRTFact(double RTFact) {
		this.RTFact = RTFact;
    }
    
    public void setPol(String pol) {
    	this.polarization = pol;
    }
    
    public void setComp(String comp) {
		this.component = comp;
    }
    
    /**
     * returns travel distance in degrees.
     */
    public double getDistDeg() {
        return RtoD * getDist();
    }

    /**
     * returns distance in radians and in the range 0-PI. Note this may not be
     * the actual distance traveled.
     */
    public double getModuloDist() {
        double moduloDist = getDist() % TWOPI;
        if(moduloDist > Math.PI) {
            moduloDist = TWOPI - moduloDist;
        }
        return moduloDist;
    }

    /**
     * returns distance in degrees and in the range 0-180. Note this may not be
     * the actual distance traveled.
     */
    public double getModuloDistDeg() {
        double moduloDist = (RtoD * getDist()) % 360;
        if(moduloDist > 180) {
            moduloDist = 360 - moduloDist;
        }
        return moduloDist;
    }
    
    public void flipDistance() {
    	this.dist = this.getDist() * -1.;
    }

    /** returns ray parameter in seconds per radian */
    public double getRayParam() {
        return rayParam;
    }

    /** returns ray parameter in seconds per deg */
    public double getRayParamDeg() {
        return getRayParam()/RtoD;
    }

    public double getIncidentAngle() {
        return incidentAngle;
    }
    
    public double getTakeoffAngle() {
        return takeoffAngle;
    }

    public int getRayParamIndex() {
        return rayParamIndex;
    }

    /** returns phase name */
    public String getName() {
        return name;
    }

    /**
     * returns purist's version of name. Depths are changed to reflect the true
     * depth of the interface.
     */
    public String getPuristName() {
        return puristName;
    }

    /** returns source depth in kilometers */
    public double getSourceDepth() {
        return sourceDepth;
    }

    /** returns pierce points as TimeDist objects. */
    public TimeDist[] getPierce() {
        if (pierce == null) {
            this.pierce = getPhase().calcPierceTimeDist(this).toArray(new TimeDist[0]);
        }
        return pierce;
    }

    /** returns pierce points as TimeDist objects. */
    public TimeDist[] getPath() {
        if (path == null) {
            this.path = getPhase().calcPathTimeDist(this).toArray(new TimeDist[0]);
        }
        return path;
    }

    public String toString() {
    	
        String desc =  Outputs.formatDistance(getModuloDistDeg()) + Outputs.formatDepth(getSourceDepth()) + "   " + getName()
                + "  " + Outputs.formatTime(getTime()) + "  " + Outputs.formatRayParam(Math.PI / 180.0 * getRayParam())
                + "  " + Outputs.formatDistance(getTakeoffAngle()) + " " + Outputs.formatDistance(getIncidentAngle())
                + " " + Outputs.formatDistance(getDistDeg())+" "+getRayParamIndex()
                + " " + Outputs.formatScientific(getDddp());
        if (getName().equals(getPuristName())) {
            desc += "   = ";
        } else {
            desc += "   * ";
        }
        desc += getPuristName();
        return desc;
    }

    public int getNumPiercePoints() {
        if(pierce != null) {
            return pierce.length;
        } else {
            return 0;
        }
    }

    public int getNumPathPoints() {
        if(path != null) {
            return path.length;
        } else {
            return 0;
        }
    }

    public TimeDist getPiercePoint(int i) {
        // don't check for i> length since we want an ArrayOutOfBounds anyway
        return pierce[i];
    }

    /**
     * finds the first pierce point at the given depth.
     * 
     * @throws ArrayIndexOutOfBoundsException
     *             if depth is not found
     */
    public TimeDist getFirstPiercePoint(double depth) {
        for(int i = 0; i < pierce.length; i++) {
            if(pierce[i].getDepth() == depth) {
                return pierce[i];
            }
        }
        throw new ArrayIndexOutOfBoundsException("No Pierce point found for depth "
                + depth);
    }

    /**
     * finds the last pierce point at the given depth.
     * 
     * @throws ArrayIndexOutOfBoundsException
     *             if depth is not found
     */
    public TimeDist getLastPiercePoint(double depth) {
        TimeDist piercepoint = null;
        for(int i = 0; i < pierce.length; i++) {
            if(pierce[i].getDepth() == depth) {
                piercepoint = pierce[i];
            }
        }
        if(piercepoint == null) {
            throw new ArrayIndexOutOfBoundsException("No Pierce point found for depth "
                    + depth);
        }
        return piercepoint;
    }

    public TimeDist getPathPoint(int i) {
        // don't check for i> length since we want an ArrayOutOfBounds anyway
        return path[i];
    }

    protected static final double TWOPI = 2.0 * Math.PI;

    protected static final double DtoR = Math.PI / 180.0;

    protected static final double RtoD = 180.0 / Math.PI;
}
