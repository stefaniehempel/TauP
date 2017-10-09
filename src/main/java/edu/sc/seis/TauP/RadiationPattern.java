package edu.sc.seis.TauP;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.InvalidClassException;
import java.io.OptionalDataException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.StreamCorruptedException;
import java.io.StreamTokenizer;
import java.util.List;

public class RadiationPattern extends TauP_Time {

	protected double strike = Double.MAX_VALUE, dip=strike, rake=strike;
	protected double radpatTakeOff = strike;
	protected double radpatAzimuth = strike;
	protected boolean tensor = true;
	protected double mrr = strike, mtt = strike, mpp=strike, mrt=strike, mrp=strike, mtp=strike;
	protected String psmecainfo = "";
	protected static boolean gmtScript = true;
	protected static boolean matlabData = false;
	protected static boolean map = false;
	protected static String mapFile = "taup_radpat_map";
	protected static String gmtFile = "taup_radpat.gmt";
	protected static String matlabFile = "taup_radpat.dat";
	protected static String psFile = "taup_radpat.ps";
	
	/**
	 * Constructors
	 */
	
	/** empty radPat */
	public RadiationPattern() {
	}
	
	/** read CMTSOLUTION from specified fileName */
	public RadiationPattern(String fileName) {
		readCMTSOLUTION(fileName);
	}
	
	/** create RadPat based on input of strike dip rake */
	public RadiationPattern(double strike, double dip, double rake) {
		this.strike = strike;
		this.dip = dip;
		this.rake = rake;
    }
	
	/** computing radiation pattern */
	public double getRadPatFactor(String phaseName) {
		if (tensor) return getRadPatFactorTensor(phaseName);
		else return getRadPatFactorSDS(phaseName);
	}
	
	public double getRadPatFactorTensor(String phaseName) {
		// code adapted after S. Chevrot, Septembre 2002
		double localAz = radpatAzimuth/180.*Math.PI;
		double localTOA = radpatTakeOff/180.*Math.PI;
		double [] gam=new double[3];
		double [] phi=new double[3];
		double xm1, xm2, xm3;
		double norm=//1;
		Math.sqrt(mrr*mrr+mtt*mtt+mpp*mpp+
				2.*(mrt*mrt+mrp*mrp+mtp*mtp))/Math.sqrt(2.);
		
		if (phaseName.equals("P") || phaseName.equals("p")) {
			  gam[0]=0.0-Math.cos(localTOA);
			  gam[1]=0.0-Math.sin(localTOA)*Math.cos(localAz);
			  gam[2]=Math.sin(localTOA)*Math.sin(localAz);
			  xm1=gam[0]*mrr+gam[1]*mrt+gam[2]*mrp;
			  xm2=gam[0]*mrt+gam[1]*mtt+gam[2]*mtp;
			  xm3=gam[0]*mrp+gam[1]*mtp+gam[2]*mpp;
			  return (gam[0]*xm1+gam[1]*xm2+gam[2]*xm3)/norm;
		}
		else if (phaseName.equals("SV") || phaseName.equals("sv")) {
			  gam[0]=0.0-Math.cos(localTOA);
			  gam[1]=0.0-Math.sin(localTOA)*Math.sin(localAz);
			  gam[2]=Math.sin(localTOA)*Math.sin(localAz);
			  phi[0]=Math.sin(localTOA);
			  phi[1]=0.0-Math.cos(localTOA)*Math.cos(localAz);
			  phi[2]=Math.cos(localTOA)*Math.sin(localAz);
			  xm1=gam[0]*mrr+gam[1]*mrt+gam[2]*mrp;
			  xm2=gam[0]*mrt+gam[1]*mtt+gam[2]*mtp;
			  xm3=gam[0]*mrp+gam[1]*mtp+gam[2]*mpp;
			  return (phi[0]*xm1+phi[1]*xm2+phi[2]*xm3)/norm;
		  }  
		  else if (phaseName.equals("SH")  || phaseName.equals("S") || phaseName.equals("s") ) {
			  gam[0]=0.0-Math.cos(localTOA);
			  gam[1]=0.0-Math.sin(localTOA)*Math.cos(localAz);
			  gam[2]=Math.sin(localTOA)*Math.sin(localAz);
			  phi[0]=0.;
			  phi[1]=Math.sin(localAz);
			  phi[2]=Math.cos(localAz);
			  xm1=gam[0]*mrr+gam[1]*mrt+gam[2]*mrp;
			  xm2=gam[0]*mrt+gam[1]*mtt+gam[2]*mtp;
			  xm3=gam[0]*mrp+gam[1]*mtp+gam[2]*mpp;
			  return (phi[0]*xm1+phi[1]*xm2+phi[2]*xm3)/norm;
		  }
		  else return 0.;
	}
	
	public double getRadPatFactorSDS(String phaseName) {
		// from Aki and Richards, 2002, p. 108 & 109
		
		double localAz = radpatAzimuth/180.*Math.PI;
		double localTOA = radpatTakeOff/180.*Math.PI;
		
		double dtr=Math.PI/180;
		double dphi = localAz-strike*dtr;
		
		if (phaseName.equals("P") || phaseName.equals("p")) {
			return Math.cos(rake*dtr)*
					Math.sin(dip*dtr)*
					Math.pow(Math.sin(localTOA),2)*
					Math.sin(2.*dphi)
					-Math.cos(rake*dtr)*Math.cos(dip*dtr)*Math.sin(2.*localTOA)*Math.cos(dphi)
					+Math.sin(rake*dtr)*Math.sin(2.*dip*dtr)*
					(Math.pow(Math.cos(localTOA),2)
					-Math.pow(Math.sin(localTOA),2)*Math.pow(Math.sin(dphi),2))
					+Math.sin(rake*dtr)*Math.cos(2.*dip*dtr)*Math.sin(2.*localTOA)*Math.sin(dphi);
		}
		else if (phaseName.equals("SV")) {
			return Math.sin(rake*dtr)*Math.cos(2.*dip*dtr)*Math.cos(2.*localTOA)*Math.sin(dphi)
					-Math.cos(rake*dtr)*Math.cos(dip*dtr)*Math.cos(2.*localTOA)*Math.cos(dphi)
					+0.5*Math.cos(rake*dtr)*Math.sin(dip*dtr)*Math.sin(2.*localTOA)*Math.sin(2.*dphi)
					-0.5*Math.sin(rake*dtr)*Math.sin(2.*dip*dtr)*Math.sin(2.*localTOA)*
					(1.+Math.pow(Math.sin(dphi),2));
		}
		else if (phaseName.equals("SH") || phaseName.equals("s") || phaseName.equals("S")) {
			return Math.cos(rake*dtr)*Math.cos(dip*dtr)*Math.cos(localTOA)*Math.sin(dphi)
					+Math.cos(rake*dtr)*Math.sin(dip*dtr)*Math.sin(localTOA)*Math.cos(2.*dphi)
					+Math.sin(rake*dtr)*Math.cos(2.*dip*dtr)*Math.cos(localTOA)*Math.cos(dphi)
					-0.5*Math.sin(rake*dtr)*Math.sin(2.*dip*dtr)*Math.sin(localTOA)*Math.sin(2.*dphi);
		}
		else return 0.;
	}

	/** read psmecavelo info */
	public void getPsMecaInfo(String fileName) {
		try {
	        String line = null;
	        
            FileReader fileReader = 
                new FileReader(fileName);

            BufferedReader bufferedReader = 
                new BufferedReader(fileReader);
            
            // just takes the last line of the file
            while((line = bufferedReader.readLine()) != null) {
            	String[] array = line.trim().split("\\s+");
            	try {
            		// test if first 10 elements are convertible to number, better ways to do this?
        			double tmp;
        			for (int i=0; i<array.length-1; i++) {
        				tmp = Double.parseDouble(array[i]);
        			}
        			// try if last element is convertible to number
        			try {
        				tmp = Double.parseDouble(array[array.length-1]);
        				psmecainfo = line;
        			} catch(NumberFormatException ex) { // otherwise add sourceDepth dummy
            			psmecainfo += array[0] +" "+ array[1] +" 0. ";
            			for (int i=2; i<array.length-1; i++) {
            				psmecainfo += array[i] + " ";
            			}
            			psmecainfo += "0. 0.";
        			}
        			if (eventLat==Double.MAX_VALUE) {
        				eventLon=Double.parseDouble(array[0]);
        				eventLat=Double.parseDouble(array[1]);
        			}
        		} catch(NumberFormatException ex2) {
        			throw new NumberFormatException("psmeca-file format not recognized - should be\n"
        					+ "lon lat depth str1 dip1 slip1 str2 dip2 slip2 mant exp plon plat... or\n"
        					+ "lon lat str1 dip1 rake1 str2 dip2 rake2 sc iexp name\n");
        		}
            }
            bufferedReader.close();
        }
        catch(FileNotFoundException ex) {
            if (DEBUG) System.out.println(
                "Unable to open file '" + 
                fileName + "'");  
            System.err.println("Caught FileNotFoundException: " + ex.getMessage());
            ex.printStackTrace();
        }
        catch(IOException ex) {
            if (DEBUG) System.out.println(
                "Error reading file '" 
                + fileName + "'");                  
            System.err.println("Caught IOException: " + ex.getMessage());
            ex.printStackTrace();
        }
	}
	
	/** create fake psmeca (new format) from CMTSOLUTION */
	public void getPsMecaNewInfo(String fileName) {
		if (fileName.isEmpty()) fileName="CMTSOLUTION";
		try {
	        String line = null;
	        double tmp, order=0.;
	        int i=0;
	        
            FileReader fileReader = 
                new FileReader(fileName);

            BufferedReader bufferedReader = 
                new BufferedReader(fileReader);
            
            while((line = bufferedReader.readLine()) != null) {
        		String[] array = line.trim().split("\\s+");
        		if (i==0) {
        			psmecainfo += array[5] +" "+ array[6] +" "+ array[7];
        			i++;
        		}
        		if (array[0].startsWith("M")) {
        			tmp = Double.parseDouble(array[1]);
        			if (i==1) order=Math.floor(Math.log10(Math.abs(tmp)/10));
        			psmecainfo += " " + tmp/Math.pow(10, order);
        			i++;
        		}
            }
            psmecainfo += " " + order + " X Y NAME";
           
            bufferedReader.close();
        }
        catch(FileNotFoundException ex) {
            if (DEBUG) System.out.println(
                "Unable to open file '" + 
                fileName + "'");  
            System.err.println("Caught FileNotFoundException: " + ex.getMessage());
            ex.printStackTrace();
        }
        catch(IOException ex) {
            if (DEBUG) System.out.println(
                "Error reading file '" 
                + fileName + "'");                  
            System.err.println("Caught IOException: " + ex.getMessage());
            ex.printStackTrace();
        }
	}
	
	public void readCMTSOLUTION() {
		readCMTSOLUTION("CMTSOLUTION");
	}
	
	public void readCMTSOLUTION(String fileName) {
		if (fileName.isEmpty()) fileName="CMTSOLUTION";
		try {

	        String line = null;
	        
            FileReader fileReader = 
                new FileReader(fileName);

            BufferedReader bufferedReader = 
                new BufferedReader(fileReader);
            
            while((line = bufferedReader.readLine()) != null) {
            		String[] array = line.split(":");
            		if (array[0].equals("Mrr")) this.mrr=Double.parseDouble(array[1]);
            		if (array[0].equals("Mtt")) this.mtt=Double.parseDouble(array[1]);
            		if (array[0].equals("Mpp")) this.mpp=Double.parseDouble(array[1]);
            		if (array[0].equals("Mrt")) this.mrt=Double.parseDouble(array[1]);
            		if (array[0].equals("Mrp")) this.mrp=Double.parseDouble(array[1]);
            		if (array[0].equals("Mtp")) this.mtp=Double.parseDouble(array[1]);
            }
           
            bufferedReader.close();   
            this.tensor = true;
        }
        catch(FileNotFoundException ex) {
            if (DEBUG) System.out.println(
                "Unable to open file '" + 
                fileName + "'");  
            System.err.println("Caught FileNotFoundException: " + ex.getMessage());
            ex.printStackTrace();
        }
        catch(IOException ex) {
            if (DEBUG) System.out.println(
                "Error reading file '" 
                + fileName + "'");                  
            System.err.println("Caught IOException: " + ex.getMessage());
            ex.printStackTrace();
        }
	}
	
	public boolean hasValue() {
		double norm = Math.sqrt(mrr*mrr+mtt*mtt+mpp*mpp+
						2.*(mrt*mrt+mrp*mrp+mtp*mtp))/Math.sqrt(2.);
		if (this.strike==0. && norm==0.) {
			return false;
		}
		else return true;
	}
	
	private void setStrikeDipRake(double strike, double dip, double rake) {
		this.strike=strike;
		this.dip=dip;
		this.rake=rake;
		this.tensor = false;
	}
	
	private void setMoment(double mrr, double mtt, double mpp, double mrt, double mrp, double mtp) {
		this.mrr = mrr;
		this.mtt = mtt;
		this.mpp = mpp;
		this.mrt = mrt;
		this.mrp = mrp;
		this.mtp = mtp;
		this.tensor = true;
	}
	
	public void setAzimuth(double radpatAzimuth) {
		this.radpatAzimuth=radpatAzimuth;
	}
	public void setTakeOffAngle(double radpatTakeOff) {
		this.radpatTakeOff=radpatTakeOff;
	}
	public void setDirections(double radpatAzimuth, double radpatTakeOff) {
		this.radpatAzimuth=radpatAzimuth;
		this.radpatTakeOff=radpatTakeOff;
	}
	
	public void clearRadPat() {
		strike = 0.;
		dip = 0.;
		rake = 0.;
		radpatTakeOff = 0.;
		radpatAzimuth = 0.;
		tensor = true;
		mrr = 0.;
		mtt = 0.;
		mpp = 0.;
		mrt = 0.;
		mrp = 0.;
		mtp = 0.;
	}
	
	// input and help
	public void printHelp() {
        Alert.info("Enter:" 
        		+ "-ph phase list                -- comma separated phase list\n"
        		+ "-sdr s d r                    -- for new strike/ dip/ rake input\n"
                + "-mom mrr mtt mpp mrt mrp mtp	 -- for new moment input\n"
                + "-mf filename                  -- to input filename with CMTsolution\n"
                + "-deg                          -- distance\n"
                + "-az                           -- for new radpatAzimuth\n"
                + "--gmt                         -- to plot 2D radpat\n"
                + "--dat                         -- to create output for 3D radpat plot e.g. in MatLab\n"
                + "--map                         -- create gmt-file with great circle arc and psmeca (needs sta input)\n"
                + "--sta lat lon                 -- station coordintes\n"
                + " or \nq to quit.\n");
    }
	
	public String[] parseCmdLineArgs(String[] leftOverArgs) throws IOException {
        int i = 0;
        int numNoComprendoArgs = 0;
        leftOverArgs = super.parseCmdLineArgs(leftOverArgs);
        String[] noComprendoArgs = new String[leftOverArgs.length];
        while(i < leftOverArgs.length) {
        	if(dashEquals("help", leftOverArgs[i])) {
                printHelp();
                noComprendoArgs[numNoComprendoArgs++] = leftOverArgs[i];
        	} else if(dashEquals("gmt", leftOverArgs[i])) {
                gmtScript = true;
            } else if(dashEquals("dat", leftOverArgs[i])) {
                matlabData = true;
            } else if(dashEquals("map", leftOverArgs[i])) {
                map = true;
            } else if(dashEquals("help", leftOverArgs[i])) {
                noComprendoArgs[numNoComprendoArgs++] = leftOverArgs[i];
            } else if(i < leftOverArgs.length - 3) {
            	if(dashEquals("sdr", leftOverArgs[i])) {
            		System.out.println("reading strike dip rake");
	            	setStrikeDipRake(Double.valueOf(leftOverArgs[i+1]).doubleValue(),Double.valueOf(leftOverArgs[i+2]).doubleValue(),
	            			Double.valueOf(leftOverArgs[i+3]).doubleValue());
	            	i = i+3;
            	} else if (i < leftOverArgs.length - 6 )
            		if(dashEquals("mom", leftOverArgs[i])) {
            			System.out.println("reading moment");
            			setMoment(Double.valueOf(leftOverArgs[i+1]).doubleValue(),Double.valueOf(leftOverArgs[i+2]).doubleValue(),
                    			Double.valueOf(leftOverArgs[i+3]).doubleValue(),Double.valueOf(leftOverArgs[i+4]).doubleValue(),
                    			Double.valueOf(leftOverArgs[i+5]).doubleValue(),Double.valueOf(leftOverArgs[i+6]).doubleValue());
                    	i = i+6;
            		} else {
		                noComprendoArgs[numNoComprendoArgs++] = leftOverArgs[i];
            		}
            	} else {
	                noComprendoArgs[numNoComprendoArgs++] = leftOverArgs[i];
	            }
            i++;
        }
        
        if(numNoComprendoArgs > 0) {
            String[] temp = new String[numNoComprendoArgs];
            System.arraycopy(noComprendoArgs, 0, temp, 0, numNoComprendoArgs);
            return temp;
        } else {
            return new String[0];
        }
    }
	
	public static void printHead( PrintWriter out ) throws IOException {
    	out.write("/bin/rm -f " + psFile + " \n");
		out.write("# draw surface and label distances.\n");
		out.write("psbasemap -K -P -R0/360/0/1 -JP6.0i/90. -B30p/500N > " + psFile + " \n");
		out.write("# draw circles for branches, note these are scaled for a\n"); 
		out.write("# map using -JP6.0i\n");
		out.write("psxy -K -O -P -R -JP -Sc -A >> " + psFile + " <<ENDLAYERS\n");
		out.write("0.0 0.0 6.0i\n");
    	out.write("ENDLAYERS\n\n");		    
	}	
	
	// test routines or all around patterns - by default based on CMTSOLUTION
	public void compute2DRadPat() throws IOException, TauModelException, SlownessModelException, VelocityModelException {
		//TBD include linestyle and legend for various phases 
		
        double[] takeOffs = new double[361];
        double[][] vals = new double[takeOffs.length][3];
        
        this.setAzimuth(azimuth);
        
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(gmtFile)));
        printHead( out );
        
        for (int i=0; i<361; i++) {
        	this.setTakeOffAngle(i);
        	vals[i][0]=this.getRadPatFactor("P");
        	vals[i][2]=this.getRadPatFactor("SV");
        	vals[i][1]=this.getRadPatFactor("SH");
        }
        printRadPatLine(out,vals,3);
        
        if (degrees!=Double.MAX_VALUE) {
        	calculate(degrees);
        	List<Arrival> arrivals = getArrivals();
	        for (Arrival currArrival : arrivals) {
	        	printTakeOffLines(out, currArrival.getTakeoffAngle());
	        }
        }
        
        out.close();
	}
	
	public void printTakeOffLines( PrintWriter out , double takeOffAngle) throws IOException {
		out.write("psxy -P -R -K -O -JP -A -W.5,black >> " + psFile + " <<END\n");
		out.write("0.,0.\n");
		out.write(takeOffAngle + " 1.\n");
		out.write("END\n\n");
	}
	
	public void compute3DRadPat() throws IOException, TauModelException, SlownessModelException, VelocityModelException {
		
        double[] takeOffs = new double[181];
        double[] azimuths = new double[361];
        double[][] vals = new double[azimuths.length*takeOffs.length][5];
        
        // open output file
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(matlabFile)));
        
        int cnt = 0;
        
        for (int i=0; i<361; i++) {
        	this.setAzimuth(i);
        	for (int j=0; j<181; j++) {
        		this.setTakeOffAngle(j);
        		vals[cnt][0]=i;
        		vals[cnt][1]=j;
        		vals[cnt][2]=this.getRadPatFactor("P");
        		vals[cnt][3]=this.getRadPatFactor("SV");
        		vals[cnt][4]=this.getRadPatFactor("SH");
        		cnt++;
        	}
        }
        printRadPatLineMatlab(out,vals,5);
        
        // compute arrivals to obtain takeOffAngle, rayParam, azimuth and degrees
        if (degrees!=Double.MAX_VALUE && azimuth!=Double.MAX_VALUE) {
	        calculate(degrees);
	        List<Arrival> arrivals = getArrivals();
	        for (Arrival currArrival : arrivals) {
	        	out.write("0 " + currArrival.getRayParamDeg() + " " + currArrival.getTakeoffAngle() + " " + azimuth + " " + degrees + "\n");
	        }
        } else if (degrees!=Double.MAX_VALUE || azimuth!=Double.MAX_VALUE) System.out.println("Warning: need both distance and azimuth information to return info on exiting rays!\n");
        
        out.close();
	}
	
	/** to be called by SeismicPhase */
	public double[] getOneRadpat(String firstleg, double takeoff, double az) {
		double[] radEnergy;
		
		this.setAzimuth(az);
		this.setTakeOffAngle(takeoff);
		
		if ( (firstleg.equalsIgnoreCase("P")) ||
				(firstleg.equalsIgnoreCase("SV")) ||
				(firstleg.equalsIgnoreCase("SH"))
				) {
    		radEnergy = new double[1];
    		radEnergy[0] = this.getRadPatFactor(firstleg);
    		}
    	else {
    		radEnergy = new double[2];
    		radEnergy[0] = this.getRadPatFactor("SV");
    		radEnergy[1] = this.getRadPatFactor("SH");
    	}
		
		return radEnergy;
	}
	
	public void mapView() throws IOException {
		
		if ( degrees == Double.MAX_VALUE && stationLat == Double.MAX_VALUE ) {
			throw new RuntimeException("Define receiver location to plot map!");
		}
		else {
			File f = new File(TauP_Time.cmtFileName);
			if(!f.exists()) { 
				f = new File("PSMECAVELO");
			}
			if(f.exists()) this.getPsMecaInfo(f.getName());
			else if (eventLat == Double.MAX_VALUE) throw new RuntimeException("Define source location to plot map!");
			//otherwise just plot star
			
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(mapFile + ".gmt")));
			// decide for min and max coordinates
			// for now a problem w date line TBD
			double minlat = Math.min(stationLat, eventLat);
			double maxlat = Math.max(stationLat, eventLat);
			double minlon = Math.min(stationLon, eventLon);
			double maxlon = Math.max(stationLon, eventLon);
			double mapWidth = Math.max(20., maxlon-minlon);
			double mapHeight = Math.max(20., maxlat-minlat);
			minlon -= .2*mapWidth; maxlon += .2*mapWidth;
			minlat -= .2*mapHeight; maxlat += .2*mapHeight;
			
			out.println("gmt psbasemap -R" + minlon+"/"+maxlon+"/"+minlat+"/"+maxlat + " -Jh0.3i -Bafg -K > " + mapFile + ".ps");
			// coast line
			out.println("gmt pscoast -R -J -W.1 -K -O >> " + mapFile + ".ps");
			// great circle arc
			out.println("gmt psxy -R -J -O -K -W.6,black << END >> " + mapFile + ".ps");
			out.println(stationLon + "," + stationLat);
			out.println(eventLon + "," + eventLat);
			out.println("END");
			// plot station
			out.println("gmt psxy -R -J -O -K -St0.4 -Gblack << END >> " + mapFile + ".ps");
			out.println(stationLon + "," + stationLat);
			out.println("END");
			// load this from file
			if (psmecainfo.isEmpty()) {
				out.println("gmt psxy -R -J -O -K -Sa0.4 -Gblack << END >> " + mapFile + ".ps");
				out.println(eventLon + "," + eventLat);
			} else {
				// plot focal mechanism
				out.println("gmt psmeca -R -J -Sc0.4 -h1 -K -O -Wblack << END >> " + mapFile + ".ps");
				out.println("lon lat depth str dip slip st dip slip mant exp plon plat");
				out.println(psmecainfo);
			}
			out.println("END");
			
			out.close();
		}
	}
	
	public void returnSpecificList() throws TauModelException, IOException, SlownessModelException, VelocityModelException {
		// produce a list of arrivals for this source receiver setup
//        if (stationLat < 200 && eventLat < 200) {
//        	degrees = SphericalCoords.distance(eventLat,eventLon,stationLat,stationLon);
//        	azimuth = SphericalCoords.azimuth(eventLat,eventLon,stationLat,stationLon);
//        }
        String firstleg;
        double[] radEnergy = new double[1];
        
        PrintWriter out = new PrintWriter(new OutputStreamWriter(System.out));
        
        int maxNameLength = 5;
        int maxPuristNameLength = 5;
        for(int j = 0; j < arrivals.size(); j++) {
            if(((Arrival)arrivals.get(j)).getName().length() > maxNameLength) {
                maxNameLength = ((Arrival)arrivals.get(j)).getName()
                        .length();
            }
            if(((Arrival)arrivals.get(j)).getPuristName().length() > maxPuristNameLength) {
                maxPuristNameLength = ((Arrival)arrivals.get(j)).getPuristName()
                        .length();
            }
        }
        Format phaseFormat = new Format("%-" + maxNameLength + "s");
        Format phasePuristFormat = new Format("%-" + maxPuristNameLength + "s");
        
        String[] polString = {"V","H"};
        
        calculate(degrees);
        List<Arrival> arrivals = getArrivals();
        for (Arrival currArrival : arrivals) {
            firstleg = currArrival.getName().substring(0,1);
        	radEnergy = getOneRadpat(firstleg, currArrival.getTakeoffAngle(), azimuth);
        	
        	for (int i=0; i<radEnergy.length; i++) {
	        	out.print(Outputs.formatDistance(currArrival.getModuloDistDeg())
	                    + Outputs.formatDepth(depth) + "   ");
	            if (radEnergy.length>1) out.print(phaseFormat.form(currArrival.getName() + " " + polString[i]));
	            else out.print(phaseFormat.form(currArrival.getName()));	            
	            out.print("  "
	                    + Outputs.formatTime(currArrival.getTime())
	                    + "  "
	                    + Outputs.formatRayParam(Math.PI / 180.0
	                            * currArrival.getRayParam()) + "  ");
	            out.print(Outputs.formatDistance(currArrival.getTakeoffAngle())+" ");
	            out.print(Outputs.formatScientific(radEnergy[i]));
	            out.println();
        	}
        }
        out.flush();
	}
	
	public static void printRadPatLine(PrintWriter out, double[][] vals, int columns) {
		double test = vals.length;
		String[] gmtColor = new String[3]; gmtColor[0]="blue"; gmtColor[1]="red"; gmtColor[2]="magenta";   
		for (int j=0; j<columns; j++) {
			out.write("psxy -P -R -K -O -JP -A -W2," + gmtColor[j] + " >> " + psFile + " <<END\n");
			for (int i=0; i<vals.length; i++) {
				out.write("  " + i + "  " + Math.abs(vals[i][j]) + "\n");
			}
			out.write("END\n\n");
		}
	}
	
	public static void printRadPatLineMatlab(PrintWriter out, double[][] vals, int columns) {
		for (int i=0; i<vals.length; i++) {
			for (int j=0; j<columns; j++) {
			out.write(vals[i][j] + " ");
			}
			out.write("\n");		
		}
	}
	
	public void start() throws TauModelException, TauPException, IOException, SlownessModelException, VelocityModelException {
		if((degrees != Double.MAX_VALUE || takeoffAngle != Double.MAX_VALUE 
                || (stationLat != Double.MAX_VALUE
                && stationLon != Double.MAX_VALUE
                && eventLat != Double.MAX_VALUE && eventLon != Double.MAX_VALUE))) {
            /* enough info given on cmd line, so just do one calc. */
            if (takeoffAngle != Double.MAX_VALUE) {
                calcTakeoff(takeoffAngle);
            } else {
                if(degrees == Double.MAX_VALUE) {
                    degrees = SphericalCoords.distance(stationLat,
                                                       stationLon,
                                                       eventLat,
                                                       eventLon);
                    azimuth = SphericalCoords.azimuth(eventLat,
                                                      eventLon,
                                                      stationLat,
                                                      stationLon);
                    backAzimuth = SphericalCoords.azimuth(stationLat,
                                                          stationLon,
                                                          eventLat,
                                                          eventLon);
                }
                depthCorrect(depth, getReceiverDepth());
            }
        }
		radpatStart();
	}
	
	public void radpatStart() throws IOException, TauModelException, SlownessModelException, VelocityModelException {

		if (strike==Double.MAX_VALUE && mrr==Double.MAX_VALUE) {
			this.clearRadPat();
			this.readCMTSOLUTION(TauP_Time.cmtFileName);
		}
		
        if (matlabData) compute3DRadPat();
        else if (map) mapView();
        else if (gmtScript) {
        	if (azimuth!=Double.MAX_VALUE) compute2DRadPat();
        	else throw new RuntimeException("Need to specify azimuth to plot 2D radiation pattern");
        }
//        else if (screen) returnSpecificList();
    }
	
	public static void main(String[] args) throws FileNotFoundException,
    IOException, StreamCorruptedException, ClassNotFoundException,
    OptionalDataException, SlownessModelException, VelocityModelException {
        RadiationPattern me = new RadiationPattern();
        try {
        	String[] noComprendoArgs = me.parseCmdLineArgs(args);
        	printNoComprendoArgs(noComprendoArgs);
        	me.init();
            me.start();
        } catch(TauModelException e) {
            System.err.println("Caught TauModelException: " + e.getMessage());
            e.printStackTrace();
        } catch(TauPException e) {
            System.err.println("Caught TauModelException: " + e.getMessage());
            e.printStackTrace();
        } catch(IOException e) {
            System.err.println("Caught IOException: " + e.getMessage());
            e.printStackTrace();
        }
    }

}