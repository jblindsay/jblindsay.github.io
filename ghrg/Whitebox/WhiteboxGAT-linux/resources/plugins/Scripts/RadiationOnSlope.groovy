/*
 * Copyright (C) 2014 Rebecca Warren
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
import java.awt.event.ActionListener
import java.awt.event.ActionEvent
import java.util.concurrent.Future
import java.util.concurrent.*
import java.util.Date
import java.util.ArrayList
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
//def name = "Radiation on a Slope"
//def descriptiveName = "Radiation on a Slope"
//def description = "Calculated the radiation at a point on a sloping surface."
//def toolboxes = ["TerrainAnalysis"]

public class RadnSlope implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public RadnSlope(WhiteboxPluginHost pluginHost, 
		String[] args, def descriptiveName) {
		this.pluginHost = pluginHost
		this.descriptiveName = descriptiveName
			
		if (args.length > 0) {
			execute(args)
		} else {
			// Create a dialog for this tool to collect user-specified
			// tool parameters.
			sd = new ScriptDialog(pluginHost, descriptiveName, this)	
		
			// Specifying the help file will display the html help
			// file in the help pane. This file should be be located 
			// in the help directory and have the same name as the 
			// class, with an html extension.
			def helpFile = "RadSlope"
			sd.setHelpFile(helpFile)
		
			// Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.
			def pathSep = File.separator
			def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "RadiationOnSlope.groovy"
			sd.setSourceFile(scriptFile)
			
			// add some components to the dialog
			sd.addDialogFile("Specify the name of the digital elevation model raster file.", "Input Digital Elevation Model (DEM):", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Specify the name of the output file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogDataInput("Julian day (1 to 365)", "Julian Day:", "1", true, false)
			sd.addDialogDataInput("Latitude between +/- 66.5", "Latitude (Degrees):", "44.5", true, false)
			sd.addDialogDataInput("Dew point temperature for the day", "Dew Point Temperature (Degrees Celsius):", "-12", true, false)
			sd.addDialogDataInput("Average albedo of the surface", "Albedo:", "0.8", true, false)
			sd.addDialogDataInput("Factor to describe thickness of the atmosphere at a certain location", "Optical Air Mass:", "3.6", true, false)
			sd.addDialogDataInput("Coefficient for dust attenuation in the atmosphere", "Ydust:", "0.05", true, false)

			// resize the dialog to the standard size and display it
			sd.setSize(800, 400)
			sd.visible = true
		}
	}

	// The CompileStatic annotation can be used to significantly
	// improve the performance of a Groovy script to nearly 
	// that of native Java code.
	@CompileStatic
	private void execute(String[] args) {
	  try {
	  	
		if (args.length != 8) {
			pluginHost.showFeedback("Incorrect number of arguments given to tool.")
			return
		}

		def wd = pluginHost.getWorkingDirectory() 
		double z
		int progress
		int oldProgress = -1
  		
		// read the input parameters
		String demFile = args[0]
		String outputFile = args[1]
		int julianDay = Integer.parseInt(args[2])
		double latitude = Math.toRadians(Double.parseDouble(args[3]))
		double dewPoint = Double.parseDouble(args[4])
		double albedo = Double.parseDouble(args[5])
		double optAirMass = Double.parseDouble(args[6])
		double yDust = Double.parseDouble(args[7])
		

		WhiteboxRaster dem = new WhiteboxRaster(demFile, "r")
		int numRows = dem.getNumberRows()
		int numColumns = dem.getNumberColumns()
		double nodata = dem.getNoDataValue()

		// create output raster for radiation on a slope
		def KcsFile = wd + "Kcs.dep" 
		WhiteboxRaster radnOnSlope = new WhiteboxRaster(KcsFile, "rw", demFile, DataType.FLOAT, nodata)

		for (int row in 0..(numRows - 1)) {
  			for (int col in 0..(numColumns - 1)) {
  				z = dem.getValue(row, col)
  				radnOnSlope.setValue(row, col, z)
  			}
  			progress = (int)(100f * row / numRows)
			if (progress > oldProgress) {
				pluginHost.updateProgress(progress)
				oldProgress = progress
				// check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
			}
  		}

  		dem.close()


  		oldProgress = -1
		// get variables from DEM
			def zFactor = "1.0"
			// calculate slope
			def slopeFile = wd + "slope.dep" 
			args = [demFile, slopeFile, zFactor] 
			pluginHost.runPlugin("Slope", args, false, true)

			// calculate aspect
			def aspectFile = wd + "aspect.dep" 
			args = [demFile, aspectFile, zFactor] 
			pluginHost.runPlugin("Aspect", args, false, true)

			WhiteboxRaster slope = new WhiteboxRaster(slopeFile, "r")

			WhiteboxRaster aspect = new WhiteboxRaster(aspectFile, "r")

			
			// calculate total atmospheric transmissivity
					double t, tsa, Wp, asa, bsa, Ket, Eo, declination, dayAngle, eqLat, eqSlope, Tsr, Tss, Thr, Ths, slopeBeta, aspectAlpha, Kdif, Ys, Ts, As, Bs, KetHorPlane
					double Isc = 4.921
					double angularVelocity = 0.2618
			
					Wp = 1.12 * Math.exp(0.0614 * dewPoint)
					asa = -0.124 - (0.0207 * Wp)
					bsa = -0.0682 - (0.0248 * Wp)
					tsa = Math.exp(asa + (bsa  * optAirMass))
					t = tsa - yDust
					dayAngle = 2*Math.PI*(julianDay-1)/365
					Eo = 1.000110 + (0.034221 * Math.cos(dayAngle)) + (0.001280 * Math.sin(dayAngle)) + (0.000719 * Math.cos(2 * dayAngle)) + (0.000077 * Math.sin(2 * dayAngle))
					declination=(0.006918-0.399912*Math.cos(dayAngle)+0.070257*Math.sin(dayAngle)-0.006758*Math.cos(2*dayAngle)+(0.000907*Math.cos(3*dayAngle)+0.00148*Math.sin(3*dayAngle)))
					As= -0.0363 - (0.0084 * Wp)
					Bs= -0.0572 - (0.0173 * Wp)
					Ts = Math.exp(As + (Bs * optAirMass))
					Ys = 1 - Ts + yDust

//					pluginHost.returnData(String.valueOf(declination))

					for (int row in 0..(numRows - 1)) {
						for (int col in 0..(numColumns - 1)) {
							
							// get beta and alpha values from slope and aspect
							slopeBeta = Math.toRadians(slope.getValue(row, col))
							aspectAlpha = Math.toRadians(aspect.getValue(row, col))
								
							// calculate Longitude Difference (eqSlope) and Equivalent Latitude (eqLat)
							eqSlope = Math.atan(Math.sin(slopeBeta)*Math.sin(aspectAlpha)/(Math.cos(slopeBeta)*Math.cos(latitude)-Math.sin(slopeBeta)*Math.sin(latitude)*Math.cos(aspectAlpha)))
							eqLat=Math.asin(Math.sin(slopeBeta)*Math.cos(aspectAlpha)*Math.cos(latitude)+Math.cos(slopeBeta)*Math.sin(latitude))
		
							// calculate sunrise and sunset times (Thr/Tsr and Ths/Tss)
							Thr = -Math.acos(-Math.tan(declination) * Math.tan(latitude)) / angularVelocity
							Ths = Math.acos(-Math.tan(declination) * Math.tan(latitude)) / angularVelocity
		
							Tss=(Math.acos(-Math.tan(eqLat) * Math.tan(declination)) - eqSlope) / angularVelocity
							Tsr=(-Math.acos(-Math.tan(eqLat) * Math.tan(declination)) - eqSlope) / angularVelocity
					
		 					Tsr=Math.max(Tsr,Thr)
		 					Tss=Math.min(Tss,Ths)
		
		 					// calculate Extraterrestrial Radiation on Slope (Ket)
							Ket = Isc*Eo*(Math.cos(eqLat)*Math.cos(declination)*((Math.sin(0.2618*Tss+eqSlope)-Math.sin(0.2618*Tsr+eqSlope))/0.2618)+Math.sin(eqLat)*Math.sin(declination)*(Tss-Tsr))
					
							// calculate Diffuse Radiation (Kdif)
							Kdif = 0.5 * Ys * Ket
					
							// calculate Backscattered Radiation (Kbs)
							double Kbs, Kg, Kdir, Kcs
							Kdir = t * Ket
							Kg = Kdir + Kdif
							Kbs = 0.5 * Ys * albedo * Kg
		
							// calculate Total Incident Radiation on Slope (Kcs)
							Kcs = 11.575 * (Kg + Kbs)
							
							// fill raster with Kcs values
							radnOnSlope.setValue(row, col, Kcs)
								
							}
						}
	
		radnOnSlope.addMetadataEntry("Created by the "
                    + descriptiveName + " tool.")
        radnOnSlope.addMetadataEntry("Created on " + new Date())
		radnOnSlope.close()

		// display the output image
		pluginHost.returnData(KcsFile)

		// reset the progress bar
		pluginHost.updateProgress(0)
	  } catch (OutOfMemoryError oe) {
            pluginHost.showFeedback("An out-of-memory error has occurred during operation.")
      } catch (Exception e) {
            pluginHost.showFeedback("An error has occurred during operation. See log file for details.")
            pluginHost.logException("Error in " + descriptiveName, e)
      }
	}
	
	@Override
    public void actionPerformed(ActionEvent event) {
    	if (event.getActionCommand().equals("ok")) {
    		final def args = sd.collectParameters()
			sd.dispose()
			final Runnable r = new Runnable() {
            	@Override
            	public void run() {
                	execute(args)
            	}
        	}
        	final Thread t = new Thread(r)
        	t.start()
    	}
    }
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def f = new RadnSlope(pluginHost, args, descriptiveName)
}
