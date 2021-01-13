/*
 * Copyright (C) 2013 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
import java.lang.Math.*
import java.util.concurrent.Future
import java.util.concurrent.*
import java.util.Date
import java.util.ArrayList
import java.util.Map
import java.util.Random
import jopensurf.*;
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import whitebox.structures.XYPoint
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "EstimateHeightsFromParallax"
def descriptiveName = "Estimate Heights From Parallax"
def description = "Estimates heights from parallax"
def toolboxes = ["Photogrammetry"]

public class EstimateHeightsFromParallax implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public EstimateHeightsFromParallax(WhiteboxPluginHost pluginHost, 
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
			def helpFile = "EstimateHeightsFromParallax"
			sd.setHelpFile(helpFile)
		
			// Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.
			def pathSep = File.separator
			def scriptFile = pluginHost.getResourcesDirectory() + 
			 "plugins" + pathSep + "Scripts" + pathSep + 
			   "EstimateHeightsFromParallax.groovy"
			sd.setSourceFile(scriptFile)
			
			// add some components to the dialog
			sd.addDialogFile("Left image pricipal point", "Left Image Principal Point:", "open", "Vector Files (*.shp), SHP", true, false)
			sd.addDialogFile("Left image conjugate pricipal point", "Left Image Conjugate Principal Point:", "open", "Vector Files (*.shp), SHP", true, false)
			sd.addDialogFile("Left image tie points", "Left Image Tie Points:", "open", "Vector Files (*.shp), SHP", true, false)
			
			sd.addDialogFile("Right image pricipal point", "Right Image Principal Point:", "open", "Vector Files (*.shp), SHP", true, false)
			sd.addDialogFile("Right image conjugate pricipal point", "Right Image Conjugate Principal Point:", "open", "Vector Files (*.shp), SHP", true, false)
			sd.addDialogFile("Right image tie points", "Right Image Tie Points:", "open", "Vector Files (*.shp), SHP", true, false)

			sd.addDialogDataInput("Flying height above sea level", "Flying Height (ASL):", "948.6", true, false)
            sd.addDialogDataInput("Average ground elevation", "Average Ground Elevation:", "338.0", true, false)
            
			sd.addDialogDataInput("Output field name", "Output Field Name:", "HEIGHT", false, false)
            
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
	  	
			if (args.length != 9) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}

			int progress, oldProgress
			double x, y, p, leftX, rightX, leftY, rightY
			
			// read the input parameters
			String leftPPFile = args[0]
			String leftCPPFile = args[1]
			String leftTiePointsFile = args[2]

			String rightPPFile = args[3]
			String rightCPPFile = args[4]
			String rightTiePointsFile = args[5]
			double A = Double.parseDouble(args[6])
			double G = Double.parseDouble(args[7])
			double H = A - G
			String fieldName = args[8]

			double gridRes = 0.084609684
			double imageUnitConversion = 2.11197E-05

			// open the various shapefiles
			ShapeFile leftPP = new ShapeFile(leftPPFile)
			ShapeFile leftCPP = new ShapeFile(leftCPPFile)
			ShapeFile leftTiePoints = new ShapeFile(leftTiePointsFile)
			
			ShapeFile rightPP = new ShapeFile(rightPPFile)
			ShapeFile rightCPP = new ShapeFile(rightCPPFile)
			ShapeFile rightTiePoints = new ShapeFile(rightTiePointsFile)
			
			// now make sure that the they are all of a POINT base shapetype
			if (leftPP.getShapeType().getBaseType() != ShapeType.POINT) {
				pluginHost.showFeedback("ERROR: One of the input files is not of a POINT ShapeType.")
				return
			}
			if (leftCPP.getShapeType().getBaseType() != ShapeType.POINT) {
				pluginHost.showFeedback("ERROR: One of the input files is not of a POINT ShapeType.")
				return
			}
			if (leftTiePoints.getShapeType().getBaseType() != ShapeType.POINT) {
				pluginHost.showFeedback("ERROR: One of the input files is not of a POINT ShapeType.")
				return
			}
			if (rightPP.getShapeType().getBaseType() != ShapeType.POINT) {
				pluginHost.showFeedback("ERROR: One of the input files is not of a POINT ShapeType.")
				return
			}
			if (rightCPP.getShapeType().getBaseType() != ShapeType.POINT) {
				pluginHost.showFeedback("ERROR: One of the input files is not of a POINT ShapeType.")
				return
			}
			if (rightTiePoints.getShapeType().getBaseType() != ShapeType.POINT) {
				pluginHost.showFeedback("ERROR: One of the input files is not of a POINT ShapeType.")
				return
			}

			// make sure that the PP and CPP files only have one record in each
			if (leftPP.getNumberOfRecords() != 1) {
				pluginHost.showFeedback("ERROR: Input principal point file should only contain one feature.")
				return
			}
			if (leftCPP.getNumberOfRecords() != 1) {
				pluginHost.showFeedback("ERROR: Input conjugate principal point file should only contain one feature.")
				return
			}
			if (rightPP.getNumberOfRecords() != 1) {
				pluginHost.showFeedback("ERROR: Input principal point file should only contain one feature.")
				return
			}
			if (rightCPP.getNumberOfRecords() != 1) {
				pluginHost.showFeedback("ERROR: Input conjugate principal point file should only contain one feature.")
				return
			}

			// make sure that the tie points files have the same number of records
			if (leftTiePoints.getNumberOfRecords() != rightTiePoints.getNumberOfRecords()) {
				pluginHost.showFeedback("ERROR: Input tie points files should have an equal number of features.")
				return
			}

			int numTiePoints = leftTiePoints.getNumberOfRecords()

			// read in the PP and CPP coordinates
			double leftPPx, leftPPy, leftCPPx
			double leftCPPy, rightPPx, rightPPy
			double rightCPPx, rightCPPy
			
			double[][] point
			
			point = leftPP.getRecord(0).getGeometry().getPoints()
			leftPPx = point[0][0]
			leftPPy = point[0][1]
			leftPP = null

			point = leftCPP.getRecord(0).getGeometry().getPoints()
			leftCPPx = point[0][0]
			leftCPPy = point[0][1]
			leftCPP = null

			point = rightPP.getRecord(0).getGeometry().getPoints()
			rightPPx = point[0][0]
			rightPPy = point[0][1]
			rightPP = null

			point = rightCPP.getRecord(0).getGeometry().getPoints()
			rightCPPx = point[0][0]
			rightCPPy = point[0][1]
			rightCPP = null

			// calculate the left/right absolute parallax
			double leftP = Math.sqrt((leftPPy - leftCPPy) * (leftPPy - leftCPPy) + (leftPPx - leftCPPx) * (leftPPx - leftCPPx))
			double rightP = Math.sqrt((rightPPy - rightCPPy) * (rightPPy - rightCPPy) + (rightPPx - rightCPPx) * (rightPPx - rightCPPx))

			// calculate the average absolute parallax and the air base (B)
			double P = (leftP + rightP) / 2
			double B = P * gridRes
			
			// figure out the angles of the flight lines in each image
			double leftFlightLineAngle = -Math.atan2((leftCPPy - leftPPy), (leftCPPx - leftPPx))
			double rightFlightLineAngle = -Math.atan2((rightPPy - rightCPPy), (rightPPx - rightCPPx))
			
			// for each tie point, calculate the parallax and heights
			ShapeFileRecord leftRecord, rightRecord

			AttributeTable leftTable = leftTiePoints.getAttributeTable()
			AttributeTable rightTable = rightTiePoints.getAttributeTable()

			DBFField field = new DBFField()
            field = new DBFField()
            field.setName(fieldName)
            field.setDataType(DBFField.DBFDataType.NUMERIC)
            field.setFieldLength(10)
            field.setDecimalCount(4)
            leftTable.addField(field)
            rightTable.addField(field)

			double refP, h
			double[] hValues = new double[numTiePoints]
			double[] xParallax = new double[numTiePoints]
			double[] yParallax = new double[numTiePoints]
			double avgHValue = 0
			for (int r = 0; r < numTiePoints; r++) {
				leftRecord = leftTiePoints.getRecord(r)
				rightRecord = rightTiePoints.getRecord(r)
				
				point = leftRecord.getGeometry().getPoints()
				x = point[0][0]
				y = point[0][1]
				leftX = (x - leftPPx) * Math.cos(leftFlightLineAngle) - (y - leftPPy) * Math.sin(leftFlightLineAngle)
				leftY = (x - leftPPx) * Math.sin(leftFlightLineAngle) + (y - leftPPy) * Math.cos(leftFlightLineAngle)
				
				point = rightRecord.getGeometry().getPoints()
				x = point[0][0]
				y = point[0][1]
				rightX = (x - rightPPx) * Math.cos(rightFlightLineAngle) - (y - rightPPy) * Math.sin(rightFlightLineAngle)
				rightY = (x - rightPPx) * Math.sin(rightFlightLineAngle) + (y - rightPPy) * Math.cos(rightFlightLineAngle)
				
				p = Math.abs(leftX - rightX)

				if (r > 0) {
//				double h = H - B * f / (p * imageUnitConversion)
				 	h = H  * (p - refP) / (P + (p - refP))
				} else {
					refP = p
					h = 0
				}
				
				hValues[r] = h
				avgHValue += h
				xParallax[r] = p
				yParallax[r] = Math.abs(leftY - rightY)
				
//				Object[] recData = leftTable.getRecord(r)
//                recData[recData.length - 1] = new Double(h)
//                leftTable.updateRecord(r, recData)
//
//                rightTable.getRecord(r)
//                recData[recData.length - 1] = new Double(h)
//                rightTable.updateRecord(r, recData)
			}

			avgHValue = avgHValue / numTiePoints

			for (int r = 0; r < numTiePoints; r++) {

				h = G + (hValues[r] - avgHValue)

				Object[] recData = leftTable.getRecord(r)
                recData[recData.length - 1] = new Double(h)
                leftTable.updateRecord(r, recData)

                rightTable.getRecord(r)
                recData[recData.length - 1] = new Double(h)
                rightTable.updateRecord(r, recData)
			}

			leftTiePoints.write()
			rightTiePoints.write()

			pluginHost.showFeedback("Operation complete. The calculated height data has been written \nto the attribute tables of the two tie point files.")
			
	  	} catch (OutOfMemoryError oe) {
            pluginHost.showFeedback("An out-of-memory error has occurred during operation.")
	    } catch (Exception e) {
	        pluginHost.showFeedback("An error has occurred during operation. See log file for details.")
	        pluginHost.logException("Error in " + descriptiveName, e)
        } finally {
        	// reset the progress bar
        	pluginHost.updateProgress(0)
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
	def f = new EstimateHeightsFromParallax(pluginHost, args, descriptiveName)
}
