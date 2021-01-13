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
import java.util.concurrent.Future
import java.util.concurrent.*
import java.util.Date
import java.util.ArrayList
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "FindSaddlePoints"
def descriptiveName = "Find Saddle Points"
def description = "Locate saddle points in a digital elevation model."
def toolboxes = ["TerrainAnalysis"]

public class FindSaddlePoints implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public FindSaddlePoints(WhiteboxPluginHost pluginHost, 
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
			def helpFile = "FindSaddlePoints"
			sd.setHelpFile(helpFile)
		
			// Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.
			def pathSep = File.separator
			def scriptFile = pluginHost.getResourcesDirectory() + 
			 "plugins" + pathSep + "Scripts" + pathSep + 
			 "FindSaddlePoints.groovy"
			sd.setSourceFile(scriptFile)
			
			// add some components to the dialog
			sd.addDialogFile("Ridge line file", "Ridge Network File:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Valley line file", "Valley Network File:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("D8 flow pointer file", "D8 Flow Pointer File:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Digital elevation network file", "DEM File:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Output file", "Output File:", "save", "Whitebox Files (*.dep; *.shp), DEP, SHP", true, false)
			
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
	  	
		if (args.length != 5) {
			pluginHost.showFeedback("Incorrect number of arguments given to tool.")
			return
		}
		// read the input parameters
		String ridgesFile = args[0]
		String valleysFile = args[1]
		String pointerFile = args[2]
		String demFile = args[3]
		String outputFile = args[4]
		boolean isRasterOutput = true
		if (!outputFile.toLowerCase().endsWith(".dep")) {
			isRasterOutput = false
		}

		System.out.println(outputFile)
		
		int i, row, col, rN, cN, n, numChannelHeads
		double z, valleyVal, ridgeVal, pointerVal
		int headVal
		double flowDir, initialElev
		final double LnOf2 = Math.log(2.0)
		def dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ]
		int[] dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ]
		double[] inflowingVals = [ 16, 32, 64, 128, 1, 2, 4, 8 ]
        
		
		WhiteboxRaster ridges = new WhiteboxRaster(ridgesFile, "r")
		int rows = ridges.getNumberRows()
		int cols = ridges.getNumberColumns()
		double ridgesNoData = ridges.getNoDataValue()
		
		WhiteboxRaster valleys = new WhiteboxRaster(valleysFile, "r")
		if (valleys.getNumberColumns() != cols || 
		  valleys.getNumberRows() != rows) {
		  	pluginHost.showFeedback("Each of the input files must have the same dimensions, i.e. rows and columns")
		  	return
		}
		double valleysNoData = valleys.getNoDataValue()
		
		WhiteboxRaster pointer = new WhiteboxRaster(pointerFile, "r")
		if (pointer.getNumberColumns() != cols || 
		  pointer.getNumberRows() != rows) {
		  	pluginHost.showFeedback("Each of the input files must have the same dimensions, i.e. rows and columns")
		  	return
		}
		double pointerNoData = pointer.getNoDataValue()
		
		WhiteboxRaster dem = new WhiteboxRaster(demFile, "r")
		if (dem.getNumberColumns() != cols || 
		  dem.getNumberRows() != rows) {
		  	pluginHost.showFeedback("Each of the input files must have the same dimensions, i.e. rows and columns")
		  	return
		}
		double demNoData = dem.getNoDataValue()

		String tempFile
		if (isRasterOutput) {
			tempFile = outputFile.replace(".dep", "_temp.dep")
		} else {
			tempFile = outputFile.replace(".shp", "_temp.dep")
		}
		System.out.println(tempFile)
		WhiteboxRaster tempGrid = new WhiteboxRaster(tempFile, "rw", 
  		  	  ridgesFile, DataType.FLOAT, ridgesNoData)
  		tempGrid.isTemporaryFile = true

  		// Find the channel heads.
  		for (row = 0; row < rows; row++) {
  			for (col = 0; col < cols; col++) {
  				valleyVal = valleys.getValue(row, col)
  				if (valleyVal != valleysNoData && valleyVal != 0 &&
  				  pointer.getValue(row, col) != pointerNoData) {
  					n = 0;
                    for (i = 0; i < 8; i++) {
                    	rN = row + dY[i]
                    	cN = col + dX[i]
                    	valleyVal = valleys.getValue(rN, cN)
                    	if (valleyVal != 0 && valleyVal != valleysNoData &&
                    	  pointer.getValue(rN, cN) == inflowingVals[i]) { 
                    	  	n++; 
                    	}
                    }

                    if (n == 0) {
                    	// it's a channel head
                    	numChannelHeads++
                    	tempGrid.setValue(row, col, numChannelHeads)
                    }
  				}
  			}
  		}

  		double[] xCoords = new double[numChannelHeads + 1]
  		double[] yCoords = new double[numChannelHeads + 1]
  		double[] minElev = new double[numChannelHeads + 1]

  		for (i = 0; i < numChannelHeads; i++) {
  			minElev[i] = Double.POSITIVE_INFINITY
  		}

//  		tempGrid.close();
//  		return

  		// descend each flowpath starting from a ridge cell
		// locating the flowpaths that terminate at channel
		// heads, keeping track of the flowpath with the 
		// lowest starting cell for each channel head.
		for (row = 0; row < rows; row++) {
  			for (col = 0; col < cols; col++) {
  				ridgeVal = ridges.getValue(row, col)
  				if (ridgeVal != ridgesNoData && ridgeVal != 0 &&
  				  pointer.getValue(row, col) != pointerNoData) {
  				  	z = dem.getValue(row, col)
  				  	
					// descend the flow path intiated at this
					// cell until you terminate at a valley
					boolean flag = true;
                    cN = col;
                    rN = row;
                    while (flag) {
                        flowDir = pointer.getValue(rN, cN);
                        if (flowDir > 0) {
                            i = (int) (Math.log(flowDir) / LnOf2);
                            cN += dX[i];
                            rN += dY[i];
                            valleyVal = valleys.getValue(rN, cN)
                            if (valleyVal != 0 && valleyVal != valleysNoData) {
                            	// you've hit a valley bottom
                                flag = false;
                                // is it a channel head?
                                headVal = (int)tempGrid.getValue(rN, cN)
                                if (headVal > 0) {
                                	if (z < minElev[headVal]) { 
                                		minElev[headVal] = z
                                		if (isRasterOutput) {
                                			xCoords[headVal] = col
                                			yCoords[headVal] = row
                                		} else {
                                			xCoords[headVal] = ridges.getXCoordinateFromColumn(col)
                                			yCoords[headVal] = ridges.getYCoordinateFromRow(row)
                                		}
                                	}
                                }
                            }
                        } else {
                        	// probably the edge of the grid or a pit
                            flag = false;
                        }
                    }
  				}
  			}
		}
		
		
		ridges.close()
		valleys.close()
		pointer.close()
		dem.close()
		tempGrid.close()

		if (isRasterOutput) {
			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", 
  		  	  ridgesFile, DataType.FLOAT, ridgesNoData)

  		  	for (i = 0; i < numChannelHeads; i++) {
				if (minElev[i] != Double.POSITIVE_INFINITY) {
					row = (int)yCoords[i]
					col = (int)xCoords[i]
					output.setValue(row, col, (double)i)
				}
  		  	}
  		  	
  		  	output.addMetadataEntry("Created by the "
                    + descriptiveName + " tool.")
        	output.addMetadataEntry("Created on " + new Date())
			output.close()
		
		} else {
			// set up the output files of the shapefile and the dbf
            DBFField[] fields = new DBFField[1];

            fields[0] = new DBFField();
            fields[0].setName("FID");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setFieldLength(10);
            fields[0].setDecimalCount(0);
            
            ShapeFile output = new ShapeFile(outputFile, ShapeType.POINT, fields);

            whitebox.geospatialfiles.shapefile.Point wbGeometry
            
			for (i = 0; i < numChannelHeads; i++) {
				if (minElev[i] != Double.POSITIVE_INFINITY) {
  					wbGeometry = new whitebox.geospatialfiles.shapefile.Point(xCoords[i], yCoords[i]);                  
                    Object[] rowData = new Object[1]
                    rowData[0] = new Double(i)
                    output.addRecord(wbGeometry, rowData);
				}
  			}

  			output.write()
		}

		// display the output image
		pluginHost.returnData(outputFile)

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
	def f = new FindSaddlePoints(pluginHost, args, descriptiveName)
}
