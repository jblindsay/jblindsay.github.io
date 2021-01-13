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
//import java.util.concurrent.Future
//import java.util.concurrent.*
import java.util.Date
import java.util.ArrayList
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "TraceDownslopeFlowpaths"
def descriptiveName = "Trace Downslope Flowpaths"
def description = "Traces downslope flowpaths from one or more target sites."
def toolboxes = ["HydroTools", "FlowpathTAs"]

public class TraceDownslopeFlowpaths implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public TraceDownslopeFlowpaths(WhiteboxPluginHost pluginHost, 
		String[] args, def name, def descriptiveName) {
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
			sd.setHelpFile(name)
		
			// Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.
			def pathSep = File.separator
			def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + name + ".groovy"
			sd.setSourceFile(scriptFile)
			
			// add some components to the dialog
			sd.addDialogFile("Input D8 flow pointer raster", "Input D8 Pointer Raster:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Input seed points file", "Input Seed Points File:", "open", "Whitebox Files (*.dep; *.shp), DEP, SHP", true, false)
            sd.addDialogFile("Output file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogCheckBox("Set the background value in the output image to NoData?", "Set background to NoData?", true)
			
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
	  		int progress, oldProgress, colN, rowN, c
	  		int[] dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ]
			int[] dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ]
        	double flowDir
        	boolean flag
        	double outputValue = 1.0
        	final double LnOf2 = 0.693147180559945;
        	
			if (args.length != 4) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String inputFile = args[0]
			String inputSeedFile = args[1]
			String outputFile = args[2]
			boolean backgroundNoData = Boolean.parseBoolean(args[3])
			
			// read the input image
			WhiteboxRaster pntr = new WhiteboxRaster(inputFile, "r")
			double nodata = pntr.getNoDataValue()
			int rows = pntr.getNumberRows()
			int cols = pntr.getNumberColumns()

			double backgroundValue
			if (backgroundNoData) {
				backgroundValue = nodata
			} else {
				backgroundValue = 0.0
			}

			WhiteboxRaster output
			if (backgroundNoData) {
				output = new WhiteboxRaster(outputFile, "rw", 
  		  	     inputFile, DataType.FLOAT, nodata)
			} else {
				output = new WhiteboxRaster(outputFile, "rw", 
  		  	     inputFile, DataType.FLOAT, 0.0)
			}

			if (inputSeedFile.toLowerCase().endsWith(".dep")) {
				// does the seed file have the same dimensions as the pointer?
				WhiteboxRasterInfo seeds = new WhiteboxRasterInfo(inputSeedFile)
				double seedNoData = seeds.getNoDataValue()
				if (rows != seeds.getNumberRows() ||
				  cols != seeds.getNumberColumns()) {
					pluginHost.showFeedback("Error: The seed file must have the same dimensions as the D8 pointer file.")
					return;
				}

				oldProgress = -1
				for (int row in 0..(rows - 1)) {
					double[] data = seeds.getRowValues(row);
	  				for (int col in 0..(cols - 1)) {
	  					if (data[col] != seedNoData && data[col] > 0) {
	  						// it's a seed point
	  						outputValue = data[col]
	  						output.setValue(row, col, outputValue)
	  						flag = false;
	                        colN = col;
	                        rowN = row;
	                        while (!flag) {
	                            // find its downslope neighbour
	                            flowDir = pntr.getValue(rowN, colN);
	                            if (flowDir > 0 && flowDir != nodata) {
	                            	//move x and y accordingly
	                                c = (int) (Math.log(flowDir) / LnOf2);
	                                colN += dX[c];
	                                rowN += dY[c];
	                                //if the new cell already has a value in the output, use that as the outletID
	                                if (output.getValue(rowN, colN) != backgroundValue) {
	                                    flag = true;
	                                } else {
	                                	output.setValue(rowN, colN, outputValue)
	                                }
	                            } else {
	                                flag = true;
	                            }
	                        }
	  					}
	  				}
	  				progress = (int)(100f * row / rows)
					if (progress > oldProgress) {
						pluginHost.updateProgress(progress)
						oldProgress = progress
					}
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
				
				seeds.close()
				
			} else if (inputSeedFile.toLowerCase().endsWith(".shp")) {
				def input = new ShapeFile(inputSeedFile)
	            ShapeType shapeType = input.getShapeType()
				if (shapeType.getBaseType() != ShapeType.POINT) {
					pluginHost.showFeedback("Error: The seed file must be of a POINT shape type.")
					return;
				}

				int featureNum = 0;
            	oldProgress = -1;
				int row, col;
            	int numPoints = input.getNumberOfRecords();
//            	double[][] points = new double[numPoints][2]
            	double[][] point

            	for (ShapeFileRecord record : input.records) {
            	 	point = record.getGeometry().getPoints()

					col = output.getColumnFromXCoordinate(point[0][0])
					row = output.getRowFromYCoordinate(point[0][1])

					outputValue = record.getRecordNumber()

					if (row >= 0 && row < rows && col >= 0 && col < cols) {
						// it's a seed point
						output.setValue(row, col, outputValue)
  						flag = false;
                        colN = col;
                        rowN = row;
                        
						while (!flag) {
                            // find its downslope neighbour
                            flowDir = pntr.getValue(rowN, colN);
                            if (flowDir > 0 && flowDir != nodata) {
                            	//move x and y accordingly
                                c = (int) (Math.log(flowDir) / LnOf2);
                                colN += dX[c];
                                rowN += dY[c];

                                //if the new cell already has a value in the output, use that as the outletID
                                if (output.getValue(rowN, colN) != backgroundValue) {
                                    flag = true;
                                } else {
                                	output.setValue(rowN, colN, outputValue)
                                }
                            } else {
                                flag = true;
                            }
                        }
					}
					
					progress = (int)(100f * featureNum / (numPoints - 1))
            		if (progress != oldProgress) {
						pluginHost.updateProgress(progress)
            			oldProgress = progress
            		}
            		// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
					featureNum++
            	}
            	
			} else {
				pluginHost.showFeedback("The seed file is not of a recognizable type.")
				return;
			}
			
			
			pntr.close()

			output.addMetadataEntry("Created by the "
	                    + descriptiveName + " tool.")
	        output.addMetadataEntry("Created on " + new Date())
			output.close()
	
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
	def tdf = new TraceDownslopeFlowpaths(pluginHost, args, name, descriptiveName)
}
