/*
 * Copyright (C) 2014 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
import java.util.Date
import java.util.ArrayList
import java.util.Queue
import java.util.LinkedList
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.WhiteboxRasterBase.DataScale
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "TopologicalStreamOrder"
def descriptiveName = "Topological Stream Order"
def description = "Assigns each link in a stream network its topological order."
def toolboxes = ["StreamAnalysis"]

public class TopologicalStreamOrder implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public TopologicalStreamOrder(WhiteboxPluginHost pluginHost, 
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
			sd.addDialogFile("Input streams raster", "Input Streams Raster:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Input D8 flow pointer raster", "Input D8 Pointer Raster:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Output file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
			
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
	  		int progress, oldProgress, x, y, nx, ny, x2, y2, c
	  		int[] dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ]
			int[] dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ]
			double[] inflowingVals = [16, 32, 64, 128, 1, 2, 4, 8]
        	double flowDir, nFlowDir, streamVal, maxFlowAccum, order
        	double outNodata = -32768
        	boolean flag, isOutlet
        	final double LnOf2 = 0.693147180559945
        	
			if (args.length != 3) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String streamsFile = args[0]
			String pointerFile = args[1]
			String outputFile = args[2]
			
			WhiteboxRaster streams = new WhiteboxRaster(streamsFile, "r")
			double nodata = streams.getNoDataValue()
			int rows = streams.getNumberRows()
			int cols = streams.getNumberColumns()
			
			WhiteboxRaster pntr = new WhiteboxRaster(pointerFile, "r")
			double pointerNodata = pntr.getNoDataValue()
			if (rows != pntr.getNumberRows() || cols != pntr.getNumberColumns()) {
				pluginHost.returnData("Error: All input files must have the same dimensions (rows and columns).")
				return
			}

			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", 
  		  	     streamsFile, DataType.INTEGER, outNodata, outNodata)
  		  	output.setPreferredPalette("spectrum.pal")

			Queue<ConfluenceCell> confluences = new LinkedList<>()

			oldProgress = -1
			for (int row in 0..(rows - 1)) {
				for (int col in 0..(cols - 1)) {
					streamVal = streams.getValue(row, col)
  					if (streamVal != nodata && streamVal > 0) {
  						isOutlet = false
                        flowDir = pntr.getValue(row, col)
                        if (flowDir == 0) { 
                            isOutlet = true
                        } else {
                            c = (int)(Math.log(flowDir) / LnOf2)
                            if (streams.getValue(row + dY[c], col + dX[c]) <= 0 
                                    || streams.getValue(row + dY[c], col + dX[c]) == nodata) {
                                isOutlet = true;
                            }
                        }

						order = 1
						if (isOutlet) {
                            x = col
                            y = row
                            flag = true
                            while (flag) {
                                output.setValue(y, x, order)
                                
                                nx = 0
                                ny = 0
                                boolean[] unsolvedInflows = new boolean[8]
                                int myInflow = -1
                                int numInflow = 0
                                for (c = 0; c < 8; c++) {
                                    x2 = x + dX[c]
                                    y2 = y + dY[c]
                                    nFlowDir = pntr.getValue(y2, x2)
                                    if (streams.getValue(y2, x2) > 0 
                                            && nFlowDir == inflowingVals[c]) {
                                        nx = x2
                                        ny = y2
                                        myInflow = c
                                        unsolvedInflows[c] = true
                                        numInflow++
                                    }
                                }
                                
                                if (numInflow > 0) {
                                	if (numInflow > 1) {
                                		order++
                                	}
                                	unsolvedInflows[myInflow] = false
                                    for (c = 0; c < 8; c++) {
                                    	if (unsolvedInflows[c]) {
                                    		confluences.add(new ConfluenceCell(y + dY[c], x + dX[c], order))
                                    	}
                                    }
                                    
                                	x = nx
                                    y = ny
                                } else {
                                    flag = false
                                }
                            }
                        }
  					}
  				}
  				progress = (int)(100f * row / rows)
				if (progress > oldProgress) {
					pluginHost.updateProgress("Loop 1 of 2:", progress)
					oldProgress = progress

					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}

			pluginHost.updateProgress("Loop 2 of 2:", 0)
			while (confluences.size() > 0) {
				ConfluenceCell confluence = confluences.poll()
				order = confluence.order
				
				x = confluence.column
                y = confluence.row
                flag = true
                while (flag) {
                    output.setValue(y, x, order)
                    
                    // find the upslope neighbouring stream cell
                    // with the highest flow accumulation value
                    maxFlowAccum = 0
                    nx = 0
                    ny = 0
                    boolean[] unsolvedInflows = new boolean[8]
                    int myInflow = -1
                    int numInflow = 0
                    for (c = 0; c < 8; c++) {
                        x2 = x + dX[c]
                        y2 = y + dY[c]
                        nFlowDir = pntr.getValue(y2, x2)
                        if (streams.getValue(y2, x2) > 0 
                                && nFlowDir == inflowingVals[c]) {
                            nx = x2
                            ny = y2
                            myInflow = c
                            unsolvedInflows[c] = true
                            numInflow++
                        }
                    }
                    
                    if (numInflow > 0) {
                    	if (numInflow > 1) {
                    		order++
                    	}
                    	unsolvedInflows[myInflow] = false
                        for (c = 0; c < 8; c++) {
                        	if (unsolvedInflows[c]) {
                        		confluences.add(new ConfluenceCell(y + dY[c], x + dX[c], order))
                        	}
                        }
                        
                    	x = nx
                        y = ny
                    } else {
                        flag = false
                    }
                }
			}

			streams.close()
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
        	pluginHost.updateProgress("Progress:", 0)
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

    class ConfluenceCell {
		int row
		int column
		double order

		public ConfluenceCell(int row, int column, double order) {
			this.row = row
			this.column = column
			this.order = order
		}
    	
    }
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def tdf = new TopologicalStreamOrder(pluginHost, args, name, descriptiveName)
}
