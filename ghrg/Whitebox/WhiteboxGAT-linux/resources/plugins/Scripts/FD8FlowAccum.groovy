/*
 * Copyright (C) 2015 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
import java.util.ArrayDeque
import java.text.DecimalFormat
import java.util.stream.IntStream
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.WhiteboxRasterBase.DataScale
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import groovy.transform.CompileStatic
import groovy.time.TimeDuration
import groovy.time.TimeCategory

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "FD8FlowAccum2"
def descriptiveName = "FD8 Flow Accumulation 2"
def description = "Calculates the FD8 flow accumulation raster."
def toolboxes = ["FlowAccum"]

public class FD8FlowAccum2 implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName

	public FD8FlowAccum2(WhiteboxPluginHost pluginHost, 
		String[] args, String name, String descriptiveName) {
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
			sd.addDialogFile("Input DEM raster", "Input DEM Raster:", "open", "Raster Files (*.dep), DEP", true, false)
            sd.addDialogFile("Output file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogDataInput("Exponent Parameter", "Exponent Parameter:", "1.1", true, false)
            sd.addDialogDataInput("Threshold for convergent flow (grid cells)", "Flow convergence threshold (grid cells):", "", true, true)
            
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

			Date start = new Date()
			
			int progress, oldProgress, i, row, col, rN, cN, dir, n
	  		int[] dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ]
			int[] dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ]
			int[] inflowingVals = [ 4, 5, 6, 7, 0, 1, 2, 3 ]
        	double z, zN, slope, maxSlope
        	DecimalFormat df = new DecimalFormat("###.#")
			
			if (args.length < 3) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String inputFile = args[0]
			String outputFile = args[1]
			double power = Double.parseDouble(args[2]);
			double convergenceThreshold = Double.POSITIVE_INFINITY
			if (!(args[3].toLowerCase().contains("not specified"))) {
				convergenceThreshold = Double.parseDouble(args[3])
			}
			
			// read the input image
			WhiteboxRaster dem = new WhiteboxRaster(inputFile, "r")
			dem.setForceAllDataInMemory(true)
			double nodata = dem.getNoDataValue()
			int rows = dem.getNumberRows()
			int rowsLessOne = rows - 1
			int cols = dem.getNumberColumns()
			int numPixels = rows * cols
			double gridResX = dem.getCellSizeX();
            double gridResY = dem.getCellSizeY();
            double diagGridRes = Math.sqrt(gridResX * gridResX + gridResY * gridResY);
            double[] gridLengths = [diagGridRes, gridResX, diagGridRes, gridResY, diagGridRes, gridResX, diagGridRes, gridResY]


			pluginHost.updateProgress("Creating output file:", 0)
			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", inputFile, DataType.FLOAT, 1d)
			output.setForceAllDataInMemory(true)
			output.setPreferredPalette("blueyellow.pal");
            output.setDataScale(DataScale.CONTINUOUS);
            output.setZUnits("dimensionless");
            output.setNonlinearity(0.2)
            
			// count the number of inflowing neighbours
			int numSolvedCells = 0;
			byte[][] numInflow = new byte[rows][cols]
			LinkedList<SimpleGridCell> stack = new LinkedList<SimpleGridCell>();
            for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					z = dem.getValue(row, col)
					if (z != nodata) {
						output.setValue(row, col, 1);
						n = 0
						for (i = 0; i < 8; i++) {
							rN = row + dY[i]
							cN = col + dX[i]
							zN = dem.getValue(rN, cN)
							if (zN != nodata && zN > z) {
								n++
							}
						}
						numInflow[row][col] = (byte)n;
						if (n == 0) {
							stack.push(new SimpleGridCell(row, col));
							numInflow[row][col] = (byte)-1;
							numSolvedCells++;
						}
					} else {
						numInflow[row][col] = (byte)-1;
						numSolvedCells++;
						output.setValue(row, col, nodata);
					}
				}
				progress = (int)(100f * row / rowsLessOne)
				if (progress != oldProgress) {
					pluginHost.updateProgress("Counting Inflowing Neighbours:", progress)
					oldProgress = progress

					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}

            SimpleGridCell gc;
            double[] weights = new double[8]
			boolean[] downslope = new boolean[8]
            int r, c
			double faValue, faValueN, totalWeights;
			oldProgress = -1
            while (!stack.isEmpty()) {
            	gc = stack.pop()
            	row = gc.row;
                col = gc.col;
                z = dem.getValue(row, col)
				faValue = output.getValue(row, col)

				// calculate the weights
				totalWeights = 0
				weights = new double[8]
				downslope = new boolean[8]
				if (faValue < convergenceThreshold) {
					for (i = 0; i < 8; i++) {
						zN = dem.getValue(row+dY[i], col+dX[i])
						if (zN < z && zN != nodata) {
							weights[i] = Math.pow(z-zN, power)
							totalWeights += weights[i]
							downslope[i] = true
						}
					}
				} else {
					// find the steepest downslope neighbour and give it all to them
            		dir = -1
					maxSlope = -99999999
					for (i = 0; i < 8; i++) {
						zN = dem.getValue(row + dY[i], col + dX[i])
						if (zN != nodata) {
							slope = (z - zN) / gridLengths[i]
							if (slope > 0) {
								downslope[i] = true
								if (slope > maxSlope) {
									maxSlope = slope
									dir = i
								}
							}
						}
					}
					weights[dir] = 1.0;
					totalWeights = 1.0;
				}
				if (totalWeights > 0) {
					// now perform the neighbour accumulation
					for (i = 0; i < 8; i++) {
						r = row + dY[i]
						c = col + dX[i]
						if (downslope[i]) {
							output.incrementValue(r, c, faValue*(weights[i]/totalWeights))
							numInflow[r][c]--; // = (byte)(numInflow[r][c] - 1);
							if (numInflow[r][c] == 0) {
								stack.push(new SimpleGridCell(r, c))
								numInflow[r][c] = (byte)-1
								
								numSolvedCells++
								progress = (int)(100f * numSolvedCells / numPixels)
								if (progress != oldProgress) {
									pluginHost.updateProgress("Accumulating flow:", progress)
									oldProgress = progress
				
									// check to see if the user has requested a cancellation
									if (pluginHost.isRequestForOperationCancelSet()) {
										pluginHost.showFeedback("Operation cancelled")
										return
									}
								}
							}
						}
					}
				}
            }

            dem.close();
            
            output.flush();

			Date stop = new Date()
			TimeDuration td = TimeCategory.minus( stop, start )
			
			output.addMetadataEntry("Created by the "
	                    + descriptiveName + " tool.")
	        output.addMetadataEntry("Created on " + new Date())
	        output.addMetadataEntry("Input file: $inputFile")
	        output.addMetadataEntry("Power: $power");
	        output.addMetadataEntry("Convergence Threshold: $convergenceThreshold")
			output.addMetadataEntry("Elapsed time: $td")
			output.close()

			pluginHost.showFeedback("Elapsed time: $td")
			
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

	@CompileStatic
    class SimpleGridCell {

        public int row;
        public int col;
        
        public SimpleGridCell(int row, int col) {
            this.row = row;
            this.col = col;
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
	def tdf = new FD8FlowAccum2(pluginHost, args, name, descriptiveName)
}
