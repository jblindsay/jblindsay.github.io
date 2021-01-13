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
import java.beans.PropertyChangeEvent
import java.beans.PropertyChangeListener
import java.util.Date
import java.util.ArrayList
import java.util.PriorityQueue
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.WhiteboxRasterBase.DataScale
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.*
import whitebox.utilities.StringUtilities
import whitebox.structures.BooleanBitArray2D
import whitebox.structures.NibbleArray2D
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "FloodOrder"
def descriptiveName = "Flood Order"
def description = "Calculates the flood order."
def toolboxes = ["HydroTools", "TerrainAnalysis"]

public class FloodOrder implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public FloodOrder(WhiteboxPluginHost pluginHost, 
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
            sd.addDialogFile("Output file", "Output File:", "save", "Raster Files (*.dep), DEP", true, false)
			
			// resize the dialog to the standard size and display it
			sd.setSize(800, 400)
			sd.visible = true
		}
	}

	@CompileStatic
	private void execute(String[] args) {
		try {
	  		int progress, oldProgress, col, row, colN, rowN, numPits, r, c
	  		int numSolvedCells = 0
	  		int dir, numCellsInPath
	  		double z, zN, zTest, zN2, lowestNeighbour
	  		boolean isPit, isEdgeCell, flag
	  		//double SMALL_NUM = 0.0001d
	  		GridCell gc
			int[] dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ]
			int[] dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ]
			int[] backLink = [5, 6, 7, 8, 1, 2, 3, 4]
			double[] outPointer = [0, 1, 2, 4, 8, 16, 32, 64, 128]
			
			if (args.length < 2) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String demFile = args[0]
			String outputFile = args[1]
			
			// read the input image
			WhiteboxRaster dem = new WhiteboxRaster(demFile, "r")
			double nodata = dem.getNoDataValue()
			int rows = dem.getNumberRows()
			int cols = dem.getNumberColumns()
			int rowsLessOne = rows - 1
			int colsLessOne = cols - 1
			int numCellsTotal = rows * cols

			double minVal = dem.getMinimumValue()
			int elevDigits = (String.valueOf((int)(dem.getMaximumValue() - minVal))).length()
			double elevMultiplier = Math.pow(10, 8-elevDigits)
			double SMALL_NUM = 1 / elevMultiplier * 10
			long priority

			double[][] output = new double[rows + 2][cols + 2]
			double[][] floodOrder = new double[rows + 2][cols + 2]
			PriorityQueue<GridCell> queue = new PriorityQueue<GridCell>((2 * rows + 2 * cols) * 2);
			BooleanBitArray2D inQueue = new BooleanBitArray2D(rows + 2, cols + 2)
			
			// initialize the grids
			oldProgress = -1
			for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					z = dem.getValue(row, col)
					output[row + 1][col + 1] = z
					if (z != nodata) {
						isPit = true
						isEdgeCell = false
						lowestNeighbour = Double.POSITIVE_INFINITY
						for (int n in 0..7) {
							zN = dem.getValue(row + dY[n], col + dX[n])
							if (zN != nodata && zN < z) {
								isPit = false
								break
							} else if (zN == nodata) {
								isEdgeCell = true
							} else {
								if (zN < lowestNeighbour) { lowestNeighbour = zN }
							}
						}
						if (isPit) { 
							if (isEdgeCell) {
								//priority = (long)(z * elevMultiplier * 100000)
								queue.add(new GridCell(row + 1, col + 1, z, 0))
								inQueue.setValue(row + 1, col + 1, true)
							}
						}
					} else {
                        numSolvedCells++
                        floodOrder[row][col] = -1
                    }
				}
				progress = (int)(100f * row / rowsLessOne)
				if (progress > oldProgress) {
					pluginHost.updateProgress("Loop 1 of 3", progress)
					oldProgress = progress

					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}
			def paletteName = dem.getPreferredPalette()
			dem.close()

			for (row = 0; row < rows + 2; row++) {
				output[row][0] = nodata
				output[row][cols + 1] = nodata
			}
			
			for (col = 0; col < cols + 2; col++) {
				output[0][col] = nodata
				output[rows + 1][col] = nodata
			}

			// now breach
			oldProgress = (int) (100f * numSolvedCells / numCellsTotal);
            pluginHost.updateProgress("Loop 2 of 3: ", oldProgress);

			int order = 0
			int flatIndex
			int k
            while (queue.isEmpty() == false) {
                gc = queue.poll();
                row = gc.row;
                col = gc.col;
                flatIndex = gc.flatIndex;
                //z = gc.z;
                floodOrder[row][col] = order
                order++
                for (int i = 0; i < 8; i++) {
                    rowN = row + dY[i];
                    colN = col + dX[i];
                    zN = output[rowN][colN];
                    if ((zN != nodata) && (!inQueue.getValue(rowN, colN))) {
                        numSolvedCells++;
                        k = 0
						if (zN == output[row][col]) {
							k = flatIndex + 1
						}
                        //priority = (long)(zN * elevMultiplier * 100000 + (k % 10000))
						queue.add(new GridCell(rowN, colN, zN, k))
						
                        //gc = new GridCell(rowN, colN, zN);
                        //queue.add(gc);
                        inQueue.setValue(rowN, colN, true)
                    }
                }
                progress = (int) (100f * numSolvedCells / numCellsTotal);
                if (progress > oldProgress) {
                    pluginHost.updateProgress("Loop 2 of 3", progress)
                    oldProgress = progress;
                    // check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
                }
            }

            // output the data
			WhiteboxRaster outputRaster = new WhiteboxRaster(outputFile, "rw", 
  		  	     demFile, DataType.FLOAT, nodata)
			outputRaster.setPreferredPalette("spectrum.plt")
  		  	
			oldProgress = -1
  		  	for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					z = floodOrder[row + 1][col + 1]
					if (z != -1) {
						outputRaster.setValue(row, col, z)
					} else {
						outputRaster.setValue(row, col, nodata)
					}
				}
  		  		progress = (int)(100f * row / rowsLessOne)
				if (progress > oldProgress) {
					pluginHost.updateProgress("Loop 3 of 3", progress)
					oldProgress = progress

					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}
			
			outputRaster.addMetadataEntry("Created by the "
	                    + descriptiveName + " tool.")
	        outputRaster.addMetadataEntry("Created on " + new Date())
			outputRaster.close()
	
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

	@CompileStatic
    class GridCell implements Comparable<GridCell> {

        public int row;
        public int col;
        public double z;
        public int flatIndex;
        //public long priority

        public GridCell(int row, int col, double z, int flatIndex) { //, long priority) { // double Z) {
            this.row = row;
            this.col = col;
            this.z = z;
            this.flatIndex = flatIndex;
            //this.priority = priority;
        }

        @Override
        public int compareTo(GridCell other) {
        	if (this.z > other.z) {
        		return 1
        	} else if (this.z < other.z) {
        		return -1
        	} else {
        		if (this.flatIndex > other.flatIndex) {
        			return 1
        		} else if (this.flatIndex < other.flatIndex) {
        			return -1
        		}
				return 0
        	}
        }
    }
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def tdf = new FloodOrder(pluginHost, args, name, descriptiveName)
}
