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
def name = "FastBreachDepressions"
def descriptiveName = "Breach Depressions (Fast)"
def description = "Calculates the distance of grid cells to the nearest downslope stream cell."
def toolboxes = ["DEMPreprocessing"]

public class FastBreachDepressions implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public FastBreachDepressions(WhiteboxPluginHost pluginHost, 
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
            sd.addDialogFile("Output DEM file", "Output DEM File:", "save", "Raster Files (*.dep), DEP", true, false)
			DialogDataInput di = sd.addDialogDataInput("Maximum breach channel length (pixels). Leave blank for none.", "Max. breach length (pixels; optional):", "", true, true)
            DialogCheckBox dcb1 = sd.addDialogCheckBox("Output flow pointer (direction) file", "Output flow pointer file", false)
            DialogCheckBox dcb2 = sd.addDialogCheckBox("Perform flow accumulation after breaching", "Perform flow accumulation", false)
            dcb2.setVisible(false)

            def lstrFile = { evt -> if (evt.getPropertyName().equals("value")) { 
					String value = dcb1.getValue()
            		if (value != null && !value.isEmpty()) { 
            			if (value.toLowerCase().equals("true")) {
            				dcb2.setVisible(true)
            			} else {
            				dcb2.setVisible(false)
            			}
            		}
            	}
            } as PropertyChangeListener
            dcb1.addPropertyChangeListener(lstrFile)
            
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
	  		double SMALL_NUM = 0.0001d
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
			int maxLength = -1
			boolean maxLengthUsed = false
			if (!(args[2].trim()).isEmpty() && !(args[2].toLowerCase().equals("not specified"))) {
				maxLength = Integer.parseInt(args[2])
				maxLengthUsed = true
			}
			String pointerFile = ""
			String flowAccumFile = ""
			boolean outputPointer = false
			if (args.length >= 4 && !(args[3].trim()).isEmpty()) {
				outputPointer = Boolean.parseBoolean(args[3])
				pointerFile = outputFile.replace(".dep", "_flow_pntr.dep")
			}

			boolean performFlowAccumulation = false
			if (args.length >= 5 && !(args[4].trim()).isEmpty() && outputPointer) {
				performFlowAccumulation = Boolean.parseBoolean(args[4])
				flowAccumFile = outputFile.replace(".dep", "_flow_accum.dep")
			}
			
			// read the input image
			WhiteboxRaster dem = new WhiteboxRaster(demFile, "r")
			double nodata = dem.getNoDataValue()
			int rows = dem.getNumberRows()
			int cols = dem.getNumberColumns()
			int rowsLessOne = rows - 1
			int colsLessOne = cols - 1
			int numCellsTotal = rows * cols

			double[][] output = new double[rows + 2][cols + 2]
			BooleanBitArray2D pits = new BooleanBitArray2D(rows + 2, cols + 2)
			BooleanBitArray2D inQueue = new BooleanBitArray2D(rows + 2, cols + 2)
			NibbleArray2D flowdir = new NibbleArray2D(rows + 2, cols + 2)
			PriorityQueue<GridCell> queue = new PriorityQueue<GridCell>((2 * rows + 2 * cols) * 2);

			// find the pit cells and initialize the grids
			numPits = 0
  		  	oldProgress = -1
			for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					z = dem.getValue(row, col)
					output[row + 1][col + 1] = z
					flowdir.setValue(row + 1, col + 1, 0)
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
								queue.add(new GridCell(row + 1, col + 1, z))
								inQueue.setValue(row + 1, col + 1, true)
								flowdir.setValue(row + 1, col + 1, 0)
							} else {
								pits.setValue(row + 1, col + 1, true)
								numPits++ 
								/* raising a pit cell to just lower than the 
								 *  elevation of its lowest neighbour will 
								 *  reduce the length and depth of the trench
								 *  that is necessary to eliminate the pit
								 *  by quite a bit on average.
								 */
								output[row + 1][col + 1] = lowestNeighbour - SMALL_NUM
							}
						}
					} else {
                        numSolvedCells++
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
				flowdir.setValue(row, 0, 0)
				flowdir.setValue(row, cols + 1, 0)
			}
			
			for (col = 0; col < cols + 2; col++) {
				output[0][col] = nodata
				output[rows + 1][col] = nodata
				flowdir.setValue(0, col, 0)
				flowdir.setValue(rows + 1, col, 0)
			}

			// now breach
			oldProgress = (int) (100f * numSolvedCells / numCellsTotal);
            pluginHost.updateProgress("Loop 2 of 3: ", oldProgress);

            if (!maxLengthUsed) {
	            while (queue.isEmpty() == false) {
	                gc = queue.poll();
	                row = gc.row;
	                col = gc.col;
	                z = gc.z;
	                
	                for (int i = 0; i < 8; i++) {
	                    rowN = row + dY[i];
	                    colN = col + dX[i];
	                    zN = output[rowN][colN];
	                    if ((zN != nodata) && (!inQueue.getValue(rowN, colN))) {
	                        flowdir.setValue(rowN, colN, backLink[i])
	                        if (pits.getValue(rowN, colN)) {
	                        	// trace the flowpath back until you find a lower cell
	                        	zTest = zN
	                        	r = rowN
	                        	c = colN
	                        	flag = true
	                        	while (flag) {
	                        		zTest -= SMALL_NUM // ensures a small increment slope
	                        		dir = flowdir.getValue(r, c)
	                        		if (dir > 0) {
	                        			r += dY[dir - 1]
		                            	c += dX[dir - 1]
		                            	zN2 = output[r][c]
		                            	if (zN2 <= zTest || zN2 == nodata) {
		                            		// a lower grid cell has been found
		                            		flag = false
		                            	} else {
		                            		output[r][c] = zTest
		                            		// this cell is already in the 
		                            		// queue but with a higher elevation
		                            		gc = new GridCell(r, c, zTest);
	                        				queue.add(gc);
	                        				numCellsTotal++
		                            	}
	                        		} else {
	                        			flag = false
	                        		}
	                        	}
	                        	pits.setValue(rowN, colN, false)
	                        }
	                        numSolvedCells++;
	                        gc = new GridCell(rowN, colN, zN);
	                        queue.add(gc);
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
            } else {
            	while (queue.isEmpty() == false) {
	                gc = queue.poll();
	                row = gc.row;
	                col = gc.col;
	                z = gc.z;
	                
	                for (int i = 0; i < 8; i++) {
	                    rowN = row + dY[i];
	                    colN = col + dX[i];
	                    zN = output[rowN][colN];
	                    if ((zN != nodata) && (!inQueue.getValue(rowN, colN))) {
	                        flowdir.setValue(rowN, colN, backLink[i])
	                        if (pits.getValue(rowN, colN)) {
	                        	// trace the flowpath back until you find a lower cell
	                        	numCellsInPath = 0
	                        	zTest = zN
	                        	r = rowN
	                        	c = colN
	                        	flag = true
	                        	while (flag) {
	                        		zTest -= SMALL_NUM // ensures a small increment slope
	                        		dir = flowdir.getValue(r, c)
	                        		if (dir > 0) {
	                        			r += dY[dir - 1]
		                            	c += dX[dir - 1]
		                            	zN2 = output[r][c]
		                            	if (zN2 <= zTest || zN2 == nodata) {
		                            		// a lower grid cell has been found
		                            		flag = false
		                            	}
	                        		} else {
	                        			flag = false
	                        		}
	                        		numCellsInPath++
	                        		if (numCellsInPath > maxLength) { flag = false }
	                        	}

								if (numCellsInPath <= maxLength) {
									zTest = zN
		                        	r = rowN
		                        	c = colN
		                        	flag = true
		                        	while (flag) {
		                        		zTest -= SMALL_NUM // ensures a small increment slope
		                        		dir = flowdir.getValue(r, c)
		                        		if (dir > 0) {
		                        			r += dY[dir - 1]
			                            	c += dX[dir - 1]
			                            	zN2 = output[r][c]
			                            	if (zN2 <= zTest || zN2 == nodata) {
			                            		// a lower grid cell has been found
			                            		flag = false
			                            	} else {
			                            		output[r][c] = zTest
			                            		// this cell is already in the 
			                            		// queue but with a higher elevation
			                            		gc = new GridCell(r, c, zTest);
		                        				queue.add(gc);
		                        				numCellsTotal++
			                            	}
		                        		} else {
		                        			flag = false
		                        		}
		                        	}
								}
	                        	pits.setValue(rowN, colN, false)
	                        }
	                        numSolvedCells++;
	                        gc = new GridCell(rowN, colN, zN);
	                        queue.add(gc);
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
            }

            // output the data
			WhiteboxRaster outputRaster = new WhiteboxRaster(outputFile, "rw", 
  		  	     demFile, DataType.FLOAT, nodata)
			outputRaster.setPreferredPalette(paletteName)
  		  	WhiteboxRaster pointer
  		  	if (outputPointer) {
  		  		pointer = new WhiteboxRaster(pointerFile, "rw", 
  		  	     demFile, DataType.FLOAT, nodata)
  		  	    pointer.setDataScale(DataScale.CATEGORICAL)
  		  	    pointer.setPreferredPalette("qual.pal")
  		  	}

			oldProgress = -1
  		  	for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					z = output[row + 1][col + 1]
					outputRaster.setValue(row, col, z)
					if (outputPointer) {
						if (z != nodata) {
							pointer.setValue(row, col, outPointer[flowdir.getValue(row + 1, col + 1)])
						} else {
							pointer.setValue(row, col, nodata)
						}
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
	
			if (outputPointer) {
				pointer.addMetadataEntry("Created by the "
	                    + descriptiveName + " tool.")
		        pointer.addMetadataEntry("Created on " + new Date())
		        pointer.close()
		
				// display the output image
				pluginHost.returnData(pointerFile)
			}

			// display the output image
			pluginHost.returnData(outputFile)
			
			if (performFlowAccumulation) {
				String[] args2 = [pointerFile, flowAccumFile, "number of upslope grid cells", "false"] 
				pluginHost.runPlugin("FlowAccumD8", args2, false)
			}
			
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

        public GridCell(int Row, int Col, double Z) {
            row = Row;
            col = Col;
            z = Z;
        }

        @Override
        public int compareTo(GridCell other) {
            final int BEFORE = -1;
            final int EQUAL = 0;
            final int AFTER = 1;

            if (this.z < other.z) {
                return BEFORE;
            } else if (this.z > other.z) {
                return AFTER;
            }

            if (this.row < other.row) {
                return BEFORE;
            } else if (this.row > other.row) {
                return AFTER;
            }

            if (this.col < other.col) {
                return BEFORE;
            } else if (this.col > other.col) {
                return AFTER;
            }

            return EQUAL;
        }
    }
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def tdf = new FastBreachDepressions(pluginHost, args, name, descriptiveName)
}
