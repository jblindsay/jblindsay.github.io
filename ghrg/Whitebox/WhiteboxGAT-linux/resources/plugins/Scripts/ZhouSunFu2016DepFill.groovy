/*
 * Copyright (C) 2016 Dr. John Lindsay <jlindsay@uoguelph.ca>
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

 /* This is an experiment with the Zhou, Sun, and Fu (2016) modified
  *  priority flood depression filling algorithm. The algorithm uses
  *  a regular queue for most of the heavy lifting instead of a 
  *  priority queue. According to the publication, the efficiency
  *  difference between a queue and a pq is so great that the
  *  algorithm fills DEMs in about half the time as the regular
  *  prority flood method. I was never able to acheive anything
  *  near that level of improvement. In fact, for many DEMs the
  *  modified algorithm was slower. I don't know if it is down to
  *  the queue implementation in Java being inefficient or the
  *  pq implementation being partiuclarly efficient but either way
  *  the gains in efficiency just didn't pan out.
  */
 
import java.awt.event.ActionListener
import java.awt.event.ActionEvent
import java.beans.PropertyChangeEvent
import java.beans.PropertyChangeListener
import java.util.Date
import java.util.ArrayList
import java.util.LinkedList
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
import groovy.transform.CompileStatic
import groovy.time.TimeDuration
import groovy.time.TimeCategory

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
//def name = "TestDepFill"
//def descriptiveName = "TestDepFill"
//def description = "Calculates the distance of grid cells to the nearest downslope stream cell."
//def toolboxes = ["DEMPreprocessing"]

public class TestDepFill implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public TestDepFill(WhiteboxPluginHost pluginHost, 
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
			
			// resize the dialog to the standard size and display it
			sd.setSize(800, 400)
			sd.visible = true
		}
	}

	@CompileStatic
	private void execute(String[] args) {
		try {
			int progress, oldProgress, col, row, colN, rowN, r, c
	  		int numSolvedCells = 0
	  		double zIn, zNIn, zOut, zNOut, z, z2;
	  		double SMALL_NUM = 0.00; //0.001;
	  		boolean addToPQ
	  		GridCell gc
			int[] dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ]
			int[] dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ]
			
			if (args.length < 2) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String demFile = args[0]
			String outputFile = args[1]
			
			// read the input image
			WhiteboxRaster dem = new WhiteboxRaster(demFile, "r")
			dem.setForceAllDataInMemory(true);
			double nodata = dem.getNoDataValue()
			int rows = dem.getNumberRows()
			int cols = dem.getNumberColumns()
			int rowsLessOne = rows - 1
			int colsLessOne = cols - 1
			int numCellsTotal = rows * cols

			double nodataOut = -32768.0;
			double notInQueueVal = -998.0;
  		  	WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", demFile, DataType.FLOAT, notInQueueVal)
  		  	//WhiteboxRaster output2 = new WhiteboxRaster(outputFile.replace(".dep", "_tmp1.dep"), "rw", demFile, DataType.FLOAT, nodataOut)
  		  	output.setForceAllDataInMemory(true);
  		  	output.setNoDataValue(nodataOut);
			def paletteName = dem.getPreferredPalette()
			output.setPreferredPalette(paletteName)

			Date start = new Date()
  		  	
			PriorityQueue<GridCell> pq = new PriorityQueue<GridCell>((2 * rows + 2 * cols) * 2);
			//Deque<GridCell> q = new LinkedList<GridCell>();
				
			pluginHost.updateProgress("Loop 1 of 3", 0);
			
//			for (row = 0; row < rows; row++) {
//				col = 0;
//				zIn = dem.getValue(row, col);
//				if (zIn == nodata) { zIn = nodataOut; }
//				output.setValue(row, col, zIn);
//				q.addLast(new GridCell(row, col, zIn));
//				col = cols - 1;
//				zIn = dem.getValue(row, col);
//				if (zIn == nodata) { zIn = nodataOut; }
//				output.setValue(row, col, zIn);
//				q.addLast(new GridCell(row, col, zIn));
//				numSolvedCells += 2;
//			}
//			
//			for (col = 1; col < cols - 1; col++) {
//				row = 0;
//				zIn = dem.getValue(row, col);
//				if (zIn == nodata) { zIn = nodataOut; }
//				output.setValue(row, col, zIn);
//				q.addLast(new GridCell(row, col, zIn));
//				row = rows - 1;
//				zIn = dem.getValue(row, col);
//				if (zIn == nodata) { zIn = nodataOut; }
//				output.setValue(row, col, zIn);
//				q.addLast(new GridCell(row, col, zIn));
//				numSolvedCells += 2;
//			}
//
//			while (!q.isEmpty()) {
//            	gc = q.poll();
//				r = gc.row;
//				c = gc.col;
//				z2 = gc.z;
//				for (int j = 0; j < 8; j++) {
//					rowN = r + dY[j];
//    				colN = c + dX[j];
//    				if (output.getValue(rowN, colN) == notInQueueVal) {
//    					zNOut = dem.getValue(rowN, colN);
//	                    if (zNOut == nodata) { zNOut = nodataOut; }
//	                    if (zNOut <= z2) {
//	                    	z3 = z2 + SMALL_NUM;
//    						if (z2 == nodataOut) { z3 = nodataOut; }
//	                    	output.setValue(rowN, colN, z3); //z2 + SMALL_NUM);
//	                    	output2.setValue(rowN, colN, 1);
//           					order++;
//            
//							q.addLast(new GridCell(rowN, colN, z3)); //z2 + SMALL_NUM));
//							numSolvedCells++;
//	                    } else {
//	                    	// There's a higher neighbour of r, c
//	                    	addToPQ = true;
//	                    }
//    				}
//				}
//				if (addToPQ) {
//					q2.addLast(new GridCell(r, c, z2));
//				}
//            }

			


			for (row = 0; row < rows; row++) {
				col = 0;
				zIn = dem.getValue(row, col);
				if (zIn == nodata) { zIn = nodataOut; }
				output.setValue(row, col, zIn);
				pq.add(new GridCell(row, col, zIn));
				col = cols - 1;
				zIn = dem.getValue(row, col);
				if (zIn == nodata) { zIn = nodataOut; }
				output.setValue(row, col, zIn);
				pq.add(new GridCell(row, col, zIn));
				numSolvedCells += 2;
			}
			
			for (col = 1; col < cols - 1; col++) {
				row = 0;
				zIn = dem.getValue(row, col);
				if (zIn == nodata) { zIn = nodataOut; }
				output.setValue(row, col, zIn);
				pq.add(new GridCell(row, col, zIn));
				row = rows - 1;
				zIn = dem.getValue(row, col);
				if (zIn == nodata) { zIn = nodataOut; }
				output.setValue(row, col, zIn);
				pq.add(new GridCell(row, col, zIn));
				numSolvedCells += 2;
			}

//			int k = 0;
//			int m = 0;
			boolean regionGrow = false;
			if (regionGrow) {
				double z3;
				//Deque<GridCell> q = new ArrayDeque<GridCell>((2 * rows + 2 * cols) * 2);
				//Deque<GridCell> q2 = new ArrayDeque<GridCell>((2 * rows + 2 * cols) * 2);
				Deque<GridCell> q = new LinkedList<GridCell>();
				Deque<GridCell> q2 = new LinkedList<GridCell>();

				double order = 0;
				while (!pq.isEmpty() && numSolvedCells < numCellsTotal) {
					gc = pq.poll();
	                row = gc.row;
	                col = gc.col;
	                zOut = gc.z;
//	                k++;
	                for (int i = 0; i < 8; i++) {
	                    rowN = row + dY[i];
	                    colN = col + dX[i];
	                    if (output.getValue(rowN, colN) == notInQueueVal) {
	                    	zNOut = dem.getValue(rowN, colN);
	                    	if (zNOut == nodata) { zNOut = nodataOut; }
	                    	if (zNOut <= zOut) {
	                    		z = zOut + SMALL_NUM;
	                    		if (zOut == nodataOut) { z = nodataOut; }
	                  
	                            output.setValue(rowN, colN, z);
	                            //output2.setValue(rowN, colN, 1);
	                            order++;
	                            numSolvedCells++;
	                            // Find all depression cells
	                            q.addLast(new GridCell(rowN, colN, z));
	                            while (!q.isEmpty()) {
	                            	gc = q.poll();
	                				r = gc.row;
	                				c = gc.col;
	                				z2 = gc.z;
	                				addToPQ = false
	                				for (int j = 0; j < 8; j++) {
	                					rowN = r + dY[j];
	                    				colN = c + dX[j];
	                    				if (output.getValue(rowN, colN) == notInQueueVal) {
	                    					zNOut = dem.getValue(rowN, colN);
						                    if (zNOut == nodata) { zNOut = nodataOut; }
						                    if (zNOut <= z2) {
						                    	z3 = z2 + SMALL_NUM;
	                    						if (z2 == nodataOut) { z3 = nodataOut; }
						                    	output.setValue(rowN, colN, z3); //z2 + SMALL_NUM);
						                    	//output2.setValue(rowN, colN, 1);
	                           					order++;
	                            
												q.addLast(new GridCell(rowN, colN, z3)); //z2 + SMALL_NUM));
												numSolvedCells++;
						                    } else {
						                    	// There's a higher neighbour of r, c
						                    	addToPQ = true;
						                    }
	                    				}
	                				}
	                				if (addToPQ) {
//	                					pq.add(new GridCell(r, c, z2));
	                					q2.addLast(new GridCell(r, c, z2));
	                				}
	                            }
	                            while (!q2.isEmpty()) {
	                            	gc = q2.poll();
	                				r = gc.row;
	                				c = gc.col;
	                				z2 = gc.z;
	                            	addToPQ = false;
	                            	for (int j = 0; j < 8; j++) {
	                            		rowN = r + dY[j];
	                    				colN = c + dX[j];
	                    				if (output.getValue(rowN, colN) == notInQueueVal) {
	                    					addToPQ = true;
	                    					break;
	                    				}
	                            	}
	                            	if (addToPQ) { pq.add(new GridCell(r, c, z2)); }
	                            }
	                        } else {
	                			output.setValue(rowN, colN, zNOut);
	                			//output2.setValue(rowN, colN, 2);
	                            numSolvedCells++;
	                            // Find all slope cells
	                            q.addLast(new GridCell(rowN, colN, zNOut));
	                            while (!q.isEmpty()) {
	                            	gc = q.poll();
	                				r = gc.row;
	                				c = gc.col;
	                				z2 = gc.z;
	                				addToPQ = false;
	                				for (int j = 0; j < 8; j++) {
	                					rowN = r + dY[j];
	                    				colN = c + dX[j];
	                    				if (output.getValue(rowN, colN) == notInQueueVal) {
	                    					zNOut = dem.getValue(rowN, colN);
						                    if (zNOut == nodata) { zNOut = nodataOut; }
						                    if (zNOut > z2) {
						                    	output.setValue(rowN, colN, zNOut);
						                    	//output2.setValue(rowN, colN, 2);
	                            				order++;
	                            
												q.addLast(new GridCell(rowN, colN, zNOut));
												numSolvedCells++;
						                    } else {
						                    	addToPQ = true;
						                    }
	                    				}
	                				}
	                				if (addToPQ) {
//	                					pq.add(new GridCell(r, c, z2));
	                					q2.addLast(new GridCell(r, c, z2));
	                				}
	                            }
	                            while (!q2.isEmpty()) {
	                            	gc = q2.poll();
	                				r = gc.row;
	                				c = gc.col;
	                				z2 = gc.z;
	                            	addToPQ = false;
	                            	for (int j = 0; j < 8; j++) {
	                            		rowN = r + dY[j];
	                    				colN = c + dX[j];
	                    				if (output.getValue(rowN, colN) == notInQueueVal) {
	                    					addToPQ = true;
	                    					break;
	                    				}
	                            	}
	                            	if (addToPQ) { pq.add(new GridCell(r, c, z2)); }
	                            }
	                        }
	                    }
	                }
					progress = (int) (100f * numSolvedCells / numCellsTotal);
	                if (progress != oldProgress) {
	                    pluginHost.updateProgress("Loop 2 of 2", progress)
	                    oldProgress = progress;
	                    // check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
	                }
				}
			} else {
				while (!pq.isEmpty() && numSolvedCells < numCellsTotal) {
					gc = pq.poll();
	                row = gc.row;
	                col = gc.col;
	                zOut = gc.z;
//	                k++;
	                for (int i = 0; i < 8; i++) {
	                    rowN = row + dY[i];
	                    colN = col + dX[i];
	                    if (output.getValue(rowN, colN) == notInQueueVal) {
	                    	zNOut = dem.getValue(rowN, colN);
	                    	if (zNOut == nodata) { zNOut = nodataOut; }
	                    	if (zNOut <= zOut) {
	                            zNOut = zOut + SMALL_NUM;
	                        }
	                        output.setValue(rowN, colN, zNOut);
                            pq.add(new GridCell(rowN, colN, zNOut));
                            numSolvedCells++;
	                    }
	                }
					progress = (int) (100f * numSolvedCells / numCellsTotal);
	                if (progress != oldProgress) {
	                    pluginHost.updateProgress("Loop 2 of 2", progress)
	                    oldProgress = progress;
	                    // check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
	                }
				}
			}
			
			Date stop = new Date();
			TimeDuration td = TimeCategory.minus(stop, start);
			
			// output the data
            dem.close();
			output.addMetadataEntry("Created by the " + descriptiveName + " tool.")
	        output.addMetadataEntry("Created on " + new Date())
			output.close()
			//output2.close()

			pluginHost.showFeedback("Elapsed time: $td");
//			pluginHost.showFeedback("k = $k")
			
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

        public GridCell(int row, int col, double z) {
            this.row = row;
            this.col = col;
            this.z = z;
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
	def tdf = new TestDepFill(pluginHost, args, name, descriptiveName)
}
