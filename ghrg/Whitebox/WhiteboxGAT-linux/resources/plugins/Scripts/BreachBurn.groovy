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
import whitebox.geospatialfiles.shapefile.ShapeFileRecord
import whitebox.ui.plugin_dialog.*
import whitebox.utilities.StringUtilities
import whitebox.structures.BooleanBitArray2D
import whitebox.structures.NibbleArray2D
import whitebox.structures.DoubleArray2D
import whitebox.structures.IntArray2D
import whitebox.structures.BoundingBox;
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "BreachBurn"
def descriptiveName = "BreachBurn"
//def description = "Burns streams into a DEM."
//def toolboxes = ["DEMPreprocessing"]

public class BreachBurn implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public BreachBurn(WhiteboxPluginHost pluginHost, 
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
            sd.addDialogFile("Input streams file", "Input Streams File:", "open", "Vector Files (*.shp), SHP", true, false)
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
	  		int dir, n
	  		double z, zN, zTest, zN2, lowestNeighbour
	  		boolean isPit, isEdgeCell, isStream, flag
	  		GridCell gc
	  		double LARGE_NUM = Float.MAX_VALUE
			int numInflowing
  		  	double s, sN
  		  	
	  		/*
	  		 * 7  8  1
	  		 * 6  X  2
	  		 * 5  4  3
	  		 */
			int[] dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ]
			int[] dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ]
			int[] backLink = [5, 6, 7, 8, 1, 2, 3, 4]
			double[] outPointer = [0, 1, 2, 4, 8, 16, 32, 64, 128]
			
			if (args.length < 3) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String demFile = args[0]
			String streamsFile = args[1]
			String outputFile = args[2]
			
			// read the input image
			pluginHost.updateProgress("Reading data...", 0);
			WhiteboxRaster dem = new WhiteboxRaster(demFile, "r")
			double nodata = dem.getNoDataValue()
			int rows = dem.getNumberRows()
			int cols = dem.getNumberColumns()
			int rowsLessOne = rows - 1
			int colsLessOne = cols - 1
			int numCellsTotal = rows * cols
			def paletteName = dem.getPreferredPalette()

			//double minVal = dem.getMinimumValue()
			int elevDigits = (String.valueOf((int)(dem.getMaximumValue()))).length()
			double elevMultiplier = Math.pow(10, 8-elevDigits)
			double SMALL_NUM = 1 / elevMultiplier * 10
			long priority

			IntArray2D linkID = new IntArray2D(rows, cols, -32768)
			IntArray2D linkPosition = new IntArray2D(rows, cols, -32768)
			BooleanBitArray2D linkEndNodes = new BooleanBitArray2D(rows, cols)
			
			//BooleanBitArray2D streams = new BooleanBitArray2D(rows, cols)
			
			// perform vector-to-raster conversion
			ShapeFile input = new ShapeFile(streamsFile)
			ShapeType shapeType = input.getShapeType()
            if (shapeType.getBaseType() != ShapeType.POLYLINE) {
            	pluginHost.showFeedback("The input shapefile should be of a POLYLINE ShapeType.")
            	return
            }
            
            int numFeatures = input.getNumberOfRecords()
        	int count = 0
			double[][] points
			int startingPointInPart, endingPointInPart
			int i
			int x1, y1, x2, y2 //xPrime, yPrime;
			BoundingBox box;
			int topRow, bottomRow, leftCol, rightCol;
			double rowYCoord, colXCoord;

			int featureNum = 0;

//			WhiteboxRaster outStream = new WhiteboxRaster(outputFile, "rw", 
//  		  	     demFile, DataType.FLOAT, nodata)
//			//outStream.setPreferredPalette("black_white.plt") //spectrum.plt")
//			outStream.setForceAllDataInMemory(true);

			oldProgress = -1
			for (ShapeFileRecord record : input.records) {
				int recNum = record.getRecordNumber()
                points = record.getGeometry().getPoints()
				int numPoints = points.length;
				int[] partData = record.getGeometry().getParts()
				int numParts = partData.length
				for (int part = 0; part < numParts; part++) {
					featureNum++
					//box = new BoundingBox();             
					startingPointInPart = partData[part];
                    if (part < numParts - 1) {
                        endingPointInPart = partData[part + 1];
                    } else {
                        endingPointInPart = numPoints;
                    }
					n = 0
					for (i = startingPointInPart; i < endingPointInPart - 1; i++) {
						x1 = dem.getColumnFromXCoordinate(points[i][0]);
						y1 = dem.getRowFromYCoordinate(points[i][1])

						x2 = dem.getColumnFromXCoordinate(points[i + 1][0]);
						y2 = dem.getRowFromYCoordinate(points[i + 1][1])

						int d = 0;
 
				        int dy = (int)Math.abs(y2 - y1);
				        int dx = (int)Math.abs(x2 - x1);
				 
				        int dy2 = (dy << 1); // slope scaling factors to avoid floating
				        int dx2 = (dx << 1); // point
				 
				        int ix = x1 < x2 ? 1 : -1; // increment direction
				        int iy = y1 < y2 ? 1 : -1;

				        if (dy <= dx) {
				            for (;;) {
				                if (linkID.getValue(y1, x1) == -32768 && dem.getValue(y1, x1) != nodata) {
				            		n++;
				            		//streams.setValue(y1, x1, true);
                                    linkID.setValue(y1, x1, featureNum);
                                    linkPosition.setValue(y1, x1, n);
                                    //outStream.setValue(y1, x1, featureNum);
				            	}
				                if (x1 == x2) {
				                    break;
				                }
				                x1 += ix;
				                d += dy2;
				                if (d > dx) {
				                    y1 += iy;
				                    d -= dx2;
				                }
				            }
				        } else {
				            for (;;) {
				            	if (linkID.getValue(y1, x1) == -32768 && dem.getValue(y1, x1) != nodata) {
				            		n++;
				            		//streams.setValue(y1, x1, true);
                                    linkID.setValue(y1, x1, featureNum);
                                    linkPosition.setValue(y1, x1, n);
				            		//outStream.setValue(y1, x1, featureNum);
				            	}
				                if (y1 == y2) {
				                    break;
				                }
				                y1 += iy;
				                d += dx2;
				                if (d > dy) {
				                    x1 += ix;
				                    d -= dy2;
				                }
				            }
				        }
					}
				}
				
				count++
                progress = (int)(100f * count / numFeatures)
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Rasterizing Streams...", progress)
            		oldProgress = progress
            		// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}
            }

            //outStream.close()
            //return
            


            /*  Scan the raster looking for stream cells
             *  identifying all grid cells that are either the 
             *  maximum or minimum link position value for 
             *  their link ID within their 3 x 3 neighborhood. 
             *  These are link end nodes (either upstream or 
             *  downstream positions). 
             */
//             WhiteboxRaster outStream = new WhiteboxRaster(outputFile, "rw", 
//  		  	     demFile, DataType.FLOAT, nodata)
//			outStream.setPreferredPalette("black_white.plt") //spectrum.plt")
//			outStream.setForceAllDataInMemory(true);

			int id, idN, position, positionN
			int minNeighbour, maxNeighbour
			oldProgress = -1
			for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					id = linkID.getValue(row, col)
					if (id > 0) { //  it's a stream cell
						position = linkPosition.getValue(row, col)
						minNeighbour = position
						maxNeighbour = position
						for (n = 0; n < 8; n++) {
                    		rowN = row + dY[n];
                    		colN = col + dX[n];
							positionN = linkPosition.getValue(rowN, colN)
							idN = linkID.getValue(rowN, colN)
							if (idN == id) {
								if (positionN < minNeighbour) { minNeighbour = positionN }
								if (positionN > maxNeighbour) { maxNeighbour = positionN }
							}
						}
						if (minNeighbour == position || maxNeighbour == position) {
							linkEndNodes.setValue(row, col, true)
							//outStream.setValue(row, col, 1)
						}
					}
				}
				progress = (int)(100f * row / rowsLessOne)
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

			//outStream.close()
            //return
			
            
			//double[][] output = new double[rows + 2][cols + 2]
			DoubleArray2D output = new DoubleArray2D(rows, cols, nodata)
			NibbleArray2D flowdir = new NibbleArray2D(rows, cols)
			PriorityQueue<GridCell> queue = new PriorityQueue<GridCell>((2 * rows + 2 * cols) * 2);
			BooleanBitArray2D inQueue = new BooleanBitArray2D(rows, cols)
			
			// initialize the grids
			oldProgress = -1
			for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					z = dem.getValue(row, col)
					output.setValue(row, col, z)
					flowdir.setValue(row, col, 0)
					if (z != nodata) {
						isStream = (linkID.getValue(row, col) > 0)
						isPit = true
						isEdgeCell = false
						for (n = 0; n < 8; n++) {
							zN = dem.getValue(row + dY[n], col + dX[n])
							if (isPit && zN != nodata && zN < z) {
								isPit = false
							} else if (zN == nodata) {
								isEdgeCell = true
							}
						}
						if ((isPit && isEdgeCell) || (isStream && isEdgeCell)) {
							queue.add(new GridCell(row, col, z, 0, isStream))
							inQueue.setValue(row, col, true)
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

			oldProgress = (int) (100f * numSolvedCells / numCellsTotal);
            pluginHost.updateProgress("Loop 2 of 3: ", oldProgress);
			int flatIndex
			int k
			int linkIDValue, linkIDValueN
            while (!queue.isEmpty()) {
                gc = queue.poll();
                row = gc.row;
                col = gc.col;
                flatIndex = gc.flatIndex;
                linkIDValue = linkID.getValue(row, col)
                isStream = (linkIDValue > 0)
                if (linkIDValue > 0) {
                	/* move upstream following the path of 
                	 *  minimum change in link position. Stop
                	 *  when the link id changes.
                	 */
                	flag = true
                	r = row
                	c = col
					while (flag) {
						int indexOfNextCell = -1
						int minPosDiff = Integer.MAX_VALUE
						int posDiff
						position = linkPosition.getValue(r, c)
						for (n = 0; n < 8; n++) {
	                    	rowN = r + dY[n];
	                    	colN = c + dX[n];
	                    	linkIDValueN = linkID.getValue(rowN, colN)
	                    	if ((linkIDValueN == linkIDValue) && (!inQueue.getValue(rowN, colN))) {
	                    		positionN = linkPosition.getValue(rowN, colN)
	                    		posDiff = (positionN - position) * (positionN - position)
	                    		if (posDiff < minPosDiff) {
	                    			minPosDiff = posDiff
	                    			indexOfNextCell = n
	                    		}
	                    	} else if (linkIDValueN == -32768 || linkEndNodes.getValue(rowN, colN)) {
	                    		zN = output.getValue(rowN, colN)
	                    		if ((zN != nodata) && (!inQueue.getValue(rowN, colN))) {
			                    	// it's a non-stream cell or a link end node and can be added to the queue
	                    			numSolvedCells++;
			                        flowdir.setValue(rowN, colN, backLink[n])
			                        k = 0
									if (zN == output.getValue(row, col)) {
										k = flatIndex + 1
									}
			                        queue.add(new GridCell(rowN, colN, zN, k, linkEndNodes.getValue(rowN, colN)))
									
			                        inQueue.setValue(rowN, colN, true)
			                    }
	                    	}
						}	
						if (indexOfNextCell > -1) {
							rowN = r + dY[indexOfNextCell];
	                    	colN = c + dX[indexOfNextCell];
							numSolvedCells++;
			                flowdir.setValue(rowN, colN, backLink[indexOfNextCell])
			                inQueue.setValue(rowN, colN, true)
	                    	r = rowN;
	                    	c = colN;
						} else {
							// there were no unvisited neighbours of the same link ID
							flag = false
						}
					}
                } else {
	                for (n = 0; n < 8; n++) {
	                    rowN = row + dY[n];
	                    colN = col + dX[n];
	                    zN = output.getValue(rowN, colN)
	                    if ((zN != nodata) && (!inQueue.getValue(rowN, colN))) {
	                    	linkIDValueN = linkID.getValue(rowN, colN)
	                    	if (linkIDValueN == -32768 || linkEndNodes.getValue(rowN, colN)) {
		                        numSolvedCells++;
		                        flowdir.setValue(rowN, colN, backLink[n])
		                        k = 0
								if (zN == output.getValue(row, col)) {
									k = flatIndex + 1
								}
		                        queue.add(new GridCell(rowN, colN, zN, k, linkEndNodes.getValue(rowN, colN)))
								
		                        inQueue.setValue(rowN, colN, true)
	                    	} // else it's a stream and not an end node and shouldn't be added to the queue.
	                    }
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
					z = output.getValue(row, col)
					if (z != nodata) {
						outputRaster.setValue(row, col, outPointer[flowdir.getValue(row, col)])
					}
				}
				progress = (int)(100f * row / rowsLessOne)
				if (progress > oldProgress) {
					pluginHost.updateProgress("Outputting data...", progress)
					oldProgress = progress

					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}
			
//
//			inQueue = new BooleanBitArray2D(rows, cols)
//  		  	oldProgress = -1
//  		  	for (row = 0; row < rows; row++) {
//				for (col = 0; col < cols; col++) {
//					if (streams.getValue(row, col)) {
//						numInflowing = 0			
//		                for (n = 0; n < 8; n++) {
//		                    rowN = row + dY[n];
//		                    colN = col + dX[n];
//		                    if (streams.getValue(rowN, colN) && 
//		                        flowdir.getValue(rowN+1, colN+1) == backLink[n]) { 
//		                        	numInflowing++ 
//		                    }
//		                }
//		                if (numInflowing == 0) {
//		                	r = row
//							c = col
//							zTest = dem.getValue(r, c) - SMALL_NUM
//							flag = true
//							while (flag) {
//								if (zTest > output[r + 1][c + 1] && inQueue.getValue(r, c)) {
//									flag = false // we've already traversed the stream
//								} else {
//									inQueue.setValue(r, c, true)
//									dir = flowdir.getValue(r + 1, c + 1) - 1
//									if (dir >= 0) {
//										// now find the lowest neighbour that isn't the downstream cell
//										lowestNeighbour = LARGE_NUM
//										for (n = 0; n < 8; n++) {
//											if (n != dir) { // we don't expect it to be lower than the cell it flows to
//												zN = output[r + dY[n] + 1][c + dX[n] + 1]
//							                    if (zN != nodata && zN < lowestNeighbour) { 
//							                        lowestNeighbour = zN 
//							                    }
//											}
//						                }
//
//						                zTest = lowestNeighbour - SMALL_NUM
//						                if (zTest < output[r + 1][c + 1]) {
//						                	output[r + 1][c + 1] = zTest
//						                } else {
//											zTest = output[r + 1][c + 1] - SMALL_NUM
//										}
//					                
//										r += dY[dir]
//										c += dX[dir]
//									} else {
//										lowestNeighbour = LARGE_NUM
//										for (n = 0; n < 8; n++) {
//											zN = output[r + dY[n] + 1][c + dX[n] + 1]
//							                if (zN != nodata && zN < lowestNeighbour) { 
//							                    lowestNeighbour = zN 
//							                }
//						                }
//
//						                zTest = lowestNeighbour - SMALL_NUM
//						                if (zTest < output[r + 1][c + 1]) {
//						                	output[r + 1][c + 1] = zTest
//						                }
//
//										flag = false
//									}
//								}
//							}
//							
//		                }
//					}// else {
//					//	outputRaster.setValue(row, col, dem.getValue(row, col))
//					//}
//				}
//  		  		progress = (int)(100f * row / rowsLessOne)
//				if (progress > oldProgress) {
//					pluginHost.updateProgress("Loop 3 of 4", progress)
//					oldProgress = progress
//
//					// check to see if the user has requested a cancellation
//					if (pluginHost.isRequestForOperationCancelSet()) {
//						pluginHost.showFeedback("Operation cancelled")
//						return
//					}
//				}
//			}
//
//			
//
//			oldProgress = -1
//  		  	for (row = 0; row < rows; row++) {
//				for (col = 0; col < cols; col++) {
//					z = output[row + 1][col + 1]
//					outputRaster.setValue(row, col, z)
//					if (streams.getValue(row, col)) {
//						outputRaster.setValue(row, col, 1.0)
//					} else {
//						outputRaster.setValue(row, col, 0.0)
//					}
//				}
//  		  		progress = (int)(100f * row / rowsLessOne)
//				if (progress > oldProgress) {
//					pluginHost.updateProgress("Loop 4 of 4", progress)
//					oldProgress = progress
//
//					// check to see if the user has requested a cancellation
//					if (pluginHost.isRequestForOperationCancelSet()) {
//						pluginHost.showFeedback("Operation cancelled")
//						return
//					}
//				}
//			}

			dem.close()

			outputRaster.setPreferredPalette(paletteName)
			//outputRaster.setDisplayMinimum(dem.getDisplayMinimum())
			//outputRaster.setDisplayMaximum(dem.getDisplayMaximum())
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
        public boolean streamVal;
        //public long priority

        public GridCell(int row, int col, double z, int flatIndex, boolean streamVal) { //, long priority) { // double Z) {
            this.row = row;
            this.col = col;
            this.z = z;
            this.flatIndex = flatIndex;
            this.streamVal = streamVal;
        }

        @Override
        public int compareTo(GridCell other) {
        	if (this.streamVal && !other.streamVal) {
        		return -1
        	} else if (!this.streamVal && other.streamVal) {
        		return 1
        	} else {
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

    // Return true if val is between theshold1 and theshold2.
    @CompileStatic
    private static boolean isBetween(double val, double threshold1, double threshold2) {
        if (val == threshold1 || val == threshold2) {
            return true;
        }
        return threshold2 > threshold1 ? val > threshold1 && val < threshold2 : val > threshold2 && val < threshold1;
    }
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def tdf = new BreachBurn(pluginHost, args, name, descriptiveName)
}
