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

 /* 
  *  This tool has been created to compare the results of the improved
  *  TopologicalBreachBurn method to the more commonly applied 
  *  FillBurn method of Saunders, 1999. This tool is not recommended
  *  for widespread use as TopologicalBreachBurn will almost always
  *  provide a better result. PLEASE DO NOT USE IT.
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
import whitebox.structures.BoundingBox;
import groovy.transform.CompileStatic
import groovy.time.TimeDuration
import groovy.time.TimeCategory

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "FillBurn"
def descriptiveName = "FillBurn"
def description = "Burns streams into a DEM using the FillBurn (Saunders, 1999) method."
def toolboxes = ["DEMPreprocessing"]

public class FillBurn implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public FillBurn(WhiteboxPluginHost pluginHost, 
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
	  		Date start = new Date()
			
			int progress, oldProgress, col, row, colN, rowN, numPits, r, c
	  		int numSolvedCells = 0;
	  		int numParts, numPoints, recNum, part, p, i, a, numFeatures
	  		double[][] points;
	  		int[] partData;
	  		BoundingBox box;
	  		int startingPointInPart, endingPointInPart
	  		int topRow, bottomRow, leftCol, rightCol;
			double rowYCoord, colXCoord;
	  		double x1, y1, x2, y2, xPrime, yPrime, d;
			int dir, numCellsInPath;
	  		double z, zN, zTest, zN2, lowestNeighbour;
	  		boolean isPit, isEdgeCell, flag;
	  		double SMALL_NUM = 0.0001d;
	  		GridCell gc;
			int[] dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ];
			int[] dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ];
			int[] backLink = [5, 6, 7, 8, 1, 2, 3, 4];
			double[] outPointer = [0, 1, 2, 4, 8, 16, 32, 64, 128];
			
			if (args.length < 2) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String demFile = args[0]
			String streamsFile = args[1]
			String outputFile = args[2]
			
			
			// read the DEM
			WhiteboxRaster dem = new WhiteboxRaster(demFile, "r")
			dem.setForceAllDataInMemory(true);
			double nodata = dem.getNoDataValue()
			int rows = dem.getNumberRows()
			int cols = dem.getNumberColumns()
			int rowsLessOne = rows - 1
			int colsLessOne = cols - 1
			int numCellsTotal = rows * cols
			String paletteName = dem.getPreferredPalette();

			// read the streams
			ShapeFile input = new ShapeFile(streamsFile)
			ShapeType shapeType = input.getShapeType()
            if (shapeType != ShapeType.POLYLINE) {
            	pluginHost.showFeedback("The input shapefile should be of a POLYLINE ShapeType.")
            	return
            }
            
			// Rasterize the streams
			BooleanBitArray2D streams = new BooleanBitArray2D(rows, cols)
			int featureNum = 0;
			int count = 0;
			oldProgress = -1
			for (ShapeFileRecord record : input.records) {
				recNum = record.getRecordNumber()
                points = record.getGeometry().getPoints()
				numPoints = points.length;
				partData = record.getGeometry().getParts()
				numParts = partData.length
				for (part = 0; part < numParts; part++) {
					featureNum++
					box = new BoundingBox();             
					startingPointInPart = partData[part];
                    if (part < numParts - 1) {
                        endingPointInPart = partData[part + 1];
                    } else {
                        endingPointInPart = numPoints;
                    }

					row = dem.getRowFromYCoordinate(points[startingPointInPart][1])
					col = dem.getColumnFromXCoordinate(points[startingPointInPart][0]);
                    streams.setValue(row, col, true);
                    
                    row = dem.getRowFromYCoordinate(points[endingPointInPart-1][1])
					col = dem.getColumnFromXCoordinate(points[endingPointInPart-1][0]);
                    streams.setValue(row, col, true);
                    
                    for (i = startingPointInPart; i < endingPointInPart; i++) {
                        if (points[i][0] < box.getMinX()) {
                            box.setMinX(points[i][0]);
                        }
                        if (points[i][0] > box.getMaxX()) {
                            box.setMaxX(points[i][0]);
                        }
                        if (points[i][1] < box.getMinY()) {
                            box.setMinY(points[i][1]);
                        }
                        if (points[i][1] > box.getMaxY()) {
                            box.setMaxY(points[i][1]);
                        }
                    }
                    topRow = dem.getRowFromYCoordinate(box.getMaxY());
                    bottomRow = dem.getRowFromYCoordinate(box.getMinY());
                    leftCol = dem.getColumnFromXCoordinate(box.getMinX());
                    rightCol = dem.getColumnFromXCoordinate(box.getMaxX());

					// find each intersection with a row.
                    for (row = topRow; row <= bottomRow; row++) {

                        rowYCoord = dem.getYCoordinateFromRow(row);
                        // find the x-coordinates of each of the line segments 
                        // that intersect this row's y coordinate

                        for (i = startingPointInPart; i < endingPointInPart - 1; i++) {
                            if (isBetween(rowYCoord, points[i][1], points[i + 1][1])) {
                                y1 = points[i][1];
                                y2 = points[i + 1][1];
                                if (y2 != y1) {
                                    x1 = points[i][0];
                                    x2 = points[i + 1][0];

                                    // calculate the intersection point
                                    xPrime = x1 + (rowYCoord - y1) / (y2 - y1) * (x2 - x1);
                                    col = dem.getColumnFromXCoordinate(xPrime);
                                    
									streams.setValue(row, col, true);
                                }
                            }
                        }
                    }

                    // find each intersection with a column.
                    for (col = leftCol; col <= rightCol; col++) {
                        colXCoord = dem.getXCoordinateFromColumn(col);
                        for (i = startingPointInPart; i < endingPointInPart - 1; i++) {
                            if (isBetween(colXCoord, points[i][0], points[i + 1][0])) {
                                x1 = points[i][0];
                                x2 = points[i + 1][0];
                                if (x1 != x2) {
                                    y1 = points[i][1];
                                    y2 = points[i + 1][1];

                                    // calculate the intersection point
                                    yPrime = y1 + (colXCoord - x1) / (x2 - x1) * (y2 - y1);
                                    row = dem.getRowFromYCoordinate(yPrime);

                                    streams.setValue(row, col, true);
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

            
			// Perform line-thinning
            long counter = 1;
	        int loopNum = 0;
	        int[][] elements = [ [6, 7, 0, 4, 3, 2], [7, 0, 1, 3, 5], 
	            [0, 1, 2, 4, 5, 6], [1, 2, 3, 5, 7], 
	            [2, 3, 4, 6, 7, 0], [3, 4, 5, 7, 1], 
	            [4, 5, 6, 0, 1, 2], [5, 6, 7, 1, 3] ];
	        boolean[][] vals = [ [false, false, false, true, true, true], [false, false, false, true, true], 
	            [false, false, false, true, true, true], [false, false, false, true, true],
	            [false, false, false, true, true, true], [false, false, false, true, true],
	            [false, false, false, true, true, true], [false, false, false, true, true] ];
	        
	        boolean[] neighbours = new boolean[8];
	        boolean patternMatch = false;

			while (counter > 0) {
                loopNum++;
                pluginHost.updateProgress("Line Thinning (Loop Number " + loopNum + "):", 0);
                counter = 0;
                for (row = 0; row < rows; row++) {
                    for (col = 0; col < cols; col++) {
                        if (streams.getValue(row, col)) {
                            // fill the neighbours array
                            for (i = 0; i < 8; i++) {
                                neighbours[i] = streams.getValue(row + dY[i], col + dX[i]);
                            }
                            
                            for (a = 0; a < 8; a++) {
                                // scan through element
                                patternMatch = true;
                                for (i = 0; i < elements[a].length; i++) {
                                    if (neighbours[elements[a][i]] != vals[a][i]) {
                                        patternMatch = false;
                                        //break;
                                    }
                                }
                                if (patternMatch) {
                                    streams.setValue(row, col, false);
                                    counter++;
                                }
                            }
                        }

                    }

                    progress = (int)(100f * row / rowsLessOne)
	            	if (progress != oldProgress) {
						pluginHost.updateProgress("Line Thinning (Loop Number " + loopNum + "):", progress)
	            		oldProgress = progress
	            		// check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
	            	}
                }
            }

			// Make a copy of the DEM where each stream cell
			//  has been lowered by 10,000 elevation units.
			WhiteboxRaster outputRaster = new WhiteboxRaster(outputFile, "rw", 
  		  	     demFile, DataType.FLOAT, nodata)
			outputRaster.setPreferredPalette(paletteName)
  		  	outputRaster.setForceAllDataInMemory(true);
  		  	outputRaster.setDisplayMinimum(dem.getDisplayMinimum());
			outputRaster.setDisplayMaximum(dem.getDisplayMaximum());
			
			oldProgress = -1
  		  	for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					z = dem.getValue(row, col)
					if (!streams.getValue(row, col)) {
						outputRaster.setValue(row, col, z)
					} else { // it's a stream
						if (z != nodata) {
							outputRaster.setValue(row, col, z - 1000.0)
						} else {
							outputRaster.setValue(row, col, nodata)
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

			// Fill the streams-decremented DEM

			// Find the minimum elevaton difference between the
			// filled DEM and the original DEM along the 
			// stream network and raise all stream cells by this
			// value less 1 m.

			//double[][] output = new double[rows + 2][cols + 2]
			BooleanBitArray2D inQueue = new BooleanBitArray2D(rows, cols)
			NibbleArray2D flowdir = new NibbleArray2D(rows, cols)
			PriorityQueue<GridCell> queue = new PriorityQueue<GridCell>((2 * rows + 2 * cols) * 2);

			// find the pit cells and initialize the grids
			numPits = 0
  		  	oldProgress = -1
			for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					z = outputRaster.getValue(row, col)
					flowdir.setValue(row + 1, col + 1, 0)
					if (z != nodata) {
						isEdgeCell = false
						for (int n = 0; n < 8; n++) {
							zN = outputRaster.getValue(row + dY[n], col + dX[n])
							if (zN == nodata) {
								isEdgeCell = true;
								break;
							}
						}
						if (isEdgeCell) {
							queue.add(new GridCell(row, col, z))
							inQueue.setValue(row, col, true)
							flowdir.setValue(row, col, 0)
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

			byte[][] numInflow = new byte[rows][cols];
            while (!queue.isEmpty()) {
                gc = queue.poll();
                row = gc.row;
                col = gc.col;
                z = gc.z;
                
                for (i = 0; i < 8; i++) {
                    rowN = row + dY[i];
                    colN = col + dX[i];
                    zN = outputRaster.getValue(rowN, colN);
                    if ((zN != nodata) && (!inQueue.getValue(rowN, colN))) {
                        flowdir.setValue(rowN, colN, backLink[i]);
                        if (zN <= z) {
                        	zN = z + 0.01;
                        	outputRaster.setValue(rowN, colN, zN);
                        }
                        numSolvedCells++;
                        gc = new GridCell(rowN, colN, zN);
                        queue.add(gc);
                        inQueue.setValue(rowN, colN, true)
                        numInflow[row][col]++;
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

            // Find the min diff in-stream value between the burned-filled DEM and original DEM
            double minDiff = Double.MAX_VALUE;
            oldProgress = -1
  		  	for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					if (streams.getValue(row, col) && dem.getValue(row, col) != nodata) {
						z = dem.getValue(row, col) - outputRaster.getValue(row, col);
						if (z < minDiff) { minDiff = z; }
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

			minDiff += 1.0
			oldProgress = -1
  		  	for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					if (streams.getValue(row, col) && dem.getValue(row, col) != nodata) {
						z = outputRaster.getValue(row, col) + minDiff;
						outputRaster.setValue(row, col, z);
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


			// Create the D8 flow pointer raster
			// output the flow direction grid
			WhiteboxRaster outputRasterFD = new WhiteboxRaster(outputFile.replace(".dep", "_D8dir.dep"), "rw", 
  		  	     demFile, DataType.FLOAT, nodata)
			outputRasterFD.setPreferredPalette("spectrum.plt")

			oldProgress = -1
			for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					z = dem.getValue(row, col)
					if (z != nodata) {
						outputRasterFD.setValue(row, col, outPointer[flowdir.getValue(row, col)])
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

			outputRasterFD.addMetadataEntry("Created by the "
	                    + descriptiveName + " tool.")
	        outputRasterFD.addMetadataEntry("Created on " + new Date())
			Date stop = new Date()
			TimeDuration td = TimeCategory.minus(stop, start)
			outputRasterFD.addMetadataEntry("Elapsed time: $td")
			outputRasterFD.close()

			pluginHost.returnData(outputFile.replace(".dep", "_D8dir.dep"))
	
			// Create the D8 flow accumulation raster
			WhiteboxRaster outputFA = new WhiteboxRaster(outputFile.replace(".dep", "_D8FA.dep"), "rw", 
  		  	     demFile, DataType.FLOAT, 1d)
			outputFA.setPreferredPalette("grey.pal");
            outputFA.setDataScale(DataScale.CONTINUOUS);
            outputFA.setZUnits("dimensionless");
            outputFA.setNonlinearity(0.2)
			int rN, cN;
            for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					if (numInflow[row][col] == 0) {
						r = row
						c = col
						numInflow[r][c] = -1
						flag = true
						while (flag) {
							dir = flowdir.getValue(r, c) - 1; //flowDir[r][c]
							if (dir >= 0) {
								rN = r + dY[dir]
								cN = c + dX[dir]
								if (flowdir.getValue(rN, cN) >= 0) {
									z = outputFA.getValue(r, c)
									outputFA.incrementValue(rN, cN, z)
									numInflow[rN][cN] -= 1
									if (numInflow[rN][cN] > 0) {
										flag = false
									} else {
										numInflow[rN][cN] = -1
										r = rN
										c = cN
									}
								} else if (flowdir.getValue(rN, cN) == -1) {
									flag = false
									z = outputFA.getValue(r, c)
									outputFA.incrementValue(rN, cN, z)
								}
							} else {
								flag = false
							}
						}
					}
					if (dem.getValue(row, col) == nodata) {
						outputFA.setValue(row, col, nodata)
					}
				}
				progress = (int)(100f * row / rowsLessOne)
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

			outputFA.addMetadataEntry("Created by the "
	                    + descriptiveName + " tool.")
	        outputFA.addMetadataEntry("Created on " + new Date())
			stop = new Date()
			td = TimeCategory.minus(stop, start)
			outputFA.addMetadataEntry("Elapsed time: $td")
			outputFA.close()

			pluginHost.returnData(outputFile.replace(".dep", "_D8FA.dep"))


			// outputs and final ops
			dem.close();
			outputRaster.addMetadataEntry("Created by the "
	                    + descriptiveName + " tool.")
	        outputRaster.addMetadataEntry("Created on " + new Date())
			stop = new Date()
			td = TimeCategory.minus(stop, start)
			outputRaster.addMetadataEntry("Elapsed time: $td")
			outputRaster.close();

			// display the output image
			pluginHost.returnData(outputFile)
			


//			for (row = 0; row < rows + 2; row++) {
//				output[row][0] = nodata
//				output[row][cols + 1] = nodata
//				flowdir.setValue(row, 0, 0)
//				flowdir.setValue(row, cols + 1, 0)
//			}
//			
//			for (col = 0; col < cols + 2; col++) {
//				output[0][col] = nodata
//				output[rows + 1][col] = nodata
//				flowdir.setValue(0, col, 0)
//				flowdir.setValue(rows + 1, col, 0)
//			}
//
//			// now breach
//
//            if (!maxLengthUsed) {

//            } else {
//            	while (queue.isEmpty() == false) {
//	                gc = queue.poll();
//	                row = gc.row;
//	                col = gc.col;
//	                z = gc.z;
//	                
//	                for (i = 0; i < 8; i++) {
//	                    rowN = row + dY[i];
//	                    colN = col + dX[i];
//	                    zN = output[rowN][colN];
//	                    if ((zN != nodata) && (!inQueue.getValue(rowN, colN))) {
//	                        flowdir.setValue(rowN, colN, backLink[i])
//	                        if (pits.getValue(rowN, colN)) {
//	                        	// trace the flowpath back until you find a lower cell
//	                        	numCellsInPath = 0
//	                        	zTest = zN
//	                        	r = rowN
//	                        	c = colN
//	                        	flag = true
//	                        	while (flag) {
//	                        		zTest -= SMALL_NUM // ensures a small increment slope
//	                        		dir = flowdir.getValue(r, c)
//	                        		if (dir > 0) {
//	                        			r += dY[dir - 1]
//		                            	c += dX[dir - 1]
//		                            	zN2 = output[r][c]
//		                            	if (zN2 <= zTest || zN2 == nodata) {
//		                            		// a lower grid cell has been found
//		                            		flag = false
//		                            	}
//	                        		} else {
//	                        			flag = false
//	                        		}
//	                        		numCellsInPath++
//	                        		if (numCellsInPath > maxLength) { flag = false }
//	                        	}
//
//								if (numCellsInPath <= maxLength) {
//									zTest = zN
//		                        	r = rowN
//		                        	c = colN
//		                        	flag = true
//		                        	while (flag) {
//		                        		zTest -= SMALL_NUM // ensures a small increment slope
//		                        		dir = flowdir.getValue(r, c)
//		                        		if (dir > 0) {
//		                        			r += dY[dir - 1]
//			                            	c += dX[dir - 1]
//			                            	zN2 = output[r][c]
//			                            	if (zN2 <= zTest || zN2 == nodata) {
//			                            		// a lower grid cell has been found
//			                            		flag = false
//			                            	} else {
//			                            		output[r][c] = zTest
//			                            		// this cell is already in the 
//			                            		// queue but with a higher elevation
//			                            		gc = new GridCell(r, c, zTest);
//		                        				queue.add(gc);
//		                        				numCellsTotal++
//			                            	}
//		                        		} else {
//		                        			flag = false
//		                        		}
//		                        	}
//								}
//	                        	pits.setValue(rowN, colN, false)
//	                        }
//	                        numSolvedCells++;
//	                        gc = new GridCell(rowN, colN, zN);
//	                        queue.add(gc);
//	                        inQueue.setValue(rowN, colN, true)
//	                    }
//	                }
//	                progress = (int) (100f * numSolvedCells / numCellsTotal);
//	                if (progress > oldProgress) {
//	                    pluginHost.updateProgress("Loop 2 of 3", progress)
//	                    oldProgress = progress;
//	                    // check to see if the user has requested a cancellation
//						if (pluginHost.isRequestForOperationCancelSet()) {
//							pluginHost.showFeedback("Operation cancelled")
//							return
//						}
//	                }
//	            }
//            }
//
//            // output the data
//			WhiteboxRaster outputRaster = new WhiteboxRaster(outputFile, "rw", 
//  		  	     demFile, DataType.FLOAT, nodata)
//			outputRaster.setPreferredPalette(paletteName)
//  		  	WhiteboxRaster pointer
//  		  	if (outputPointer) {
//  		  		pointer = new WhiteboxRaster(pointerFile, "rw", 
//  		  	     demFile, DataType.FLOAT, nodata)
//  		  	    pointer.setDataScale(DataScale.CATEGORICAL)
//  		  	    pointer.setPreferredPalette("qual.pal")
//  		  	}
//
//			oldProgress = -1
//  		  	for (row = 0; row < rows; row++) {
//				for (col = 0; col < cols; col++) {
//					z = output[row + 1][col + 1]
//					outputRaster.setValue(row, col, z)
//					if (outputPointer) {
//						if (z != nodata) {
//							pointer.setValue(row, col, outPointer[flowdir.getValue(row + 1, col + 1)])
//						} else {
//							pointer.setValue(row, col, nodata)
//						}
//					}
//				}
//  		  		progress = (int)(100f * row / rowsLessOne)
//				if (progress > oldProgress) {
//					pluginHost.updateProgress("Loop 3 of 3", progress)
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
//			outputRaster.addMetadataEntry("Created by the "
//	                    + descriptiveName + " tool.")
//	        outputRaster.addMetadataEntry("Created on " + new Date())
//			outputRaster.close()
//	
//			if (outputPointer) {
//				pointer.addMetadataEntry("Created by the "
//	                    + descriptiveName + " tool.")
//		        pointer.addMetadataEntry("Created on " + new Date())
//		        pointer.close()
//		
//				// display the output image
//				pluginHost.returnData(pointerFile)
//			}
//
//			// display the output image
//			pluginHost.returnData(outputFile)
//			
//			if (performFlowAccumulation) {
//				String[] args2 = [pointerFile, flowAccumFile, "number of upslope grid cells", "false"] 
//				pluginHost.runPlugin("FlowAccumD8", args2, false)
//			}
			
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
	def tdf = new FillBurn(pluginHost, args, name, descriptiveName)
}
