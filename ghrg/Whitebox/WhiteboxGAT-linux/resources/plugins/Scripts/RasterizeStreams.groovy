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
import java.beans.PropertyChangeEvent
import java.beans.PropertyChangeListener
import java.util.Date
import java.util.ArrayList
import java.util.PriorityQueue
import java.util.HashSet;
import java.text.DecimalFormat
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.WhiteboxRasterBase.DataScale
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.ShapeFileRecord
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.ui.plugin_dialog.*
import whitebox.utilities.StringUtilities;
import whitebox.structures.BooleanBitArray2D;
import whitebox.structures.NibbleArray2D;
import whitebox.structures.DoubleArray2D;
import whitebox.structures.IntArray2D;
import whitebox.structures.BoundingBox;
import whitebox.structures.XYPoint;
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "RasterizeStreams"
def descriptiveName = "Rasterize Streams"
def description = "Rasterizes vector streams."
def toolboxes = ["StreamAnalysis"]

public class RasterizeStreams implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public RasterizeStreams(WhiteboxPluginHost pluginHost, 
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
			sd.addDialogFile("Input streams file", "Input Streams File:", "open", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogFile("Input base raster", "Input Base Raster:", "open", "Raster Files (*.dep), DEP", true, false)
            sd.addDialogFile("Output file", "Output File:", "save", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogCheckBox("Use NoData value for background?", "Use NoData for background?", true)
			sd.addDialogCheckBox("Use feature number as output value?", "Use feature number as output value?", false)
			sd.addDialogCheckBox("Output stream collisions/adjacencies?", "Output stream collisions/adjacencies?", false)
			
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
	  		//GridCell gc
	  		double LARGE_NUM = Float.MAX_VALUE
			int numInflowing;
  		  	double s, sN;
  		  	double backgroundValue = 0.0;
  		  	double nodata = -32768.0;
  		  	boolean outputFeatureNum = false;
  		  	boolean outputCollisions = false;
  		  	//double outValue;
  		  	DecimalFormat df = new DecimalFormat("###,###,###,###")
			DecimalFormat df2 = new DecimalFormat("##0.0#%")
			
			
	  		/*
	  		 * 7  8  1
	  		 * 6  X  2
	  		 * 5  4  3
	  		 */
//			int[] dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ]
//			int[] dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ]
			int[] dX4 = [ 1, 0, -1, 0 ]
			int[] dY4 = [ 0, 1, 0, -1 ]
			int[] backLink = [5, 6, 7, 8, 1, 2, 3, 4]
			double[] outPointer = [0, 1, 2, 4, 8, 16, 32, 64, 128]
			
			if (args.length < 3) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String streamsFile = args[0]
			String baseFile = args[1]
			String outputFile = args[2]
			if (args.length > 3 && !args[3].toLowerCase().equals("not specified")) {
				if (args[3].toLowerCase().contains("t")) {
					backgroundValue = nodata;
				}
			}

			if (args.length > 4 && !args[4].toLowerCase().equals("not specified")) {
				if (args[4].toLowerCase().contains("t")) {
					outputFeatureNum = true;
				}
			}

			if (args.length > 5 && !args[5].toLowerCase().equals("not specified")) {
				if (args[5].toLowerCase().contains("t")) {
					outputCollisions = true;
				}
			}
			
			// read the input image
			WhiteboxRaster outputRaster = new WhiteboxRaster(outputFile, "rw", 
  		  	     baseFile, DataType.FLOAT, nodata)
  		  	outputRaster.setForceAllDataInMemory(true)
			outputRaster.setPreferredPalette("spectrum.plt")
			outputRaster.setNoDataValue(nodata)
			int rows = outputRaster.getNumberRows()
			int cols = outputRaster.getNumberColumns()
			int rowsLessOne = rows - 1
			int colsLessOne = cols - 1

			//IntArray2D linkID = new IntArray2D(rows, cols, -32768)
			//DoubleArray2D linkPosition = new DoubleArray2D(rows, cols, nodata)
			BooleanBitArray2D linkEndNodes = new BooleanBitArray2D(rows, cols)
//			BooleanBitArray2D linkCollisions = new BooleanBitArray2D(rows, cols)
//			BooleanBitArray2D linkAdjacencies = new BooleanBitArray2D(rows, cols)
			
			// perform vector-to-raster conversion
			ShapeFile input = new ShapeFile(streamsFile)
			ShapeType shapeType = input.getShapeType()
            if (shapeType != ShapeType.POLYLINE) {
            	pluginHost.showFeedback("The input shapefile should be of a POLYLINE ShapeType.")
            	return
            }
            
            int numFeatures = input.getNumberOfRecords()
        	int count = 0
			double[][] points
			int startingPointInPart, endingPointInPart
			int i
			double x1, y1, x2, y2, xPrime, yPrime, d;
			BoundingBox box;
			int topRow, bottomRow, leftCol, rightCol;
			double rowYCoord, colXCoord;

			int numLinkCollisions = 0;
			int numStreamCells = 0;
			int recNum, numPoints, numParts, part;
//			int featureNum = 0;
			int[] partData;

			oldProgress = -1
			for (ShapeFileRecord record : input.records) {
				recNum = record.getRecordNumber()
//				if (outputFeatureNum) {
//					outValue = recNum;
//				} else {
//					outValue = 1.0;
//				}
                points = record.getGeometry().getPoints()
				numPoints = points.length;
				partData = record.getGeometry().getParts()
				numParts = partData.length
				for (part = 0; part < numParts; part++) {
//					featureNum++
					box = new BoundingBox();             
					startingPointInPart = partData[part];
                    if (part < numParts - 1) {
                        endingPointInPart = partData[part + 1];
                    } else {
                        endingPointInPart = numPoints;
                    }

					row = outputRaster.getRowFromYCoordinate(points[startingPointInPart][1])
					col = outputRaster.getColumnFromXCoordinate(points[startingPointInPart][0]);
					if (outputRaster.getValue(row, col) == nodata) {
                        linkEndNodes.setValue(row, col, true);
                        outputRaster.setValue(row, col, recNum); //outValue);
                        numStreamCells++;
					}
                    
                    row = outputRaster.getRowFromYCoordinate(points[endingPointInPart-1][1])
					col = outputRaster.getColumnFromXCoordinate(points[endingPointInPart-1][0]);
					if (outputRaster.getValue(row, col) == nodata) {
                        linkEndNodes.setValue(row, col, true);
                        outputRaster.setValue(row, col, recNum); //outValue);
                        numStreamCells++;
					}
                    
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
                    topRow = outputRaster.getRowFromYCoordinate(box.getMaxY());
                    bottomRow = outputRaster.getRowFromYCoordinate(box.getMinY());
                    leftCol = outputRaster.getColumnFromXCoordinate(box.getMinX());
                    rightCol = outputRaster.getColumnFromXCoordinate(box.getMaxX());

					// find each intersection with a row.
                    for (row = topRow; row <= bottomRow; row++) {

                        rowYCoord = outputRaster.getYCoordinateFromRow(row);
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
                                    col = outputRaster.getColumnFromXCoordinate(xPrime);

									if (outputRaster.getValue(row, col) == nodata) {
	                                    outputRaster.setValue(row, col, recNum); //outValue);
	                                    numStreamCells++;
									} else if (outputRaster.getValue(row, col) != recNum && !linkEndNodes.getValue(row, col)) {
										numLinkCollisions++;
//										linkCollisions.setValue(row, col, true);
									}
                                }
                            }
                        }
                    }

                    // find each intersection with a column.
                    for (col = leftCol; col <= rightCol; col++) {
                        colXCoord = outputRaster.getXCoordinateFromColumn(col);
                        for (i = startingPointInPart; i < endingPointInPart - 1; i++) {
                            if (isBetween(colXCoord, points[i][0], points[i + 1][0])) {
                                x1 = points[i][0];
                                x2 = points[i + 1][0];
                                if (x1 != x2) {
                                    y1 = points[i][1];
                                    y2 = points[i + 1][1];

                                    // calculate the intersection point
                                    yPrime = y1 + (colXCoord - x1) / (x2 - x1) * (y2 - y1);

                                    row = outputRaster.getRowFromYCoordinate(yPrime);
                                    if (outputRaster.getValue(row, col) == nodata) {
	                                    outputRaster.setValue(row, col, recNum); //outValue);
	                                    numStreamCells++;
                                    } else if (outputRaster.getValue(row, col) != recNum && !linkEndNodes.getValue(row, col)) {
										numLinkCollisions++;
//										linkCollisions.setValue(row, col, true);
									}
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

            // now count the number of stream cell adjacencies
            int numAdjacencies = 0;
            double id, idN;
            oldProgress = -1
			for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					id = outputRaster.getValue(row, col)
					if (id != nodata && !linkEndNodes.getValue(row, col)) { //  it's a stream cell
						HashSet<Integer> hs = new HashSet<>(4);
						for (n = 0; n < 4; n++) {
                    		rowN = row + dY4[n];
                    		colN = col + dX4[n];
							idN = outputRaster.getValue(rowN, colN)
							if (idN != id && idN != nodata && !linkEndNodes.getValue(rowN, colN)) {
								hs.add(idN);
							}
						}
						if (hs.size() > 0) {
							numAdjacencies += hs.size();
//							linkAdjacencies.setValue(row, col, true);
						}
					}
				}
				progress = (int)(100f * row / rowsLessOne)
				if (progress > oldProgress) {
					pluginHost.updateProgress("Counting stream adjacencies...", progress)
					oldProgress = progress

					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}

			if (!outputFeatureNum) {
				oldProgress = -1
				for (row = 0; row < rows; row++) {
					for (col = 0; col < cols; col++) {
						z = outputRaster.getValue(row, col);
						if (z != nodata && z > 0) { //  it's a stream cell
							outputRaster.setValue(row, col, 1);
						} else {
							outputRaster.setValue(row, col, backgroundValue);
						}
						
//						if (linkAdjacencies.getValue(row, col)) {
//							outputRaster.setValue(row, col, 1.0);
//						} else if (linkCollisions.getValue(row, col)) {
//							outputRaster.setValue(row, col, 2.0);
//						} else {
//							outputRaster.setValue(row, col, nodata);
//						}
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
			}

			outputRaster.setPreferredPalette("qual.plt")
			outputRaster.setDataScale(DataScale.CATEGORICAL)
			outputRaster.findMinAndMaxVals()
			outputRaster.addMetadataEntry("Created by the "
	                    + descriptiveName + " tool.")
	        outputRaster.addMetadataEntry("Created on " + new Date())
			outputRaster.close()
	
			// display the output image
			pluginHost.returnData(outputFile)


			if (outputCollisions) {
				numAdjacencies = (int)(numAdjacencies / 2)
	
				String str = "There were ${df.format(numLinkCollisions)} (${df2.format((double)numLinkCollisions/numStreamCells)}) stream collisions \nand ${df.format(numAdjacencies)} (${df2.format((double)numAdjacencies/numStreamCells)}) erroneous stream adjacencies \nduring rasterization."
				pluginHost.returnData(str)
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

//	@CompileStatic
//    class GridCell implements Comparable<GridCell> {
//
//        public int row;
//        public int col;
//        public double z;
//        public int flatIndex;
//        public boolean streamVal;
//        //public long priority
//
//        public GridCell(int row, int col, double z, int flatIndex, boolean streamVal) { //, long priority) { // double Z) {
//            this.row = row;
//            this.col = col;
//            this.z = z;
//            this.flatIndex = flatIndex;
//            this.streamVal = streamVal;
//        }
//
//        @Override
//        public int compareTo(GridCell other) {
//        	if (this.streamVal && !other.streamVal) {
//        		return -1
//        	} else if (!this.streamVal && other.streamVal) {
//        		return 1
//        	} else {
//	        	if (this.z > other.z) {
//	        		return 1
//	        	} else if (this.z < other.z) {
//	        		return -1
//	        	} else {
//	        		if (this.flatIndex > other.flatIndex) {
//	        			return 1
//	        		} else if (this.flatIndex < other.flatIndex) {
//	        			return -1
//	        		}
//					return 0
//	        	}
//        	}
//        }
//    }

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
	def tdf = new RasterizeStreams(pluginHost, args, name, descriptiveName)
}
