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
import java.text.DecimalFormat
import java.io.File
import java.util.Date
import java.util.ArrayList
import java.util.PriorityQueue
import java.util.Arrays
import java.util.Collections
import java.beans.PropertyChangeEvent
import java.beans.PropertyChangeListener
import java.util.concurrent.Future
import java.util.concurrent.*
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.ShapeFileInfo
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.*
import whitebox.utilities.FileUtilities;
import com.vividsolutions.jts.geom.*
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase
import whitebox.geospatialfiles.VectorLayerInfo
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.geospatialfiles.shapefile.ShapeFileRecord
import whitebox.utilities.Topology
import whitebox.structures.BoundingBox
import whitebox.structures.RowPriorityGridCell
import whitebox.structures.KdTree
import whitebox.structures.SimpleGridCell
import whitebox.structures.DoubleArray2D
import whitebox.structures.IntArray2D
import groovy.transform.CompileStatic

def name = "DiscreteSibsonInterpolation"
def descriptiveName = "Discrete Sibson (Natural Neighbour) Interpolation"
def description = "Performs a discrete Sibson (natural neighbour) interpolation"
def toolboxes = ["Interpolation"]

public class DiscreteSibsonInterpolation implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
	
    public DiscreteSibsonInterpolation(WhiteboxPluginHost pluginHost, 
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
            def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "InterpolationIDW.groovy"
            sd.setSourceFile(scriptFile)
			
            // add some components to the dialog
        	//DialogFile dfIn = sd.addDialogFile("Input file", "Input Vector File:", "open", "Vector Files (*.shp), SHP", true, false)
            DialogFieldSelector dfs = sd.addDialogFieldSelector("Input file and height field.", "Input Height Field:", false)
            DialogCheckBox dcb = sd.addDialogCheckBox("Use z-values", "Use z-values", false)
            dcb.setVisible(false)
            sd.addDialogFile("Output file", "Output Raster File:", "saveAs", "Whitebox Raster Files (*.dep), DEP", true, false)
            sd.addDialogDataInput("Output raster cell size.", "Cell Size (optional):", "", true, true)
            sd.addDialogFile("Input base file", "Base Raster File (optional):", "open", "Whitebox Raster Files (*.dep), DEP", true, true)
            sd.addDialogDataInput("Discretization factor (5-20)", "Discretization factor (5-20)", "10", true, false)
			
            def listener = { evt -> if (evt.getPropertyName().equals("value")) { 
            		String value = dfs.getValue()
            		if (value != null && !value.isEmpty()) {
            			value = value.trim()
            			String[] strArray = dfs.getValue().split(";")
            			String fileName = strArray[0]
            			File file = new File(fileName)
            			if (file.exists()) {
	            			ShapeFileInfo shapefile = new ShapeFileInfo(fileName)
		            		if (shapefile.getShapeType().getDimension() == ShapeTypeDimension.Z) {
		            			dcb.setVisible(true)
		            		} else {
		            			dcb.setVisible(false)
		            		}
		            	} else {
		            		if (dcb.isVisible()) {
		            			dcb.setVisible(false)
		            		}
		            	}
            		}
            	} 
            } as PropertyChangeListener
            dfs.addPropertyChangeListener(listener)
            
            // resize the dialog to the standard size and display it
            sd.setSize(800, 400)
            sd.visible = true
        }
    }

    @CompileStatic
    private void execute(String[] args) {
        try {
        	int progress, oldProgress, rows, cols, row, col
        	double x, y, z
        	double north, south, east, west
        	double cellSize = -1.0
        	double nodata = -32768
        	ArrayList<Double> xList = new ArrayList<>()
            ArrayList<Double> yList = new ArrayList<>()
            ArrayList<Double> zList = new ArrayList<>()
	            
			if (args.length != 6) {
                pluginHost.showFeedback("Incorrect number of arguments given to tool.")
                return
            }
            
        	// read the input parameters
            String[] inputData = args[0].split(";")
            String inputFile = inputData[0]
			boolean useZValues = Boolean.parseBoolean(args[1])
            String outputFile = args[2]
			if (!args[3].toLowerCase().contains("not specified")) {
	            cellSize = Double.parseDouble(args[3]);
	        }
	        String baseFileHeader = args[4]
	        if (baseFileHeader == null || baseFileHeader.isEmpty()) {
	        	baseFileHeader = "not specified"
	        }
	        int discretizationFactor = (int)(Double.parseDouble(args[5]));
		
	        
            ShapeFile input = new ShapeFile(inputFile)

			if (cellSize < 0 && baseFileHeader.toLowerCase().equals("not specified")) {
				cellSize = Math.max((input.getxMax() - input.getxMin()) / 998, (input.getyMax() - input.getyMin()) / 998)
			}
            
			AttributeTable table = new AttributeTable(inputFile.replace(".shp", ".dbf"))
			ShapeType shapeType = input.getShapeType()
            if (shapeType.getDimension() != ShapeTypeDimension.Z && useZValues) {
            	useZValues = false
            }
			int heightField = -1
			String heightFieldName = ""
            if (inputData.length == 2 && !inputData[1].trim().isEmpty()) {
            	heightFieldName = inputData[1].trim();
            } else if (!useZValues) {
            	pluginHost.showFeedback("A field within the input file's attribute table must be selected to assign point heights.")
            	return
            }
			
			double[][] point
			Object[] recData
			Coordinate c
			GeometryFactory geomFactory = new GeometryFactory()
			int i = 0
			int numFeatures = input.getNumberOfRecords()
			oldProgress = -1
			if (!useZValues) {
				for (ShapeFileRecord record : input.records) {
					recData = table.getRecord(i)
					z = (Double)table.getValue(i, heightFieldName);
					point = record.getGeometry().getPoints()
					for (int p = 0; p < point.length; p++) {
						x = point[p][0]
						y = point[p][1]
						xList.add(x)
						yList.add(y)
						zList.add(z)
					}
					i++
	                progress = (int)(100f * i / numFeatures)
	            	if (progress != oldProgress) {
						pluginHost.updateProgress("Reading Points:", progress)
	            		oldProgress = progress
	            	}
	            	// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			} else {
				for (ShapeFileRecord record : input.records) {
					if (shapeType.getBaseType() == ShapeType.POINT) {
						PointZ ptz = (PointZ)(record.getGeometry())
                		z = ptz.getZ()
                		x = ptz.getX()
						y = ptz.getY()
						xList.add(x)
						yList.add(y)
						zList.add(z)
					} else if (shapeType.getBaseType() == ShapeType.MULTIPOINT) {
						MultiPointZ mptz = (MultiPointZ)(record.getGeometry())
						point = record.getGeometry().getPoints()
						double[] zArray = mptz.getzArray()
						for (int p = 0; p < point.length; p++) {
							x = point[p][0]
							y = point[p][1]
							z = zArray[p]
							xList.add(x)
							yList.add(y)
							zList.add(z)

							progress = (int)(100f * p / point.length)
			            	if (progress != oldProgress) {
								pluginHost.updateProgress("Reading Points:", progress)
			            		oldProgress = progress

			            		// check to see if the user has requested a cancellation
								if (pluginHost.isRequestForOperationCancelSet()) {
									pluginHost.showFeedback("Operation cancelled")
									return
								}
			            	}
						}
					} else if (shapeType.getBaseType() == ShapeType.POLYLINE) {
						PolyLineZ plz = (PolyLineZ)(record.getGeometry())
						point = record.getGeometry().getPoints()
						double[] zArray = plz.getzArray()
						for (int p = 0; p < point.length; p++) {
							x = point[p][0]
							y = point[p][1]
							z = zArray[p]
							xList.add(x)
							yList.add(y)
							zList.add(z)
							
							progress = (int)(100f * p / point.length)
			            	if (progress != oldProgress) {
								pluginHost.updateProgress("Reading Points:", progress)
			            		oldProgress = progress

			            		// check to see if the user has requested a cancellation
								if (pluginHost.isRequestForOperationCancelSet()) {
									pluginHost.showFeedback("Operation cancelled")
									return
								}
			            	}
						}
					} else if (shapeType.getBaseType() == ShapeType.POLYGON) {
						PolygonZ pz = (PolygonZ)(record.getGeometry())
						point = record.getGeometry().getPoints()
						double[] zArray = pz.getzArray()
						for (int p = 0; p < point.length; p++) {
							x = point[p][0]
							y = point[p][1]
							z = zArray[p]
							xList.add(x)
							yList.add(y)
							zList.add(z)

							progress = (int)(100f * p / point.length)
			            	if (progress != oldProgress) {
								pluginHost.updateProgress("Reading Points:", progress)
			            		oldProgress = progress

			            		// check to see if the user has requested a cancellation
								if (pluginHost.isRequestForOperationCancelSet()) {
									pluginHost.showFeedback("Operation cancelled")
									return
								}
			            	}
						}
					}
					
					i++
	                progress = (int)(100f * i / numFeatures)
	            	if (progress != oldProgress) {
						pluginHost.updateProgress("Reading Points:", progress)
	            		oldProgress = progress
	            	}
	            	// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}
			
			int numSamples = zList.size()
			KdTree<Double> pointsTree = new KdTree.SqrEuclid<Double>(2, new Integer(numSamples));
			pointsTree = new KdTree.SqrEuclid<Double>(2, new Integer(numSamples));
			oldProgress = -1;
			for (i = 0; i < numSamples; i++) {
				double[] entry = new double[2]
				entry[0] = yList.get(i)
				entry[1] = xList.get(i);
				pointsTree.addPoint(entry, zList.get(i));

				progress = (int)(100f * i / numSamples)
            	if (progress > oldProgress) {
            		oldProgress = progress
            		pluginHost.updateProgress("Building the kd-tree:", progress)
            	
            		// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}
			}

			xList.clear();
			yList.clear();
			zList.clear();
			
			north = input.getyMax() + cellSize / 2.0;
	        south = input.getyMin() - cellSize / 2.0;
	        east = input.getxMax() + cellSize / 2.0;
	        west = input.getxMin() - cellSize / 2.0;
                
			// initialize the output raster
            WhiteboxRaster output;
            if ((cellSize > 0) || ((cellSize < 0) & (baseFileHeader.toLowerCase().contains("not specified")))) {
                if ((cellSize < 0) & (baseFileHeader.toLowerCase().contains("not specified"))) {
                    cellSize = Math.min((input.getyMax() - input.getyMin()) / 500.0,
                            (input.getxMax() - input.getxMin()) / 500.0);
                }
                rows = (int) (Math.ceil((north - south) / cellSize));
                cols = (int) (Math.ceil((east - west) / cellSize));

                // update west and south
                east = west + cols * cellSize;
                south = north - rows * cellSize;

                output = new WhiteboxRaster(outputFile, north, south, east, west,
                    rows, cols, WhiteboxRasterBase.DataScale.CONTINUOUS,
                    WhiteboxRasterBase.DataType.FLOAT, nodata, nodata);
            } else {
                output = new WhiteboxRaster(outputFile, "rw",
                        baseFileHeader, WhiteboxRasterBase.DataType.FLOAT, nodata);
                rows = output.getNumberRows()
                cols = output.getNumberColumns()

                north = output.getNorth();
		        south = output.getSouth();
		        east = output.getEast();
		        west = output.getWest();
            }

            int rowsLessOne = rows - 1;
            int colsLessOne = cols - 1;
            double EWRange = east - west;
            double NSRange = north - south;
            
            output.setPreferredPalette("spectrum.pal")

            double resolutionX = output.getCellSizeX()
			double halfResolutionX = resolutionX / 2.0;

			double resolutionY = output.getCellSizeY()
			double halfResolutionY = resolutionY / 2.0;
			
			double[] entry
//            int numCells = rows * cols;
//			List<KdTree.Entry<SimpleGridCell>> results;
//			KdTree<SimpleGridCell> gridTree = new KdTree.SqrEuclid<SimpleGridCell>(2, new Integer(numCells))
//			oldProgress = -1
//			for (row = 0; row < rows; row++) {
//                y = (north - halfResolutionY) - (row * resolutionY); //output.getYCoordinateFromRow(row);
//                for (col = 0; col < cols; col++) {
//					x = (west + halfResolutionX) + (col * resolutionX); //output.getXCoordinateFromColumn(col);
//                    entry = [y, x]
//					gridTree.addPoint(entry, new SimpleGridCell(row, col));    
//                }
//
//                progress = (int)(100f * row / rowsLessOne)
//            	if (progress > oldProgress) {
//            		oldProgress = progress
//            		pluginHost.updateProgress("Initializing:", progress)
//            	
//            		// check to see if the user has requested a cancellation
//					if (pluginHost.isRequestForOperationCancelSet()) {
//						pluginHost.showFeedback("Operation cancelled")
//						return
//					}
//            	}
//			}

			IntArray2D nArray = new IntArray2D(rows, cols, 0);
			DoubleArray2D sumArray = new DoubleArray2D(rows, cols, nodata);

			int fineRows = rows * discretizationFactor;
	        int fineRowsLessOne = fineRows - 1;
			int fineCols = cols * discretizationFactor;
			double fineResolutionX = EWRange / fineCols;
			double fineHalfResolutionX = fineResolutionX / 2.0;
			double fineResolutionY = NSRange / fineRows ;
			double fineHalfResolutionY = fineResolutionY / 2.0;
			int row2, col2;
			double value;
			double distance;
			//pluginHost.showFeedback("$rows $cols $fineRows $fineCols \n $north $south $east $west");
						
			KdTree.Entry<Double> result;
			SimpleGridCell cellValue;
			double startX, startY, endX, endY;
			int startCol, startRow, endCol, endRow;
			double x2, y2, distance2;
			oldProgress = -1
			for (row = 0; row < fineRows; row++) {
                y = (north - fineHalfResolutionY) - (row * fineResolutionY);
            	for (col = 0; col < fineCols; col++) {
                	x = (west + fineHalfResolutionX) + (col * fineResolutionX);
                    entry = [y, x];
                    result = pointsTree.nearestNeighbor(entry);
                    value = (double)result.value; // the value of the nearest point
                    distance = Math.sqrt(result.distance) // the distance to the nearest point

					// define the bounding box
                    startX = x - distance;
                    endX = x + distance;
                    startY = y + distance;
                    endY = y - distance;

                    // find all the grid cell centres within the bounding box
                    startCol = (int)(Math.round(colsLessOne * (startX - west - halfResolutionX) / EWRange));
					endCol = (int)(Math.round(colsLessOne * (endX - west - halfResolutionX) / EWRange));
					startRow = (int)(Math.round(rowsLessOne * (north - halfResolutionY - startY) / NSRange));
					endRow = (int)(Math.round(rowsLessOne * (north - halfResolutionY - endY) / NSRange));

					// are any of the grid cell centres within distance of point x, y?
					for (row2 = startRow; row2 <= endRow; row2++) {
						y2 = (north - halfResolutionY) - (row2 * resolutionY);
						for (col2 = startCol; col2 <= endCol; col2++) {
							x2 = (col2 * resolutionX) + (west + halfResolutionX);
							distance2 = Math.sqrt((y2 - y) * (y2 - y) + (x2 - x) * (x2 - x));
							if (distance2 <= distance) {
								if (nArray.getValue(row2, col2) > 0) {
	                            	sumArray.incrementValue(row2, col2, value)
	                            } else {
	                            	sumArray.setValue(row2, col2, value)
	                            }
	                            nArray.incrementValue(row2, col2)
							}
						}
					}

//                    // see if there are any grid cells within this distance of the point
//					results = gridTree.neighborsWithinRange(entry, distance);
//					if (results.size() > 0) {
//						for (i = 0; i < results.size(); i++) {
//                        	cellValue = (SimpleGridCell)results.get(i).value;
//                            row2 = cellValue.row
//                            col2 = cellValue.column
//                            if (nArray.getValue(row2, col2) > 0) {
//                            	sumArray.incrementValue(row2, col2, value)
//                            } else {
//                            	sumArray.setValue(row2, col2, value)
//                            }
//                            nArray.incrementValue(row2, col2)
//                        }
//					}
                }
                progress = (int)(100f * row / fineRowsLessOne)
            	if (progress > oldProgress) {
            		oldProgress = progress
            		pluginHost.updateProgress("Interpolating:", progress)

            		// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}
			}

			oldProgress = -1
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                	output.setValue(row, col, nArray.getValue(row, col));
//                	if (nArray.getValue(row, col) > 0) {
//                		output.setValue(row, col, sumArray.getValue(row, col) / nArray.getValue(row, col));
//                	} else {
//                		// use the nearest neighbour
//                		x = (col * resolutionX) + (west + halfResolutionX); //output.getXCoordinateFromColumn(col);
//                		y = (north - halfResolutionY) - (row * resolutionY); //output.getYCoordinateFromRow(row);
//                		entry = [y, x]
//						result = pointsTree.nearestNeighbor(entry);
//                		value = (double)result.value;
//                
//                		output.setValue(row, col, value);
//                	}
                }

                progress = (int)(100f * row / rowsLessOne)
            	if (progress > oldProgress) {
            		oldProgress = progress
            		pluginHost.updateProgress("Saving data:", progress)
            	
            		// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}
            }
			
			output.addMetadataEntry("Created by the "
	                    + descriptiveName + " tool.")
	        output.addMetadataEntry("Created on " + new Date())
			output.addMetadataEntry("Discretization Factor = $discretizationFactor")
			output.close()

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
    def f = new DiscreteSibsonInterpolation(pluginHost, args, name, descriptiveName)
}
