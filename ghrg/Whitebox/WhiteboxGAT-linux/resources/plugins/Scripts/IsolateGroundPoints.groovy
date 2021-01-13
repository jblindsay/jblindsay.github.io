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
import java.util.concurrent.Future
import java.util.concurrent.*
import java.util.concurrent.atomic.AtomicInteger
import java.util.Date
import java.util.ArrayList
import java.beans.PropertyChangeEvent
import java.beans.PropertyChangeListener
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.LASReader
import whitebox.geospatialfiles.LASReader.PointRecord
import whitebox.geospatialfiles.LASReader.PointRecColours
import whitebox.structures.BoundingBox
import whitebox.structures.KdTree
import whitebox.structures.BooleanBitArray1D
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.*;
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import whitebox.geospatialfiles.shapefile.attributes.*;
import whitebox.ui.plugin_dialog.*
import groovy.transform.CompileStatic
import groovy.time.TimeDuration
import groovy.time.TimeCategory

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "IsolateGroundPoints"
def descriptiveName = "Isolate Ground Points (LiDAR)"
def description = "Seperates ground points from off-terrain points."
def toolboxes = ["LidarTools"]

public class IsolateGroundPoints implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	private AtomicInteger numSolvedTiles = new AtomicInteger(0)
	
	public IsolateGroundPoints(WhiteboxPluginHost pluginHost, 
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
			//DialogMultiFile dmf = sd.addDialogMultiFile("Select the input point files", "Input Point Files:", "Point Files (*.las; *.shp), LAS, SHP")
			DialogFile df = sd.addDialogFile("Select the input point file", "Input Point Files:", "open", "Point Files (*.las; *.shp), LAS, SHP", true, false)
            String[] listItems = []
            DialogComboBox zKey = sd.addDialogComboBox("Select the field containing the elevation (Z) data.", "Elevation (Z) field", listItems, 0)
			zKey.setVisible(false)
			sd.addDialogFile("Output file", "Output Vector File:", "close", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogDataInput("Search Distance (m)", "Search Distance (m)", "", true, false)
			sd.addDialogDataInput("Minimum Number of neighbours in subsample.", "Min. Number of Points:", "0", true, false)
            sd.addDialogDataInput("Threshold in the slope between points to define an off-terrain point.", "Inter-point Slope Threshold:", "30.0", true, false)
			sd.addDialogDataInput("Minimum off-ground point height.", "Min. Off-Ground Point Height:", "0.15", true, false)
            sd.addDialogDataInput("Minimum elevation (used to exclude drop-out points; optional).", "Minimum Elevation:", "", true, true)
			sd.addDialogCheckBox("Display output file?", "Display output file?", false)
            
			def lstr = { evt -> if (evt.getPropertyName().equals("value")) { 
            		if (df.getValue() != null && !(df.getValue().isEmpty())) {
	            		String value = df.getValue()
	            		if (value != null && !value.isEmpty()) {
	            			if (value.contains(".shp")) {
	            				def fileName = value.replace(".shp", ".dbf")
	            				def file = new File(fileName)
		            			if (file.exists()) {
		            				def shapefile = new ShapeFile(value)
			            			AttributeTable table = shapefile.getAttributeTable() //new AttributeTable(fileName)
			            			if (shapefile.getShapeType() == ShapeType.MULTIPOINTZ) {
		            					zKey.setVisible(false)
			            			} else {
			            				String[] fields = table.getAttributeTableFieldNames()
			            				zKey.setListItems(fields)
			            				zKey.setVisible(true)
			            			}
		            			}
	            			} else {
	            				zKey.setVisible(false)
	            			}
	            		}
            		}
            	}
            } as PropertyChangeListener
            df.addPropertyChangeListener(lstr)

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
//		long start = System.currentTimeMillis()
		Date start = new Date();
		int progress = 0
	    int oldProgress = -1
	    double x, y, z, zN, dist, slope, maxDist, maxSlope
	    double higherZ, lowerZ, slopeThreshold
		int a, i, p, higherPointIndex
	    double[] entry
	    String zField = ""
	    String suffix = ""
	    String outputFile = ""
	    PointRecord point;
	    List<KdTree.Entry<InterpolationRecord>> results
	    InterpolationRecord value
	    boolean displayOutput = false
	    int numNeighbours = 0
		double minElevDiff = 0.15
		double elevDiff
		double minElev = Double.NEGATIVE_INFINITY

	  	try {
		  	// read the input parameters
			final String inputFile = args[0].trim()
			if (args.length == 9) {
				zField = args[1].trim()
				outputFile = args[2]
		        maxDist = Double.parseDouble(args[3]);
		        numNeighbours = Integer.parseInt(args[4]);
				maxSlope = Double.parseDouble(args[5]);
		        slopeThreshold = Math.tan(Math.toRadians(maxSlope));
		        minElevDiff = Double.parseDouble(args[6]);
		        if (!args[7].toLowerCase().contains("not specified")) {
		        	minElev = Double.parseDouble(args[7]);
		        }
		        displayOutput = Boolean.parseBoolean(args[8])
			} else if (args.length == 8) {
				outputFile = args[1]
		        maxDist = Double.parseDouble(args[2]);
		        numNeighbours = Integer.parseInt(args[3]);
				maxSlope = Double.parseDouble(args[4]);
		        slopeThreshold = Math.tan(Math.toRadians(maxSlope));
		        minElevDiff = Double.parseDouble(args[5]);
		        if (!args[6].toLowerCase().contains("not specified")) {
		        	minElev = Double.parseDouble(args[6]);
		        }
		        displayOutput = Boolean.parseBoolean(args[7])
			} else {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}

			BooleanBitArray1D isOffTerrain
			int numPoints = -1
			if (inputFile.toLowerCase().endsWith(".las")) {
				// create the LAS file reader
				LASReader las = new LASReader(inputFile);
	            numPoints = (int) las.getNumPointRecords();
	            isOffTerrain = new BooleanBitArray1D(numPoints);
	            
	            // Read the points into the k-dimensional tree.
	            byte classVal
	            KdTree<InterpolationRecord> pointsTree = new KdTree.SqrEuclid<InterpolationRecord>(2, new Integer(numPoints))
				for (a = 0; a < numPoints; a++) {
	                point = las.getPointRecord(a);
	                classVal = point.getClassification();
	                if (!point.isPointWithheld()
	                     && point.getReturnNumber() == point.getNumberOfReturns()
				 		 && classVal != 4 && classVal != 5 
				 		 && classVal != 6 && classVal != 7) {
	                    x = point.getX();
	                    y = point.getY();
	                    z = point.getZ();
	                    if (z >= minElev) {
	                    	entry = [x, y];
	                    	pointsTree.addPoint(entry, new InterpolationRecord(x, y, z, a));
	                    } else {
	                    	isOffTerrain.setValue(a, true)
	                    }
	                } else {
	                	isOffTerrain.setValue(a, true)
	                }
	                progress = (int) (100f * (a + 1) / numPoints);
	                if (progress != oldProgress) {
	                    oldProgress = progress;
	                    pluginHost.updateProgress("Reading point data:", progress);
	                    if (pluginHost.isRequestForOperationCancelSet()) {
	                        pluginHost.showFeedback("Operation cancelled")
							return
	                    }
	                }
	            }

	            /* visit the neighbourhood around each point hunting for 
	               points that are above other points with a slope greater 
	               than the threshold. */
	            oldProgress = -1
				for (a = 0; a < numPoints; a++) {
	            	if (!isOffTerrain.getValue(a)) {
						point = las.getPointRecord(a);
						classVal = point.getClassification();
		                if (!point.isPointWithheld()
	                     && point.getReturnNumber() == point.getNumberOfReturns()
				 		 && classVal != 4 && classVal != 5 
				 		 && classVal != 6 && classVal != 7) {
		                    x = point.getX();
		                    y = point.getY();
		                    z = point.getZ();
		                    
		                    entry = [x, y];
		                    results = pointsTree.neighborsWithinRange(entry, maxDist);
							int numResults = results.size()
		                    if (numResults < numNeighbours) {
		                    	results = pointsTree.nearestNeighbor(entry, numNeighbours, true, false)
		                    	numResults = results.size()
		                    }
	                		for (p = 0; p < numResults; p++) {
		                		dist = Math.sqrt(results.get(p).distance);
		                		if (dist > 0) {
		                			value = (InterpolationRecord)results.get(p).value;
		                			zN = value.z
		                			if (z > zN) {
		                				higherZ = z
		                				lowerZ = zN
		                				higherPointIndex = a
		                			} else {
		                				higherZ = zN
		                				lowerZ = z
		                				higherPointIndex = value.index
		                			}
		                			elevDiff = higherZ - lowerZ
									slope = elevDiff / dist //(higherZ - lowerZ) / dist
		                			if (slope > slopeThreshold && elevDiff > minElevDiff) {
		                				isOffTerrain.setValue(higherPointIndex, true)
									}
		                		}
		                	}
		                }
					}
					progress = (int) (100f * (a + 1) / numPoints);
	                if (progress != oldProgress) {
	                    oldProgress = progress;
	                    pluginHost.updateProgress("Finding ground points:", progress);
	                    if (pluginHost.isRequestForOperationCancelSet()) {
	                        pluginHost.showFeedback("Operation cancelled")
							return
	                    }
	                }
	            }

	            // how many ground points are there?
	            int numOutputPoints = 0;
	            for (a = 0; a < numPoints; a++) {
					if (!isOffTerrain.getValue(a)) {
						numOutputPoints++;
					}
	            }

	            // Create the output file
				DBFField[] fields = new DBFField[1];
	
	            fields[0] = new DBFField();
	            fields[0].setName("ID");
	            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
	            fields[0].setFieldLength(10);
	            fields[0].setDecimalCount(0);
	            
				ShapeFile output = new ShapeFile(outputFile, ShapeType.MULTIPOINTZ, fields);

				double[][] xyData = new double[numOutputPoints][2]
	            double[] zData = new double[numOutputPoints]
				double[] mData = new double[numOutputPoints]
				int j = 0;
				oldProgress = -1;
				for (a = 0; a < numPoints; a++) {
					if (!isOffTerrain.getValue(a)) {
						point = las.getPointRecord(a);
						xyData[j][0] = point.getX();
						xyData[j][1] = point.getY();
						zData[j] = point.getZ();
						mData[j] = point.getIntensity();
						j++;
					}
					progress = (int) (100f * (a + 1) / numPoints);
	                if (progress != oldProgress) {
	                    oldProgress = progress;
	                    pluginHost.updateProgress("Saving output:", progress);
	                    if (pluginHost.isRequestForOperationCancelSet()) {
	                        pluginHost.showFeedback("Operation cancelled")
							return
	                    }
	                }
	            }
	
	            MultiPointZ wbPoint = new MultiPointZ(xyData, zData, mData);
	            Object[] rowData = new Object[1];
	            rowData[0] = new Double(1.0);
	            output.addRecord(wbPoint, rowData);

//				// output the points
//				DBFField[] fields = new DBFField[2];
//
//	            fields[0] = new DBFField();
//	            fields[0].setName("Z");
//	            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
//	            fields[0].setFieldLength(10);
//	            fields[0].setDecimalCount(3);
//	
//	            fields[1] = new DBFField();
//	            fields[1].setName("I");
//	            fields[1].setDataType(DBFField.DBFDataType.NUMERIC);
//	            fields[1].setFieldLength(8);
//	            fields[1].setDecimalCount(0);
//
//		        // see if the output files already exist, and if so, delete them.
//		        File outFile = new File(outputFile);
//		        if (outFile.exists()) {
//		            outFile.delete();
//		            (new File(outputFile.replace(".shp", ".dbf"))).delete();
//		            (new File(outputFile.replace(".shp", ".shx"))).delete();
//		        }
//
//		        ShapeFile output = new ShapeFile(outputFile, ShapeType.POINT, fields);
//	
//				for (a = 0; a < numPoints; a++) {
//	            	if (!isOffTerrain.getValue(a)) {
//						point = las.getPointRecord(a);
//		                if (!point.isPointWithheld()) {
//		                    x = point.getX();
//		                    y = point.getY();
//		                    z = point.getZ();
//		                    whitebox.geospatialfiles.shapefile.Point wbGeometry = new whitebox.geospatialfiles.shapefile.Point(x, y);
//	
//		                    Object[] rowData = new Object[2];
//		                    rowData[0] = z;
//		                    rowData[1] = (double) point.getIntensity();
//							
//		                    output.addRecord(wbGeometry, rowData);
//		                }
//					}
//					progress = (int) (100f * (a + 1) / numPoints);
//	                if (progress != oldProgress) {
//	                    oldProgress = progress;
//	                    pluginHost.updateProgress("Outputting point data:", progress);
//	                    if (pluginHost.isRequestForOperationCancelSet()) {
//	                        pluginHost.showFeedback("Operation cancelled")
//							return
//	                    }
//	                }
//	            }
				
	            output.write()
	            
			} else if (inputFile.toLowerCase().endsWith(".shp")) {
				// create the ShapeFile
				ShapeFile input = new ShapeFile(inputFile);
				AttributeTable table = input.getAttributeTable()
				// make sure that it is of a POINT base shapetype
				if (input.getShapeType().getBaseType() != ShapeType.POINT &&
				    input.getShapeType() != ShapeType.MULTIPOINTZ) {
					pluginHost.showFeedback("The input file must be of a POINT shape type.")
					return
				}
				
				double[][] points
					
				if (input.getShapeType() == ShapeType.MULTIPOINTZ) {
					// how many points are there?
					for (ShapeFileRecord record : input.records) {
		                a = record.getRecordNumber() - 1
		                MultiPointZ mpz = (MultiPointZ)(record.getGeometry())
		                numPoints += mpz.getNumPoints()
		                
		            }
					isOffTerrain = new BooleanBitArray1D(numPoints);

		            
		            // Read the points into the k-dimensional tree.
		            KdTree<InterpolationRecord> pointsTree = new KdTree.SqrEuclid<InterpolationRecord>(2, new Integer(numPoints))
		            int n = 0
					for (ShapeFileRecord record : input.records) {
		                a = record.getRecordNumber() - 1
		                MultiPointZ mpz = (MultiPointZ)(record.getGeometry())
		                points = mpz.getPoints()
		                double[] zArray = mpz.getzArray()

		                for (int m = 0; m < points.length; m++) {
			              	x = points[m][0]
		                    y = points[m][1]
							z = zArray[m]  
			                
		                    entry = [x, y];
		                    pointsTree.addPoint(entry, new InterpolationRecord(x, y, z, n));
		                    n++
		                    progress = (int) (100f * n / numPoints);
			                if (progress != oldProgress) {
			                    oldProgress = progress;
			                    pluginHost.updateProgress("Reading point data:", progress);
			                    if (pluginHost.isRequestForOperationCancelSet()) {
			                        pluginHost.showFeedback("Operation cancelled")
									return
			                    }
			                }
		                }
		            }

					/* visit the neighbourhood around each point hunting for 
		               points that are above other points with a slope greater 
		               than the threshold. */
					n = 0
					int r = 0
		            for (ShapeFileRecord record : input.records) {
		                a = record.getRecordNumber() - 1
		                MultiPointZ mpz = (MultiPointZ)(record.getGeometry())
		                points = mpz.getPoints()
		                double[] zArray = mpz.getzArray()
						
		                for (int m = 0; m < points.length; m++) {
		                	
			                if (!isOffTerrain.getValue(n)) {
			                	points = record.getGeometry().getPoints()
			              		x = points[0][0]
			                    y = points[0][1]
								z = zArray[m]
			                    entry = [x, y];
			                    
			                    results = pointsTree.neighborsWithinRange(entry, maxDist);
								int numResults = results.size()
			                    if (numResults < numNeighbours) {
			                    	results = pointsTree.nearestNeighbor(entry, numNeighbours, true, false)
			                    	numResults = results.size()
			                    }
			                	for (p = 0; p < numResults; p++) {
			                		dist = Math.sqrt(results.get(p).distance);
			                		if (dist > 0) {
			                			value = (InterpolationRecord)results.get(p).value;
			                			zN = value.z
			                			if (z > zN) {
			                				higherZ = z
			                				lowerZ = zN
			                				higherPointIndex = n
			                			} else {
			                				higherZ = zN
			                				lowerZ = z
			                				higherPointIndex = value.index
			                			}
			                			elevDiff = higherZ - lowerZ
										slope = elevDiff / dist
			                			if (slope > slopeThreshold && elevDiff > minElevDiff) {
											isOffTerrain.setValue(higherPointIndex, true)
										}
			                		}
			                	}
			                }
			                n++
			                progress = (int) (100f * n / numPoints);
			                if (progress != oldProgress) {
			                    oldProgress = progress;
			                    pluginHost.updateProgress("Finding ground points:", progress);
			                    if (pluginHost.isRequestForOperationCancelSet()) {
			                        pluginHost.showFeedback("Operation cancelled")
									return
			                    }
			                }
		                }
		            }

		            // output the points
		            DBFField[] fields = table.getAllFields()
	
					//outputFile = inputFile.replace(".shp", suffix + ".shp");
	
			        // see if the output files already exist, and if so, delete them.
			        File outFile = new File(outputFile);
			        if (outFile.exists()) {
			            outFile.delete();
			            (new File(outputFile.replace(".shp", ".dbf"))).delete();
			            (new File(outputFile.replace(".shp", ".shx"))).delete();
			        }
		            
		            ShapeFile output = new ShapeFile(outputFile, ShapeType.MULTIPOINTZ, fields);
					n = 0
					for (ShapeFileRecord record : input.records) {
		                a = record.getRecordNumber() - 1

		                MultiPointZ mpz = (MultiPointZ)(record.getGeometry())
		                points = mpz.getPoints()
		                double[] zArray = mpz.getzArray()
						double[] mArray = mpz.getmArray()

						// how many valid points are there in this feature?
		                int numValidPoints = 0
			            for (int m = 0; m < points.length; m++) {
		                	if (!isOffTerrain.getValue(n)) {
			                	numValidPoints++
							}
		                }

		                double[][] xyData = new double[numValidPoints][2]
		                double[] zData = new double[numValidPoints]
						double[] mData = new double[numValidPoints]
						int q = 0
		                for (int m = 0; m < points.length; m++) {
		                	if (!isOffTerrain.getValue(n)) {
			                	xyData[q][0] = points[m][0]
								xyData[q][1] = points[m][1]
								zData[q] = zArray[m]
								mData[q] = mArray[m]
								q++
							}
							
							n++
							progress = (int) (100f * n / numPoints);
			                if (progress != oldProgress) {
			                    oldProgress = progress;
			                    pluginHost.updateProgress("Outputting point data:", progress);
			                    if (pluginHost.isRequestForOperationCancelSet()) {
			                        pluginHost.showFeedback("Operation cancelled")
									return
			                    }
			                }
		                }
		                
		                if (numValidPoints > 0) {
		                	MultiPointZ wbPoint = new MultiPointZ(xyData, zData, mData);
		            		Object[] rowData = table.getRecord(a);
		                    
		                    output.addRecord(wbPoint, rowData);
				        }
		            }
		            
		            output.write()
				} else {				
		            numPoints = (int) input.getNumberOfRecords();
		            isOffTerrain = new BooleanBitArray1D(numPoints);
		            
		            // Read the points into the k-dimensional tree.
		            KdTree<InterpolationRecord> pointsTree = new KdTree.SqrEuclid<InterpolationRecord>(2, new Integer(numPoints))
					for (ShapeFileRecord record : input.records) {
		                a = record.getRecordNumber() - 1
		                
		                points = record.getGeometry().getPoints()
		              	z = (double)table.getValue(a, zField)  
		                x = points[0][0]
	                    y = points[0][1]
						
	                    entry = [x, y];
	                    pointsTree.addPoint(entry, new InterpolationRecord(x, y, z, a));
		                
		                progress = (int) (100f * (a + 1) / numPoints);
		                if (progress != oldProgress) {
		                    oldProgress = progress;
		                    pluginHost.updateProgress("Reading point data:", progress);
		                    if (pluginHost.isRequestForOperationCancelSet()) {
		                        pluginHost.showFeedback("Operation cancelled")
								return
		                    }
		                }
		            }
	
	
		            /* visit the neighbourhood around each point hunting for 
		               points that are above other points with a slope greater 
		               than the threshold. */
		
		            for (ShapeFileRecord record : input.records) {
		                a = record.getRecordNumber() - 1
		                if (!isOffTerrain.getValue(a)) {
		                	points = record.getGeometry().getPoints()
		              		x = points[0][0]
		                    y = points[0][1]
							z = (double)table.getValue(a, zField)  
		                    entry = [x, y];
		                    
		                    results = pointsTree.neighborsWithinRange(entry, maxDist);
							int numResults = results.size()
		                    if (numResults < numNeighbours) {
		                    	results = pointsTree.nearestNeighbor(entry, numNeighbours, true, false)
		                    	numResults = results.size()
		                    }
		                	for (p = 0; p < numResults; p++) {
		                		dist = Math.sqrt(results.get(p).distance);
		                		if (dist > 0) {
		                			value = (InterpolationRecord)results.get(p).value;
		                			zN = value.z
		                			if (z > zN) {
		                				higherZ = z
		                				lowerZ = zN
		                				higherPointIndex = a
		                			} else {
		                				higherZ = zN
		                				lowerZ = z
		                				higherPointIndex = value.index
		                			}
		                			elevDiff = higherZ - lowerZ
									slope = elevDiff / dist //(higherZ - lowerZ) / dist
		                			if (slope > slopeThreshold && elevDiff > minElevDiff) {
										isOffTerrain.setValue(higherPointIndex, true)
									}
		                		}
		                	}
		                }
						progress = (int) (100f * (a + 1) / numPoints);
		                if (progress != oldProgress) {
		                    oldProgress = progress;
		                    pluginHost.updateProgress("Finding ground points:", progress);
		                    if (pluginHost.isRequestForOperationCancelSet()) {
		                        pluginHost.showFeedback("Operation cancelled")
								return
		                    }
		                }
		            }
	
		            // output the points
		            DBFField[] fields = table.getAllFields()
	
					//outputFile = inputFile.replace(".shp", suffix + ".shp");
	
			        // see if the output files already exist, and if so, delete them.
			        File outFile = new File(outputFile);
			        if (outFile.exists()) {
			            outFile.delete();
			            (new File(outputFile.replace(".shp", ".dbf"))).delete();
			            (new File(outputFile.replace(".shp", ".shx"))).delete();
			        }
		            
		            ShapeFile output = new ShapeFile(outputFile, ShapeType.POINT, fields);
		
					for (ShapeFileRecord record : input.records) {
		                a = record.getRecordNumber() - 1
		                if (!isOffTerrain.getValue(a)) {
		                	points = record.getGeometry().getPoints()
		            		x = points[0][0]
		                    y = points[0][1]
							whitebox.geospatialfiles.shapefile.Point wbGeometry = new whitebox.geospatialfiles.shapefile.Point(x, y);
	
		                    Object[] rowData = table.getRecord(a);
		                    
		                    output.addRecord(wbGeometry, rowData);
						}
						progress = (int) (100f * (a + 1) / numPoints);
		                if (progress != oldProgress) {
		                    oldProgress = progress;
		                    pluginHost.updateProgress("Outputting point data:", progress);
		                    if (pluginHost.isRequestForOperationCancelSet()) {
		                        pluginHost.showFeedback("Operation cancelled")
								return
		                    }
		                }
		            }
		            
		            output.write()
				}
			}

			Date stop = new Date()
			TimeDuration td = TimeCategory.minus(stop, start)
			pluginHost.showFeedback("Elapsed time: $td")

			if (displayOutput) {
				pluginHost.returnData(outputFile)
			} else {
            	pluginHost.showFeedback("Operation complete")
			}
		
	  	} catch (OutOfMemoryError oe) {
            pluginHost.showFeedback("An out-of-memory error has occurred during operation.")
	    } catch (Exception e) {
	        pluginHost.showFeedback("An error has occurred during operation. See log file for details.")
	        pluginHost.logException("Error in " + descriptiveName, e)
	    } finally {
	    	// reset the progress bar
	    	pluginHost.updateProgress("Progress", 0)
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

	class InterpolationRecord {
        
        double x;
        double y;
        double z;
        int index;
        
        InterpolationRecord(double x, double y, double z, int index) {
            this.x = x;
            this.y = y;
            this.z = z;
            this.index = index;
        }

//        double getZ() {
//            return z;
//        }
//        
//        int getIndex() {
//        	return index;
//        }
    }
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def myClass = new IsolateGroundPoints(pluginHost, args, name, descriptiveName)
}
