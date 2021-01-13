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
import java.util.concurrent.Future
import java.util.concurrent.*
import java.util.Collections;
import java.util.Comparator;
import java.io.File
import java.util.Date
import java.util.Collections
import java.util.ArrayList
import java.util.Arrays
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.ShapeTypeDimension
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.ui.plugin_dialog.*
import whitebox.utilities.FileUtilities
import whitebox.geospatialfiles.shapefile.attributes.AttributeTable
import whitebox.geospatialfiles.shapefile.attributes.DBFField
import whitebox.utilities.Topology
import com.vividsolutions.jts.geom.Coordinate
import com.vividsolutions.jts.geom.Geometry
import com.vividsolutions.jts.geom.GeometryFactory
import com.vividsolutions.jts.geom.LineString
import com.vividsolutions.jts.geom.MultiPolygon
import groovy.transform.CompileStatic


// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "SinglepartsToMultiparts"
def descriptiveName = "Single-parts to Multi-parts"
def description = "Converts a vector to consist of multi-part features."
def toolboxes = ["VectorTools", "ConversionTools"]

public class SinglepartsToMultiparts implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
	
    public SinglepartsToMultiparts(WhiteboxPluginHost pluginHost, 
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
            DialogFile df = sd.addDialogFile("Input shapefile", "Input Shapefile:", "open", "Shapefiles (*.shp), SHP", true, false)
            String[] listItems = []
            DialogComboBox idKey = sd.addDialogComboBox("Select the field that is the ID key.", "ID field", listItems, 0)
			sd.addDialogFile("Output file", "Output Vector File:", "save", "Vector Files (*.shp), SHP", true, false)
			DialogCheckBox cb = sd.addDialogCheckBox("Search for polygon holes?", "Search for polygon holes?", false)
			cb.setVisible(false)
			def lstrDF = { evt -> if (evt.getPropertyName().equals("value")) { 
            		String fileName = df.getValue()
            		if (fileName != null && !fileName.isEmpty()) { 
            			def file = new File(fileName)
            			if (file.exists()) {
            				ShapeFile sf = new ShapeFile(fileName)
            				if (sf.getShapeType().getBaseType() == ShapeType.POLYGON) {
            					cb.setVisible(true)
            				} else {
            					cb.setVisible(false)
            				}
	            			//AttributeTable table = new AttributeTable(fileName)
	            			String[] fields = sf.getAttributeTable().getAttributeTableFieldNames()
	            			String[] allFields = new String[fields.length + 1]
	            			allFields[fields.length] = "merge all"
	            			for (int i in 0..<fields.length) {
	            				allFields[i] = fields[i]
	            			}
	            			idKey.setListItems(allFields)
            			}
            		}
            	}
            } as PropertyChangeListener
            df.addPropertyChangeListener(lstrDF)

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
            int i, progress, oldProgress, numPoints
            
            if (args.length < 3) {
                pluginHost.showFeedback("Incorrect number of arguments given to tool.")
                return
            }
            
            // read the input parameters
			String inputFile = args[0]
			String fieldName = ""
			if (!args[1].toLowerCase().trim().contains("merge all")) {
				fieldName = args[1]
			}
			String outputFile = args[2]
			boolean searchForHoles = false
			if (args.length >= 4) {
				searchForHoles = Boolean.parseBoolean(args[3])
			}
			
			def input = new ShapeFile(inputFile)
			ShapeType shapeType = input.getShapeType()

			if (fieldName.trim().isEmpty()) {
				// there is no ID field to base the merging on
				int numFeatures = input.getNumberOfRecords()
	            DBFField[] fields = new DBFField[1];

	            fields[0] = new DBFField();
	            fields[0].setName("FID");
	            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
	            fields[0].setFieldLength(10);
	            fields[0].setDecimalCount(0);

				ShapeFile output

				if (shapeType.getBaseType() == ShapeType.POINT) {
					// points need cannot handle multipart features. Have to use multipoints
					if (shapeType.getDimension() == ShapeTypeDimension.XY) {
						output = new ShapeFile(outputFile, ShapeType.MULTIPOINT, fields);
					} else if (shapeType.getDimension() == ShapeTypeDimension.Z) {
						output = new ShapeFile(outputFile, ShapeType.MULTIPOINTZ, fields);
					} else { // ShapeTypeDimension.M
						output = new ShapeFile(outputFile, ShapeType.MULTIPOINTM, fields);
					}
					
					int FID = 1
					int j = 0
					oldProgress = -1
						
					int numParts = 0
					def partsList = new ArrayList<Integer>()
					int lastPartStartingValue = 0
					PointsList points = new PointsList()
					for (int f in 0..<numFeatures) {
						whitebox.geospatialfiles.shapefile.Geometry g = input.getRecord(f).getGeometry()
						
						double[][] xyPoints = g.getPoints()
						if (shapeType.getDimension() == ShapeTypeDimension.XY) {
							for (int p in 0..<xyPoints.length) {
								points.addPoint(xyPoints[p][0], xyPoints[p][1]) 
							}
						} else if (shapeType.getDimension() == ShapeTypeDimension.Z) {
							double[] zVals = ((GeometryZ)g).getzArray()
							double[] mVals = ((GeometryZ)g).getmArray()
							for (int p in 0..<xyPoints.length) {
								points.addZPoint(xyPoints[p][0], xyPoints[p][1], zVals[p], mVals[p]) 
							}
						} else if (shapeType.getDimension() == ShapeTypeDimension.M) {
							double[] mVals = ((GeometryM)g).getmArray()
							for (int p in 0..<xyPoints.length) {
								points.addMPoint(xyPoints[p][0], xyPoints[p][1], mVals[p]) 
							}
						}
						
						int[] parts = g.getParts()
						numParts += parts.length
						for (int p in 0..<parts.length) {
							partsList.add(parts[p] + lastPartStartingValue)									
						}
						lastPartStartingValue += xyPoints.length

						progress = (int)(100f * f / (numFeatures - 1))
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
					
					int[] parts = new int[numParts]
					for (int p in 0..<numParts) {
						parts[p] = partsList.get(p)
					}

					Object[] recData = new Object[1]
        			recData[0] = new Double(1)
        		
					whitebox.geospatialfiles.shapefile.Geometry g
					
					switch (shapeType.getDimension()) {
						case ShapeTypeDimension.XY:
							g = new MultiPoint(points.getPointsArray())
        					break
        				case ShapeTypeDimension.Z:
							g = (whitebox.geospatialfiles.shapefile.Geometry)(new MultiPointZ(points.getPointsArray(), points.getZArray(), points.getMArray()))
        					break
        				case ShapeTypeDimension.M:
							g = new MultiPointM(points.getPointsArray(), points.getMArray())
        					break
					}

					output.addRecord(g, recData)

				} else if (shapeType.getBaseType() == ShapeType.POLYGON && 
				      searchForHoles) {

				    output = new ShapeFile(outputFile, shapeType, fields);
					int FID = 1
					int j = 0
					oldProgress = -1
					
					// retrieve the JTS geometries
					ArrayList<com.vividsolutions.jts.geom.Polygon> polyList = new ArrayList<>()
					com.vividsolutions.jts.geom.Polygon poly
					for (int f in 0..<numFeatures) {
						whitebox.geospatialfiles.shapefile.Geometry g = input.getRecord(f).getGeometry()
						com.vividsolutions.jts.geom.Geometry[] jtsG = g.getJTSGeometries()
						for (com.vividsolutions.jts.geom.Geometry jg : jtsG) {
							poly = (com.vividsolutions.jts.geom.Polygon)jg
							polyList.add(poly)
							poly.setSRID(polyList.size())
						}

						progress = (int)(100f * f / (numFeatures - 1))
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

					int numParts = polyList.size()

					// sort the polygons by area
					Collections.sort(polyList, new Comparator<com.vividsolutions.jts.geom.Polygon>() {

		                @Override
		                public int compare(com.vividsolutions.jts.geom.Polygon o1, com.vividsolutions.jts.geom.Polygon o2) {
		                    Double area1 = o1.getArea();
		                    Double area2 = o2.getArea();
		                    return area2.compareTo(area1);
		                }
		                
		            });

					int[] containingPoly = new int[numParts]
					// initialize containingPoly with -1s
					for (int m = 0; m < numParts; m++) {
		            	containingPoly[m] = -1
		            }
		            // iterate through, finding polys that are contained in larger polys
		            for (int m = numParts - 1; m >= 0; m--) {
		                com.vividsolutions.jts.geom.Polygon item1 = polyList.get(m);
		               	for (int n = m - 1; n >= 0; n--) {
		                    com.vividsolutions.jts.geom.Polygon item2 = polyList.get(n);
		                    if (item2.contains(item1)) {
								containingPoly[m] = n
		                        break
		                    }
		                }
		            }
					// figure out if they are holes or islands
		            boolean[] isHole = new boolean[numParts]
		            for (int m = 0; m < numParts; m++) {
		            	if (containingPoly[m] > -1) {
		            		if (!isHole[containingPoly[m]]) {
		            			isHole[m] = true
		            		}
		            	}
		            }
		            
					def partsList = new ArrayList<Integer>()
					int[] parts = new int[numParts]
					int lastPartStartingValue = 0
					PointsList points = new PointsList()
					for (int f in 0..<numParts) {
						Coordinate[] coords = polyList.get(f).getExteriorRing().getCoordinates()
						if (shapeType.getDimension() == ShapeTypeDimension.XY) {
							if (!isHole[f]) {
			                    if (!Topology.isClockwisePolygon(coords)) {
			                        for (int k = coords.length - 1; k >= 0; k--) {
			                            points.addPoint(coords[k].x, coords[k].y);
			                        }
			                    } else {
			                        for (Coordinate coord : coords) {
			                            points.addPoint(coord.x, coord.y);
			                        }
			                    }
							} else {
								if (Topology.isClockwisePolygon(coords)) {
			                        for (int k = coords.length - 1; k >= 0; k--) {
			                            points.addPoint(coords[k].x, coords[k].y);
			                        }
			                    } else {
			                        for (Coordinate coord : coords) {
			                            points.addPoint(coord.x, coord.y);
			                        }
			                    }
							}
						} else if (shapeType.getDimension() == ShapeTypeDimension.Z) {
							if (!isHole[f]) {
			                    if (!Topology.isClockwisePolygon(coords)) {
			                        for (int k = coords.length - 1; k >= 0; k--) {
			                            points.addZPoint(coords[k].x, coords[k].y, coords[k].z);
			                        }
			                    } else {
			                        for (Coordinate coord : coords) {
			                            points.addZPoint(coord.x, coord.y, coord.z);
			                        }
			                    }
							} else {
								if (Topology.isClockwisePolygon(coords)) {
			                        for (int k = coords.length - 1; k >= 0; k--) {
			                            points.addZPoint(coords[k].x, coords[k].y, coords[k].z);
			                        }
			                    } else {
			                        for (Coordinate coord : coords) {
			                            points.addZPoint(coord.x, coord.y, coord.z);
			                        }
			                    }
							}
						} else if (shapeType.getDimension() == ShapeTypeDimension.M) {
							if (!isHole[f]) {
			                    if (!Topology.isClockwisePolygon(coords)) {
			                        for (int k = coords.length - 1; k >= 0; k--) {
			                            points.addMPoint(coords[k].x, coords[k].y, coords[k].z);
			                        }
			                    } else {
			                        for (Coordinate coord : coords) {
			                            points.addMPoint(coord.x, coord.y, coord.z);
			                        }
			                    }
							} else {
								if (Topology.isClockwisePolygon(coords)) {
			                        for (int k = coords.length - 1; k >= 0; k--) {
			                            points.addMPoint(coords[k].x, coords[k].y, coords[k].z);
			                        }
			                    } else {
			                        for (Coordinate coord : coords) {
			                            points.addMPoint(coord.x, coord.y, coord.z);
			                        }
			                    }
							}
						}
						
						parts[f] = lastPartStartingValue
	                    lastPartStartingValue = points.size()
	                    
					}
					
					Object[] recData = new Object[1]
        			recData[0] = new Double(FID)
        			FID++
        			
					whitebox.geospatialfiles.shapefile.Geometry g
					
					switch (shapeType) {
						case ShapeType.POLYGON:
							g = new Polygon(parts, points.getPointsArray())
        					break
        				case ShapeType.POLYGONZ:
							g = (whitebox.geospatialfiles.shapefile.Geometry)(new PolygonZ(parts, points.getPointsArray(), points.getZArray()))
        					break
        				case ShapeType.POLYGONM:
							g = new PolygonM(parts, points.getPointsArray(), points.getMArray())
        					break
					}

					output.addRecord(g, recData)
					
				} else {
					output = new ShapeFile(outputFile, shapeType, fields);
					oldProgress = -1
					
					int numParts = 0
					def partsList = new ArrayList<Integer>()
					int lastPartStartingValue = 0
					PointsList points = new PointsList()
					for (int f in 0..<numFeatures) {
						whitebox.geospatialfiles.shapefile.Geometry g = input.getRecord(f).getGeometry()
						
						double[][] xyPoints = g.getPoints()
						if (shapeType.getDimension() == ShapeTypeDimension.XY) {
							for (int p in 0..<xyPoints.length) {
								points.addPoint(xyPoints[p][0], xyPoints[p][1]) 
							}
						} else if (shapeType.getDimension() == ShapeTypeDimension.Z) {
							double[] zVals = ((GeometryZ)g).getzArray()
							double[] mVals = ((GeometryZ)g).getmArray()
							for (int p in 0..<xyPoints.length) {
								points.addZPoint(xyPoints[p][0], xyPoints[p][1], zVals[p], mVals[p]) 
							}
						} else if (shapeType.getDimension() == ShapeTypeDimension.M) {
							double[] mVals = ((GeometryM)g).getmArray()
							for (int p in 0..<xyPoints.length) {
								points.addMPoint(xyPoints[p][0], xyPoints[p][1], mVals[p]) 
							}
						}
						
						int[] parts = g.getParts()
						numParts += parts.length
						for (int p in 0..<parts.length) {
							partsList.add(parts[p] + lastPartStartingValue)									
						}
						lastPartStartingValue += xyPoints.length

						progress = (int)(100f * f / (numFeatures - 1))
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
					
					int[] parts = new int[numParts]
					for (int p in 0..<numParts) {
						parts[p] = partsList.get(p)
					}

					Object[] recData = new Object[1]
        			recData[0] = new Double(1)
        			
					whitebox.geospatialfiles.shapefile.Geometry g
					
					switch (shapeType) {
						case ShapeType.POLYLINE:
							g = new PolyLine(parts, points.getPointsArray())
        					break
        				case ShapeType.POLYLINEZ:
							g = new PolyLineZ(parts, points.getPointsArray(), points.getZArray(), points.getMArray())
        					break
        				case ShapeType.POLYLINEM:
							g = new PolyLineM(parts, points.getPointsArray(), points.getMArray())
        					break
        				case ShapeType.POLYGON:
							g = new Polygon(parts, points.getPointsArray())
        					break
        				case ShapeType.POLYGONZ:
							g = (whitebox.geospatialfiles.shapefile.Geometry)(new PolygonZ(parts, points.getPointsArray(), points.getZArray(), points.getMArray()))
        					break
        				case ShapeType.POLYGONM:
							g = new PolygonM(parts, points.getPointsArray(), points.getMArray())
        					break
        				case ShapeType.MULTIPOINT:
							g = new MultiPoint(points.getPointsArray())
        					break
        				case ShapeType.MULTIPOINTZ:
							g = (whitebox.geospatialfiles.shapefile.Geometry)(new MultiPointZ(points.getPointsArray(), points.getZArray(), points.getMArray()))
        					break
        				case ShapeType.MULTIPOINTM:
							g = new MultiPointM(points.getPointsArray(), points.getMArray())
        					break
					}

					output.addRecord(g, recData)
				}
				
            	output.write()
				pluginHost.returnData(outputFile)




				
			} else { // there is a merging ID field
	            int numFeatures = input.getNumberOfRecords()
	            AttributeTable table = input.getAttributeTable()
	            DBFField[] inputFields = table.getAllFields();

	            int keyFieldNum = table.getFieldColumnNumberFromName(fieldName)
	            if (keyFieldNum < 0) {
	            	pluginHost.showFeedback("Could not locate field name.")
	            	return
	            }

				DBFField[] fields = new DBFField[2];

	            fields[0] = new DBFField();
	            fields[0].setName("FID");
	            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
	            fields[0].setFieldLength(10);
	            fields[0].setDecimalCount(0);

	            fields[1] = inputFields[keyFieldNum]

	
				Map<Object, ArrayList<Integer>> hm = new HashMap<Object, ArrayList<Integer>>()
            	for (i = 0; i < numFeatures; i++) {
                	Object val = table.getValue(i, keyFieldNum)
                	if (hm.containsKey(val)) {
                		hm.get(val).add(i)
                	} else {
                		ArrayList<Integer> myList = new ArrayList<Integer>()
                		myList.add(i)
                		hm.put(val, myList)
                	}
            	}
				
				ShapeFile output

				if (shapeType.getBaseType() == ShapeType.POINT) {
					// points need cannot handle multipart features. Have to use multipoints
					if (shapeType.getDimension() == ShapeTypeDimension.XY) {
						output = new ShapeFile(outputFile, ShapeType.MULTIPOINT, fields);
					} else if (shapeType.getDimension() == ShapeTypeDimension.Z) {
						output = new ShapeFile(outputFile, ShapeType.MULTIPOINTZ, fields);
					} else { // ShapeTypeDimension.M
						output = new ShapeFile(outputFile, ShapeType.MULTIPOINTM, fields);
					}
					
					int FID = 1
					int j = 0
					oldProgress = -1
					hm.each() { it -> 
						ArrayList<Integer> myList = it.value
						int numParts = 0
						def partsList = new ArrayList<Integer>()
						int lastPartStartingValue = 0
						PointsList points = new PointsList()
						for (int f in 0..<myList.size()) {
							whitebox.geospatialfiles.shapefile.Geometry g = input.getRecord(myList.get(f)).getGeometry()
							
							double[][] xyPoints = g.getPoints()
							if (shapeType.getDimension() == ShapeTypeDimension.XY) {
								for (int p in 0..<xyPoints.length) {
									points.addPoint(xyPoints[p][0], xyPoints[p][1]) 
								}
							} else if (shapeType.getDimension() == ShapeTypeDimension.Z) {
								double[] zVals = ((GeometryZ)g).getzArray()
								double[] mVals = ((GeometryZ)g).getmArray()
								for (int p in 0..<xyPoints.length) {
									points.addZPoint(xyPoints[p][0], xyPoints[p][1], zVals[p], mVals[p]) 
								}
							} else if (shapeType.getDimension() == ShapeTypeDimension.M) {
								double[] mVals = ((GeometryM)g).getmArray()
								for (int p in 0..<xyPoints.length) {
									points.addMPoint(xyPoints[p][0], xyPoints[p][1], mVals[p]) 
								}
							}
							
							int[] parts = g.getParts()
							numParts += parts.length
							for (int p in 0..<parts.length) {
								partsList.add(parts[p] + lastPartStartingValue)									
							}
							lastPartStartingValue += xyPoints.length
							
						}
						
						int[] parts = new int[numParts]
						for (int p in 0..<numParts) {
							parts[p] = partsList.get(p)
						}

						Object[] recData = new Object[2]
            			recData[0] = new Double(FID)
            			FID++
            			recData[1] = it.key
            			
						whitebox.geospatialfiles.shapefile.Geometry g
						
						switch (shapeType.getDimension()) {
							case ShapeTypeDimension.XY:
								g = new MultiPoint(points.getPointsArray())
            					break
            				case ShapeTypeDimension.Z:
								g = (whitebox.geospatialfiles.shapefile.Geometry)(new MultiPointZ(points.getPointsArray(), points.getZArray(), points.getMArray()))
            					break
            				case ShapeTypeDimension.M:
								g = new MultiPointM(points.getPointsArray(), points.getMArray())
            					break
						}

						output.addRecord(g, recData)

            			j++
            			progress = (int)(100f * j / (hm.size() - 1))
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
				} else if (shapeType.getBaseType() == ShapeType.POLYGON && 
				      searchForHoles) {
				    output = new ShapeFile(outputFile, shapeType, fields);
					int FID = 1
					int j = 0
					oldProgress = -1
					hm.each() { it -> 
						ArrayList<Integer> myList = it.value

						// retrieve the JTS geometries
						ArrayList<com.vividsolutions.jts.geom.Polygon> polyList = new ArrayList<>()
						com.vividsolutions.jts.geom.Polygon poly
						for (int f in 0..<myList.size()) {
							whitebox.geospatialfiles.shapefile.Geometry g = input.getRecord(myList.get(f)).getGeometry()
							com.vividsolutions.jts.geom.Geometry[] jtsG = g.getJTSGeometries()
							for (com.vividsolutions.jts.geom.Geometry jg : jtsG) {
								poly = (com.vividsolutions.jts.geom.Polygon)jg
								polyList.add(poly)
								poly.setSRID(polyList.size())
							}
						}

						int numParts = polyList.size()

						// sort the polygons by area
						Collections.sort(polyList, new Comparator<com.vividsolutions.jts.geom.Polygon>() {

			                @Override
			                public int compare(com.vividsolutions.jts.geom.Polygon o1, com.vividsolutions.jts.geom.Polygon o2) {
			                    Double area1 = o1.getArea();
			                    Double area2 = o2.getArea();
			                    return area2.compareTo(area1);
			                }
			                
			            });

						int[] containingPoly = new int[numParts]
						// initialize containingPoly with -1s
						for (int m = 0; m < numParts; m++) {
			            	containingPoly[m] = -1
			            }
			            // iterate through, finding polys that are contained in larger polys
			            for (int m = numParts - 1; m >= 0; m--) {
			                com.vividsolutions.jts.geom.Polygon item1 = polyList.get(m);
			               	for (int n = m - 1; n >= 0; n--) {
			                    com.vividsolutions.jts.geom.Polygon item2 = polyList.get(n);
			                    if (item2.contains(item1)) {
									containingPoly[m] = n
			                        break
			                    }
			                }
			            }
						// figure out if they are holes or islands
			            boolean[] isHole = new boolean[numParts]
			            for (int m = 0; m < numParts; m++) {
			            	if (containingPoly[m] > -1) {
			            		if (!isHole[containingPoly[m]]) {
			            			isHole[m] = true
			            		}
			            	}
			            }
			            
						def partsList = new ArrayList<Integer>()
						int[] parts = new int[numParts]
						int lastPartStartingValue = 0
						PointsList points = new PointsList()
						for (int f in 0..<numParts) {
							Coordinate[] coords = polyList.get(f).getExteriorRing().getCoordinates()
							if (shapeType.getDimension() == ShapeTypeDimension.XY) {
								if (!isHole[f]) {
				                    if (!Topology.isClockwisePolygon(coords)) {
				                        for (int k = coords.length - 1; k >= 0; k--) {
				                            points.addPoint(coords[k].x, coords[k].y);
				                        }
				                    } else {
				                        for (Coordinate coord : coords) {
				                            points.addPoint(coord.x, coord.y);
				                        }
				                    }
								} else {
									if (Topology.isClockwisePolygon(coords)) {
				                        for (int k = coords.length - 1; k >= 0; k--) {
				                            points.addPoint(coords[k].x, coords[k].y);
				                        }
				                    } else {
				                        for (Coordinate coord : coords) {
				                            points.addPoint(coord.x, coord.y);
				                        }
				                    }
								}
							} else if (shapeType.getDimension() == ShapeTypeDimension.Z) {
								if (!isHole[f]) {
				                    if (!Topology.isClockwisePolygon(coords)) {
				                        for (int k = coords.length - 1; k >= 0; k--) {
				                            points.addZPoint(coords[k].x, coords[k].y, coords[k].z);
				                        }
				                    } else {
				                        for (Coordinate coord : coords) {
				                            points.addZPoint(coord.x, coord.y, coord.z);
				                        }
				                    }
								} else {
									if (Topology.isClockwisePolygon(coords)) {
				                        for (int k = coords.length - 1; k >= 0; k--) {
				                            points.addZPoint(coords[k].x, coords[k].y, coords[k].z);
				                        }
				                    } else {
				                        for (Coordinate coord : coords) {
				                            points.addZPoint(coord.x, coord.y, coord.z);
				                        }
				                    }
								}
							} else if (shapeType.getDimension() == ShapeTypeDimension.M) {
								if (!isHole[f]) {
				                    if (!Topology.isClockwisePolygon(coords)) {
				                        for (int k = coords.length - 1; k >= 0; k--) {
				                            points.addMPoint(coords[k].x, coords[k].y, coords[k].z);
				                        }
				                    } else {
				                        for (Coordinate coord : coords) {
				                            points.addMPoint(coord.x, coord.y, coord.z);
				                        }
				                    }
								} else {
									if (Topology.isClockwisePolygon(coords)) {
				                        for (int k = coords.length - 1; k >= 0; k--) {
				                            points.addMPoint(coords[k].x, coords[k].y, coords[k].z);
				                        }
				                    } else {
				                        for (Coordinate coord : coords) {
				                            points.addMPoint(coord.x, coord.y, coord.z);
				                        }
				                    }
								}
							}
							
							parts[f] = lastPartStartingValue
		                    lastPartStartingValue = points.size()
		                    
						}
						
						Object[] recData = new Object[2]
            			recData[0] = new Double(FID)
            			FID++
            			recData[1] = it.key
            			
						whitebox.geospatialfiles.shapefile.Geometry g
						
						switch (shapeType) {
							case ShapeType.POLYGON:
								g = new Polygon(parts, points.getPointsArray())
            					break
            				case ShapeType.POLYGONZ:
								g = (whitebox.geospatialfiles.shapefile.Geometry)(new PolygonZ(parts, points.getPointsArray(), points.getZArray()))
            					break
            				case ShapeType.POLYGONM:
								g = new PolygonM(parts, points.getPointsArray(), points.getMArray())
            					break
						}

						output.addRecord(g, recData)

            			j++
            			progress = (int)(100f * j / (hm.size() - 1))
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
				} else {
					output = new ShapeFile(outputFile, shapeType, fields);
					int FID = 1
					int j = 0
					oldProgress = -1
					hm.each() { it -> 
						ArrayList<Integer> myList = it.value
						int numParts = 0
						def partsList = new ArrayList<Integer>()
						int lastPartStartingValue = 0
						PointsList points = new PointsList()
						for (int f in 0..<myList.size()) {
							whitebox.geospatialfiles.shapefile.Geometry g = input.getRecord(myList.get(f)).getGeometry()
							
							double[][] xyPoints = g.getPoints()
							if (shapeType.getDimension() == ShapeTypeDimension.XY) {
								for (int p in 0..<xyPoints.length) {
									points.addPoint(xyPoints[p][0], xyPoints[p][1]) 
								}
							} else if (shapeType.getDimension() == ShapeTypeDimension.Z) {
								double[] zVals = ((GeometryZ)g).getzArray()
								double[] mVals = ((GeometryZ)g).getmArray()
								for (int p in 0..<xyPoints.length) {
									points.addZPoint(xyPoints[p][0], xyPoints[p][1], zVals[p], mVals[p]) 
								}
							} else if (shapeType.getDimension() == ShapeTypeDimension.M) {
								double[] mVals = ((GeometryM)g).getmArray()
								for (int p in 0..<xyPoints.length) {
									points.addMPoint(xyPoints[p][0], xyPoints[p][1], mVals[p]) 
								}
							}
							
							int[] parts = g.getParts()
							numParts += parts.length
							for (int p in 0..<parts.length) {
								partsList.add(parts[p] + lastPartStartingValue)									
							}
							lastPartStartingValue += xyPoints.length
							
						}
						
						int[] parts = new int[numParts]
						for (int p in 0..<numParts) {
							parts[p] = partsList.get(p)
						}

						Object[] recData = new Object[2]
            			recData[0] = new Double(FID)
            			FID++
            			recData[1] = it.key
            			
						whitebox.geospatialfiles.shapefile.Geometry g
						
						switch (shapeType) {
							case ShapeType.POLYLINE:
								g = new PolyLine(parts, points.getPointsArray())
            					break
            				case ShapeType.POLYLINEZ:
								g = new PolyLineZ(parts, points.getPointsArray(), points.getZArray(), points.getMArray())
            					break
            				case ShapeType.POLYLINEM:
								g = new PolyLineM(parts, points.getPointsArray(), points.getMArray())
            					break
            				case ShapeType.POLYGON:
								g = new Polygon(parts, points.getPointsArray())
            					break
            				case ShapeType.POLYGONZ:
								g = (whitebox.geospatialfiles.shapefile.Geometry)(new PolygonZ(parts, points.getPointsArray(), points.getZArray(), points.getMArray()))
            					break
            				case ShapeType.POLYGONM:
								g = new PolygonM(parts, points.getPointsArray(), points.getMArray())
            					break
            				case ShapeType.MULTIPOINT:
								g = new MultiPoint(points.getPointsArray())
            					break
            				case ShapeType.MULTIPOINTZ:
								g = (whitebox.geospatialfiles.shapefile.Geometry)(new MultiPointZ(points.getPointsArray(), points.getZArray(), points.getMArray()))
            					break
            				case ShapeType.MULTIPOINTM:
								g = new MultiPointM(points.getPointsArray(), points.getMArray())
            					break
						}

						output.addRecord(g, recData)

            			j++
            			progress = (int)(100f * j / (hm.size() - 1))
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
				
            	output.write()
				pluginHost.returnData(outputFile)
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
}

if (args == null) {
    pluginHost.showFeedback("Plugin arguments not set.")
} else {
    def f = new SinglepartsToMultiparts(pluginHost, args, name, descriptiveName)
}
