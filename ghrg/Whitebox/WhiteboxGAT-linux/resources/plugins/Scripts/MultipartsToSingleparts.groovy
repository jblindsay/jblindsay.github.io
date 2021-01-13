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
def name = "MultipartsToSingleparts"
def descriptiveName = "Multi-parts to Single-parts"
def description = "Converts a vector to consist of multi-part features."
def toolboxes = ["VectorTools", "ConversionTools"]

public class MultipartsToSingleparts implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
	
    public MultipartsToSingleparts(WhiteboxPluginHost pluginHost, 
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
            sd.addDialogFile("Output file", "Output Vector File:", "save", "Vector Files (*.shp), SHP", true, false)
			DialogCheckBox cb = sd.addDialogCheckBox("Do not separate polygon holes.", "Do not separate polygon holes", true)
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
            
            if (args.length < 2) {
                pluginHost.showFeedback("Incorrect number of arguments given to tool.")
                return
            }
            
            // read the input parameters
			String inputFile = args[0]
			String outputFile = args[1]
			boolean excludePolyHoles = true
			if (args.length >= 3) {
				excludePolyHoles = Boolean.parseBoolean(args[2])
			}
			def input = new ShapeFile(inputFile)
			ShapeType shapeType = input.getShapeType()
			if (shapeType.getBaseType() == ShapeType.POINT) {
				pluginHost.showFeedback("The input file does not contain multi-part features. Operation cancelled.")
				return
			}
			
			AttributeTable table = input.getAttributeTable()
			DBFField[] fields = table.getAllFields()
			def output = new ShapeFile(outputFile, shapeType, fields)

			int numRecs = input.getNumberOfRecords()
			int[] parts
			boolean[] partHoleData
			whitebox.geospatialfiles.shapefile.Geometry g
			whitebox.geospatialfiles.shapefile.Geometry gOut
			oldProgress = -1
			progress = 0;

			if (excludePolyHoles && shapeType.getBaseType() == ShapeType.POLYGON) {
				for (i = 0; i < numRecs; i++) {
					ShapeFileRecord record = input.getRecord(i)
				    if (record.getShapeType() != ShapeType.NULLSHAPE) {
				    	Object[] recData = table.getRecord(i)
				    	g = record.getGeometry()
				    	parts = g.getParts()
				    	
				    	if (parts.length == 1) {
				    		output.addRecord(g, recData);
				    	} else { // there is more than one part

							partHoleData = new boolean[parts.length]
					    	if (shapeType == ShapeType.POLYGON) {
								partHoleData = ((Polygon)g).getPartHoleData()
							} else if (shapeType == ShapeType.POLYGONZ) {
								partHoleData = ((PolygonZ)g).getPartHoleData()
							} else if (shapeType == ShapeType.POLYGONM) {
								partHoleData = ((PolygonM)g).getPartHoleData()
							}

							// get the JTS geometries
							ArrayList<com.vividsolutions.jts.geom.Polygon> polyList = new ArrayList<>()
							com.vividsolutions.jts.geom.Polygon poly
							com.vividsolutions.jts.geom.Geometry[] jtsG = g.getJTSGeometries()
							for (com.vividsolutions.jts.geom.Geometry jg : jtsG) {
								poly = (com.vividsolutions.jts.geom.Polygon)jg
								poly.setSRID(polyList.size())
								polyList.add(poly)
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
				            // iterate through, finding hole polys that are contained in larger polys
				            for (int m = numParts - 1; m >= 0; m--) {
				            	com.vividsolutions.jts.geom.Polygon item1 = polyList.get(m);
				               	if (partHoleData[item1.getSRID()]) {
					               	for (int n = m - 1; n >= 0; n--) {
					                    com.vividsolutions.jts.geom.Polygon item2 = polyList.get(n);
					                    if (item2.contains(item1)) {
											containingPoly[m] = n
					                        break
					                    }
					                }
				               	}
				            }

				            for (int m = numParts - 1; m >= 0; m--) {
				            	com.vividsolutions.jts.geom.Polygon item1 = polyList.get(m);
				               	if (containingPoly[m] == -1) {
				               		ArrayList<com.vividsolutions.jts.geom.Polygon> myHoles = new ArrayList<>()
					               	for (int n = m - 1; n >= 0; n--) {
					                    if (containingPoly[n] == m) {
											myHoles.add(polyList.get(n))
					                    }
					                }

									def partsList = new ArrayList<Integer>()
									int[] outParts = new int[myHoles.size() + 1]
									//int lastPartStartingValue = 0
									PointsList points = new PointsList()
									Coordinate[] coords = item1.getExteriorRing().getCoordinates()
				                	if (shapeType.getDimension() == ShapeTypeDimension.XY) {
				                		if (!Topology.isClockwisePolygon(coords)) {
					                        for (int k = coords.length - 1; k >= 0; k--) {
					                            points.addPoint(coords[k].x, coords[k].y);
					                        }
					                    } else {
					                        for (Coordinate coord : coords) {
					                            points.addPoint(coord.x, coord.y);
					                        }
					                    }
				                	} else if (shapeType.getDimension() == ShapeTypeDimension.Z) {
				                		if (!Topology.isClockwisePolygon(coords)) {
					                        for (int k = coords.length - 1; k >= 0; k--) {
					                            points.addZPoint(coords[k].x, coords[k].y, coords[k].z);
					                        }
					                    } else {
					                        for (Coordinate coord : coords) {
					                            points.addZPoint(coord.x, coord.y, coord.z);
					                        }
					                    }
				                	} else if (shapeType.getDimension() == ShapeTypeDimension.M) {
				                		if (!Topology.isClockwisePolygon(coords)) {
					                        for (int k = coords.length - 1; k >= 0; k--) {
					                            points.addMPoint(coords[k].x, coords[k].y, coords[k].z);
					                        }
					                    } else {
					                        for (Coordinate coord : coords) {
					                            points.addMPoint(coord.x, coord.y, coord.z);
					                        }
					                    }
				                	}
					                	
					                for (int q = 0; q < myHoles.size(); q++) {
					                	outParts[q + 1] = points.size()
					                	
					                	coords = myHoles.get(q).getExteriorRing().getCoordinates()

					                	if (shapeType.getDimension() == ShapeTypeDimension.XY) {
					                		if (Topology.isClockwisePolygon(coords)) {
						                        for (int k = coords.length - 1; k >= 0; k--) {
						                            points.addPoint(coords[k].x, coords[k].y);
						                        }
						                    } else {
						                        for (Coordinate coord : coords) {
						                            points.addPoint(coord.x, coord.y);
						                        }
						                    }
					                	} else if (shapeType.getDimension() == ShapeTypeDimension.Z) {
					                		if (Topology.isClockwisePolygon(coords)) {
						                        for (int k = coords.length - 1; k >= 0; k--) {
						                            points.addZPoint(coords[k].x, coords[k].y, coords[k].z);
						                        }
						                    } else {
						                        for (Coordinate coord : coords) {
						                            points.addZPoint(coord.x, coord.y, coord.z);
						                        }
						                    }
					                	} else if (shapeType.getDimension() == ShapeTypeDimension.M) {
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

									switch (shapeType) {
										case ShapeType.POLYGON:
											gOut = new Polygon(outParts, points.getPointsArray())
			            					break
			            				case ShapeType.POLYGONZ:
											gOut = (whitebox.geospatialfiles.shapefile.Geometry)(new PolygonZ(outParts, points.getPointsArray(), points.getZArray(), points.getMArray()))
			            					break
			            				default: //case ShapeType.POLYGONM:
											gOut = new PolygonM(outParts, points.getPointsArray(), points.getMArray())
			            					break
									}
							
									output.addRecord(gOut, recData)

				               	}
				            }
							
				    	}
				    }
				    	
					progress = (int)(100f * i / (numRecs - 1))
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
				for (i = 0; i < numRecs; i++) {
					ShapeFileRecord record = input.getRecord(i)
				    if (record.getShapeType() != ShapeType.NULLSHAPE) {
				    	Object[] recData = table.getRecord(i)
				    	g = record.getGeometry()
				    	parts = g.getParts()
				    	partHoleData = new boolean[parts.length]
				    	if (shapeType == ShapeType.POLYGON) {
							partHoleData = ((Polygon)g).getPartHoleData()
						} else if (shapeType == ShapeType.POLYGONZ) {
							partHoleData = ((PolygonZ)g).getPartHoleData()
						} else if (shapeType == ShapeType.POLYGONM) {
							partHoleData = ((PolygonM)g).getPartHoleData()
						}
				    	if (parts.length == 1) {
				    		output.addRecord(g, recData);
				    	} else { // there is more than one part
				    		double[][] xyPoints = g.getPoints()
				    		double[] zVals
				    		double[] mVals
							if (shapeType.getDimension() == ShapeTypeDimension.Z) {
								zVals = ((GeometryZ)g).getzArray()
								mVals = ((GeometryZ)g).getmArray()
							} else if (shapeType.getDimension() == ShapeTypeDimension.M) {
								mVals = ((GeometryM)g).getmArray()
							}	
				    		for (int p in 0..<parts.length) {
				    			//if (!partHoleData[p]) {
				    			PointsList points = new PointsList()
				    			int startingPointInPart = parts[p]
								int endingPointInPart
					            if (p < parts.length - 1) {
					                endingPointInPart = parts[p + 1]
					            } else {
					                endingPointInPart = xyPoints.length
					            }
					            int numPointsInPart = endingPointInPart - startingPointInPart
	
								if (!partHoleData[p]) {
									for (int k = startingPointInPart; k < endingPointInPart; k++) {
										if (shapeType.getDimension() == ShapeTypeDimension.XY) {
											points.addPoint(xyPoints[k][0], xyPoints[k][1])
										} else if (shapeType.getDimension() == ShapeTypeDimension.Z) {
											points.addZPoint(xyPoints[k][0], xyPoints[k][1], zVals[k], mVals[k])
										} else if (shapeType.getDimension() == ShapeTypeDimension.M) {
											points.addMPoint(xyPoints[k][0], xyPoints[k][1], mVals[k])
										}
						            }
								} else {
									for (int k = endingPointInPart - 1; k >= startingPointInPart; k--) {
										if (shapeType.getDimension() == ShapeTypeDimension.XY) {
											points.addPoint(xyPoints[k][0], xyPoints[k][1])
										} else if (shapeType.getDimension() == ShapeTypeDimension.Z) {
											points.addZPoint(xyPoints[k][0], xyPoints[k][1], zVals[k], mVals[k])
										} else if (shapeType.getDimension() == ShapeTypeDimension.M) {
											points.addMPoint(xyPoints[k][0], xyPoints[k][1], mVals[k])
										}
						            }
								}
	
					            int[] outParts = [0]
	
					            switch (shapeType) {
									case ShapeType.POLYLINE:
										gOut = new PolyLine(outParts, points.getPointsArray())
		            					break
		            				case ShapeType.POLYLINEZ:
										gOut = new PolyLineZ(outParts, points.getPointsArray(), points.getZArray(), points.getMArray())
		            					break
		            				case ShapeType.POLYLINEM:
										gOut = new PolyLineM(outParts, points.getPointsArray(), points.getMArray())
		            					break
		            				case ShapeType.POLYGON:
										gOut = new Polygon(outParts, points.getPointsArray())
		            					break
		            				case ShapeType.POLYGONZ:
										gOut = (whitebox.geospatialfiles.shapefile.Geometry)(new PolygonZ(outParts, points.getPointsArray(), points.getZArray(), points.getMArray()))
		            					break
		            				case ShapeType.POLYGONM:
										gOut = new PolygonM(outParts, points.getPointsArray(), points.getMArray())
		            					break
		            				case ShapeType.MULTIPOINT:
										gOut = new MultiPoint(points.getPointsArray())
		            					break
		            				case ShapeType.MULTIPOINTZ:
										gOut = (whitebox.geospatialfiles.shapefile.Geometry)(new MultiPointZ(points.getPointsArray(), points.getZArray(), points.getMArray()))
		            					break
		            				default: //case ShapeType.MULTIPOINTM:
										gOut = new MultiPointM(points.getPointsArray(), points.getMArray())
		            					break
								}
		
								output.addRecord(gOut, recData)
				    			//}
				    		}
				    	}
				    }
					
					progress = (int)(100f * i / (numRecs - 1))
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
    def f = new MultipartsToSingleparts(pluginHost, args, name, descriptiveName)
}
