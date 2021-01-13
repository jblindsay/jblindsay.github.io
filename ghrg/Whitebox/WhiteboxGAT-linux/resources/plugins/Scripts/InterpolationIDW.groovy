/*
 * Copyright (C) 2013 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D
import org.apache.commons.math3.geometry.euclidean.threed.Plane
import org.apache.commons.math3.geometry.euclidean.twod.Vector2D
import org.apache.commons.math3.linear.SingularValueDecomposition
import org.apache.commons.math3.linear.QRDecomposition
import org.apache.commons.math3.linear.Array2DRowRealMatrix
import org.apache.commons.math3.linear.ArrayRealVector
import org.apache.commons.math3.linear.RealVector
import org.apache.commons.math3.linear.RealMatrix
import org.apache.commons.math3.linear.DecompositionSolver
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "InterpolationIDW"
def descriptiveName = "Inverse Distance Weighted (IDW) Interpolation"
def description = "Performs an inverse distance weighted (IDW) interpolation"
def toolboxes = ["Interpolation"]

public class InterpolationIDW implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName

//    KdTree<Integer> pointsTree
//	ArrayList<Double> xList = new ArrayList<>()
//    ArrayList<Double> yList = new ArrayList<>()
//    ArrayList<Double> zList = new ArrayList<>()
//    WhiteboxRaster output
//    double maxSearchDist
//    int numNeighbours
//    boolean useQuadrandSearch
//    NodalFunction[] nodalFunc
//    WeightFunction weightFunction
	
    public InterpolationIDW(WhiteboxPluginHost pluginHost, 
        String[] args, def descriptiveName) {
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
            sd.setHelpFile("InterpolationIDW")
		
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
            sd.addDialogComboBox("Weight function", "Weight Function:", ["<html><i>d</i><sup>-<i>p</i></sup> (Shepard's method)</html>", "<html>[(<i>d<sub>max</sub> - d</i>) / (<i>d<sub>max</sub>d</i>)]<sup><i>p</i></sup> (Franke & Nielson)</html>", "<html>[(<i>d<sub>max</sub> - d</i>) / (<i>d<sub>max</sub></i>)]<sup>4<i>p</i></sup> (Lindsay)</html>"], 0)
            sd.addDialogDataInput("Exponent parameter value.", "<html>Exponent Parameter (<i>p</i>):</html>", "2.0", true, false)
			
            //sd.addDialogLabel("<html><b>Nodal Function Parameters</b></html>")
            DialogOption dopt = sd.addDialogOption("Nodal function", "Nodal Function:", "Constant (Shepard's method)", "Quadratic")
//			def optListener = { evt -> if (evt.getPropertyName().equals("value")) { 
//            		String value = dopt.getValue()
//            		if (value != null && !value.isEmpty()) {
//            			println(value)
//            		}
//            	} 
//            } as PropertyChangeListener
//            dopt.addPropertyChangeListener(optListener)
            
            
            sd.addDialogLabel("<html><b>Neighbourhood Parameters</b></html>")
            sd.addDialogDataInput("Search distance.", "<html>Search Distance (<i>d<sub>max</sub> ):</html>", "", true, false)
            sd.addDialogDataInput("Minimum Number of neighbours in subsample.", "Min. Number of Points:", "", true, false)
            sd.addDialogCheckBox("Use quadrant-based search", "Use quadrant-based search", false)
            
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
        	List<KdTree.Entry<Integer>> results
            ArrayList<Double> xList = new ArrayList<>()
            ArrayList<Double> yList = new ArrayList<>()
            ArrayList<Double> zList = new ArrayList<>()
	            
			if (args.length != 11) {
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
	        int weightType = 0
	        if (args[5].toLowerCase().contains("shepard")) {
	        	weightType = 0
	        } else if (args[5].toLowerCase().contains("franke")) {
	        	weightType = 1
	        } else if (args[5].toLowerCase().contains("lindsay")) {
	        	weightType = 2
	        }
	        double exponent = Double.parseDouble(args[6])
	        String nodalFunction = "constant"
	        if (args[7].toLowerCase().contains("quad")) {
	        	nodalFunction = "quadratic"
	        }
			double maxSearchDist = Double.parseDouble(args[8])
			if (maxSearchDist < 0) {
				maxSearchDist = 0
			}
			int numNeighbours = Integer.parseInt(args[9])
			if (numNeighbours < 0) {
				numNeighbours = 0
			}
			boolean useQuadrandSearch = Boolean.parseBoolean(args[10])

			WeightFunction weightFunction = new WeightFunction(weightType, exponent)
	        
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
            	heightFieldName = inputData[1].trim()
//            	String[] fieldNames = input.getAttributeTableFields()
//            	for (int i = 0; i < fieldNames.length; i++) {
//            		if (fieldNames[i].trim().equals(heightFieldName)) {
//            			heightField = i
//            			break
//            		}
//            	}
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
			KdTree<Integer> pointsTree = new KdTree.SqrEuclid<Integer>(2, new Integer(numSamples));
			pointsTree = new KdTree.SqrEuclid<Integer>(2, new Integer(numSamples));
			for (i = 0; i < numSamples; i++) {
				double[] entry = new double[2]
				entry[0] = xList.get(i)
				entry[1] = yList.get(i);
				pointsTree.addPoint(entry, i);
			}
			
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
            }
            
            output.setPreferredPalette("spectrum.pal")

			// calculate the nodal functions
			NodalFunction[] nodalFunc = new NodalFunction[numSamples]
			if (nodalFunction.contains("constant")) {
				byte functionType = 0
				for (i = 0; i < numSamples; i++) {
					double[] p = new double[1]
					p[0] = zList.get(i)
					nodalFunc[i] = new NodalFunction(xList.get(i), yList.get(i), functionType, p)
				}
			} else { // QUADRATIC
				byte functionType = 1
				int targetNumNeighbours = Math.max(8, numNeighbours)
				oldProgress = -1
				for (i = 0; i < numSamples; i++) {
					x = xList.get(i)
                    y = yList.get(i)
                    
					double[] entry = new double[2]
                    entry[0] = x
                    entry[1] = y
                    
                    results = pointsTree.neighborsWithinRange(entry, maxSearchDist)
					int numResults = results.size()
					if (numResults < targetNumNeighbours) {
                    	results = pointsTree.nearestNeighbor(entry, targetNumNeighbours, true, useQuadrandSearch)
                    	numResults = results.size()
                    }
                    double[] p = new double[6]
                    if (numResults >= 6) { // need at least 5 points to solve the quadratic
						double[] xVals = new double[numResults - 1]
						double[] yVals = new double[numResults - 1]
						double[] zVals = new double[numResults - 1]
						int m = 0
						for (int j = 0; j < numResults; j++) {
							KdTree.Entry result = results.get(j)
							if (result.distance > 0) {
								int k = result.value
								xVals[m] = xList.get(k)
								yVals[m] = yList.get(k)
								zVals[m] = zList.get(k)
								m++
							}
						}
	                    p = solveQuadratic(x, y, zList.get(i), xVals, yVals, zVals)
                    } else {
                    	pluginHost.showFeedback("Error during script execution")
                    	return
                    }
					nodalFunc[i] = new NodalFunction(x, y, functionType, p)
					
					progress = (int)(100f * i / (numSamples - 1))
					if (progress > oldProgress) {
						pluginHost.updateProgress("Calculating Nodal Functions:", progress)
						oldProgress = progress
					}
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}

//			// leave-one-out cross validation
//			boolean loocv = true
//			if (loocv) {
//				double dMax = maxSearchDist
//				double testP
//				double sumP = 0
//				int N = 0
//				oldProgress = -1
//				for (i = 0; i < numSamples; i++) {
//					x = xList.get(i)
//                    y = yList.get(i)
//                    double ptZ = zList.get(i)
//                    double minError = Double.POSITIVE_INFINITY
//                    double minErrorP = -1
//					double[] entry = new double[2]
//                    entry[0] = x
//                    entry[1] = y
//                    
//                    results = pointsTree.neighborsWithinRange(entry, maxSearchDist)
//					int numResults = results.size()
//					if (numResults < numNeighbours) {
//                    	results = pointsTree.nearestNeighbor(entry, numNeighbours, true, useQuadrandSearch)
//                    	numResults = results.size()
//                    }
//                    if (numResults >= 6) { // need at least 5 points to solve the quadratic
//						int[] neighbourIDs = new int[numResults - 1]
//						double[] distances = new double[numResults - 1]
//						int m = 0
//						for (int j = 0; j < numResults; j++) {
//							KdTree.Entry result = results.get(j)
//							if (result.distance > 0) {
//								neighbourIDs[m] = result.value
//								distances[m] = Math.sqrt(result.distance)
//								m++
//							}
//						}
//	                    for (int p = 1; p <= 1000; p++) {
//	                    	testP = p / 10.0
//	                    	weightFunction.setP(testP)
//	                    	double[] weights = new double[numResults - 1]
//	                		double[] values = new double[numResults - 1]
//	                		double sumWeights = 0
//	                		boolean useExactValue = false
//		                	for (int k = 0; k < neighbourIDs.length; k++) {
//		                    	int j = neighbourIDs[k]
//		                    	values[k] = nodalFunc[j].getValue(x, y)
//		                    	double distance = distances[k]
//		                    	double weight = weightFunction.getPartialWeight(distance, dMax)
//			                    weights[k] = weight
//			                    sumWeights += weight
//		                    }
//	                    	z = 0
//	                    	for (int k = 0; k < neighbourIDs.length; k++) {
//	                    		z += values[k] * weights[k] / sumWeights
//	                    	}
//	                    	double errorMetric = Math.abs(z - ptZ)
//	                    	if (errorMetric < minError) {
//	                    		minError = errorMetric
//	                    		minErrorP = testP
//	                    	}
//	                    }
//	                    sumP += minErrorP
//	                    N++
//	                    println(minError + "\t" + minErrorP)
//                    } else {
//                    	pluginHost.showFeedback("Error during script execution")
//                    	return
//                    }
//					
//					progress = (int)(100f * i / (numSamples - 1))
//					if (progress > oldProgress) {
//						pluginHost.updateProgress("Optimizing Exponent Value:", progress)
//						oldProgress = progress
//					}
//					// check to see if the user has requested a cancellation
//					if (pluginHost.isRequestForOperationCancelSet()) {
//						pluginHost.showFeedback("Operation cancelled")
//						return
//					}
//				}
//				exponent = sumP / N
//				weightFunction.setP(exponent)
//			}

			
			
//			pluginHost.updateProgress("Please wait...", 0)
//			ArrayList<DoWorkOnRow> tasks = new ArrayList<>();
//			for (int r in 0..(rows - 1)) {
//				double[] data = output.getRowValues(r)
//				tasks.add(new DoWorkOnRow(data, nodata, r, pointsTree))
//			}
//	    
//			ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
//	  	    // the only reason for the getExecutorResults method 
//	  	    // is that Groovy throws a compilation type mis-match
//	  	    // error when compiled statically. I think it's a bug.
//	  	    List<Future<RowNumberAndData>> rowResults = getExecutorResults(executor, tasks);
//	        executor.shutdown();
//	
//	        i = 0
//	        oldProgress = -1
//		    for (Future<RowNumberAndData> result : rowResults) {
//	    		RowNumberAndData data = result.get()
//	    		row = data.getRow()
//	    		output.setRowValues(row, data.getData())
//	        	i++
//				// update progress bar
//				progress = (int)(100f * i / rows)
//				if (progress > oldProgress) {
//					pluginHost.updateProgress("Progress", progress)
//					oldProgress = progress
//					// check to see if the user has requested a cancellation
//					if (pluginHost.isRequestForOperationCancelSet()) {
//						pluginHost.showFeedback("Operation cancelled")
//						return
//					}
//				}
//		    }
	    
			// perform the interpolation
			double weight, distance, dMax
			oldProgress = -1;
			for (row = 0; row < rows; row++) {
				y = output.getYCoordinateFromRow(row);
				for (col = 0; col < cols; col++) {
                	x = output.getXCoordinateFromColumn(col);
                    double[] entry = new double[2]
                    entry[0] = x
                    entry[1] = y
                    
                    results = pointsTree.neighborsWithinRange(entry, maxSearchDist)
                    dMax = maxSearchDist
                    int numResults = results.size()
                    if (numResults < numNeighbours) {
                    	results = pointsTree.nearestNeighbor(entry, numNeighbours, true, useQuadrandSearch)
                    	dMax = results.get(0).distance
                    	numResults = results.size()
                    }
                	if (numResults >= 1) {
                		double[] weights = new double[numResults]
                		double[] values = new double[numResults]
                		double sumWeights = 0
                		boolean useExactValue = false
	                	for (i = 0; i < numResults; i++) {
	                    	int j = results.get(i).value
	                    	values[i] = nodalFunc[j].getValue(x, y)
	                    	distance = Math.sqrt(results.get(i).distance)
	                    	if (distance > 0) {
	                    		weight = weightFunction.getPartialWeight(distance, dMax)
		                    	weights[i] = weight
		                    	sumWeights += weight
	                    	} else {
	                    		z = values[i]
	                    		useExactValue = true
	                    		break
	                    	}
	                    }
	                    if (!useExactValue) {
	                    	z = 0
	                    	for (i = 0; i < numResults; i++) {
	                    		z += values[i] * weights[i] / sumWeights
	                    	}
	                    }
                	} else {
                		z = nodata
                	}
                    
                    output.setValue(row, col, z)
                }
                progress = (int)(100f * row / (rows - 1))
				if (progress > oldProgress) {
					pluginHost.updateProgress("Interpolating:", progress)
					oldProgress = progress
				}
				// check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
            }

			output.addMetadataEntry("Created by the "
	                    + descriptiveName + " tool.")
	        output.addMetadataEntry("Created on " + new Date())
			
			if (weightType == 0) {
				output.addMetadataEntry("Weight type = " + "shepard")
	        } else if (weightType == 1) {
	        	output.addMetadataEntry("Weight type = " + "franke")
	        } else if (weightType == 2) {
	        	output.addMetadataEntry("Weight type = " + "lindsay")
	        }
	        output.addMetadataEntry("Nodal function = " + nodalFunction)
			output.addMetadataEntry("Exponent parameter = " + exponent)
			output.addMetadataEntry("Search distance = " + maxSearchDist)
			output.addMetadataEntry("Minimum number of points = " + numNeighbours)
			output.addMetadataEntry("Quadrant-based search = " + useQuadrandSearch)
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

	@CompileStatic
    private double[] solveQuadratic(double centerX, double centerY, double centerZ, 
       double[] x, double[] y, double[] z) {
       	int numCoefficients = 5
       	int numSamples = x.length
   		if (numSamples < numCoefficients) {
   			return null
   		}
       	RealMatrix m = new Array2DRowRealMatrix(numSamples, numCoefficients)
	    for (int i = 0; i < numSamples; i++) {
	    	m.addToEntry(i, 0, (x[i] - centerX))
	        m.addToEntry(i, 1, (y[i] - centerY))
	        m.addToEntry(i, 2, (x[i] - centerX) * (x[i] - centerX))
	        m.addToEntry(i, 3, (y[i] - centerY) * (y[i] - centerY))
	        m.addToEntry(i, 4, (x[i] - centerX) * (y[i] - centerY))
	    }

	    //DecompositionSolver solver = new SingularValueDecomposition(m).getSolver();
	    DecompositionSolver solver = new QRDecomposition(m).getSolver();
	    double[] adjustedZ = new double[numSamples]
	    for (int i = 0; i < numSamples; i++) {
	    	adjustedZ[i] = z[i] - centerZ
	    }
	    RealVector values = new ArrayRealVector(adjustedZ)
	    RealVector solution = solver.solve(values);
        double[] coefficients = new double[numCoefficients + 1];
        coefficients[0] = centerZ
        for (int i = 0; i < numCoefficients; i++) {
            coefficients[i + 1] = solution.getEntry(i);
        }
        return coefficients;
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
	class WeightFunction {
		private int weightType = 0
		private double p = 2
		
		WeightFunction(int weightType, double p) {
			this.weightType = weightType
			this.p = p
		}

		double getP() {
			return p
		}

		void setP(double p) {
			this.p = p
		}

		double getPartialWeight(double distance, double maxDist) {
			switch (weightType) {
				case 0: // Shepard's
					return 1 / Math.pow(distance, p)
				case 1:
					if (distance < maxDist) {
                		return Math.pow((maxDist - distance) / (maxDist * distance), p)
                	} else {
                		return 0
                	}
                case 2:
                	return Math.pow((maxDist - distance) / maxDist, 4 * p)
			}
			return 0 // should never reach here
		}
	}
	
	@CompileStatic
    class NodalFunction {
    	private double centerX
    	private double centerY
    	private byte functionType = 0
		private double[] p
		
    	NodalFunction(double x, double y, byte functionType, double[] parameters) {
    		this.centerX = x
    		this.centerY = y
    		this.functionType = functionType // 0 = constant, 1 = quadratic
    		this.p = new double[parameters.length]
    		System.arraycopy(parameters, 0, p, 0, parameters.length) 
    		
    	}

    	double getValue(double x, double y) {
    		if (functionType == 0) { // constant
    			return p[0]
    		} else { // quadratic
    			// z = centerZ + p1 * dX + p2 * dY + p3 * dX^2 + p4 * dY^2 + p5 * dX * dY
    			double dX = x - centerX
    			double dY = y - centerY
    			return (p[0] + p[1] * dX + p[2] * dY + p[3] * dX * dX + p[4] * dY * dY + p[5] * dX * dY)
    		}
    	}
    }

//	public List<Future<RowNumberAndData>> getExecutorResults(ExecutorService executor, ArrayList<DoWorkOnRow> tasks) {
//    	List<Future<RowNumberAndData>> results = executor.invokeAll(tasks);
//		return results
//    }
//
//	//@CompileStatic
//    class DoWorkOnRow implements Callable<RowNumberAndData> {
//		private double[] data
//		private int row
//		private double nodata
//		private KdTree<Integer> pointsTree
//		
//	    DoWorkOnRow(double[] data, double nodata, int row, KdTree<Integer> pointsTree) {
//        	this.data = data
//        	this.row = row
//            this.nodata = nodata
//            this.pointsTree = pointsTree
//       	}
//        	
//        @Override
//	    public RowNumberAndData call() {
//	    	// check to see if the user has requested a cancellation
//			if (pluginHost.isRequestForOperationCancelSet()) {
//				return null
//			}
//			List<KdTree.Entry<Integer>> results
//			try {
//	    	int cols = data.length
//	       	double[] retData = new double[cols]
//   	        double distance, dMax, weight, x, y, z
//   	        int i
//
//	        for (int col in 0..(cols - 1)) {
//            	x = output.getXCoordinateFromColumn(col);
//                y = output.getYCoordinateFromRow(row);
//                double[] entry = new double[2]
//                entry[0] = x
//                entry[1] = y
//                
//                results = pointsTree.neighborsWithinRange(entry, maxSearchDist)
//                dMax = maxSearchDist
//                int numResults = results.size()
//                if (numResults < numNeighbours) {
//                	results = pointsTree.nearestNeighbor(entry, numNeighbours, true, useQuadrandSearch)
//                	dMax = results.get(0).distance
//                	numResults = results.size()
//                }
//            	if (numResults >= 1) {
//            		double[] weights = new double[numResults]
//            		double[] values = new double[numResults]
//            		double sumWeights = 0
//            		boolean useExactValue = false
//                	for (i = 0; i < numResults; i++) {
//                    	int j = results.get(i).value
//                    	values[i] = nodalFunc[j].getValue(x, y)
//                    	distance = Math.sqrt(results.get(i).distance)
//                    	if (distance > 0) {
//                    		weight = weightFunction.getPartialWeight(distance, dMax)
//	                    	weights[i] = weight
//	                    	sumWeights += weight
//                    	} else {
//                    		z = values[i]
//                    		useExactValue = true
//                    		break
//                    	}
//                    }
//                    if (!useExactValue) {
//                    	z = 0
//                    	for (i = 0; i < numResults; i++) {
//                    		z += values[i] * weights[i] / sumWeights
//                    	}
//                    }
//            	} else {
//            		z = nodata
//            	}
//                
//                retData[col] = z
//
//                
//            }
//
//            // check to see if the user has requested a cancellation
//			if (pluginHost.isRequestForOperationCancelSet()) {
//				pluginHost.showFeedback("Operation cancelled")
//				return
//			}
//            
//			def ret = new RowNumberAndData(row, retData)
//			return ret
//        } catch (Exception e) {
//	    	pluginHost.showFeedback("An error has occurred during operation. See log file for details.")
//	        pluginHost.logException("Error in " + descriptiveName, e) 
//	        return null
//	    }
//	    }
//    }
//
//	@CompileStatic
//    class RowNumberAndData {
//		private int row
//		private double[] data
//		private int numDataEntries = 0
//		
//    	RowNumberAndData(int row, double[] data) {
//    		this.row = row
//    		this.data = data
//    		this.numDataEntries = data.length
//    	}
//
//    	public int getRow() {
//    		return row
//    	}
//
//    	public void setRow(int row) {
//    		this.row = row
//    	}
//
//    	public double[] getData() {
//    		return data
//    	}
//
//		public void setData(double[] data) {
//			this.data = data
//		}
//		
//    	public double getDataAt(int n) {
//    		if (n < numDataEntries) {
//    			return data[n]
//    		} else {
//    			return null
//    		}
//    	}
//    }
}

if (args == null) {
    pluginHost.showFeedback("Plugin arguments not set.")
} else {
    def f = new InterpolationIDW(pluginHost, args, descriptiveName)
}
