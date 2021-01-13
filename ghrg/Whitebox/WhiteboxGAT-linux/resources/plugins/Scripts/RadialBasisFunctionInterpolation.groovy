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
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.ShapeFile
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
import java.util.Arrays
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D
import org.apache.commons.math3.geometry.euclidean.threed.Plane
import org.apache.commons.math3.geometry.euclidean.twod.Vector2D
import org.apache.commons.math3.linear.SingularValueDecomposition
import org.apache.commons.math3.linear.Array2DRowRealMatrix
import org.apache.commons.math3.linear.ArrayRealVector
import org.apache.commons.math3.linear.RealVector
import org.apache.commons.math3.linear.RealMatrix
import org.apache.commons.math3.linear.DecompositionSolver
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
//def name = "RadialBasisFunctionInterpolation"
//def descriptiveName = "Radial Basis Function Interpolation"
//def description = "Performs a radial basis function (RBF) interpolation"
//def toolboxes = ["Interpolation"]

public class RadialBasisFunctionInterpolation implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
    private String basisFunctionType = ""
	
    public RadialBasisFunctionInterpolation(WhiteboxPluginHost pluginHost, 
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
            sd.setHelpFile("RadialBasisFunctionInterpolation")
		
            // Specifying the source file allows the 'view code' 
            // button on the tool dialog to be displayed.
            def pathSep = File.separator
            def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "RadialBasisFunctionInterpolation.groovy"
            sd.setSourceFile(scriptFile)
			
            // add some components to the dialog
        	//DialogFile dfIn = sd.addDialogFile("Input file", "Input Vector File:", "open", "Vector Files (*.shp), SHP", true, false)
            DialogFieldSelector dfs = sd.addDialogFieldSelector("Input file and height field.", "Input Height Field:", false)
            DialogCheckBox dcb = sd.addDialogCheckBox("Use z-values", "Use z-values", false)
            dcb.setVisible(false)
            sd.addDialogFile("Output file", "Output Raster File:", "saveAs", "Whitebox Raster Files (*.dep), DEP", true, false)
            sd.addDialogDataInput("Output raster cell size.", "Cell Size (optional):", "", true, false)
            sd.addDialogFile("Input base file", "Base Raster File (optional):", "open", "Whitebox Raster Files (*.dep), DEP", true, true)
            sd.addDialogComboBox("Select a radial basis function type", "Radial Basis Function Type:", ["Multiquadric", "Thin Plate Spline", "Natural Cubic Spline", "Gaussian", "Inverse Multiquadric"], 0)
			DialogDataInput ddi = sd.addDialogDataInput("Smoothing parameter.", "Smoothing Parameter (R2 optional):", "0.1", true, false)

            sd.addDialogLabel("<html><b>Neighbourhood Parameters</b></html>")
            sd.addDialogDataInput("Search distance.", "<html>Search Distance (<i>d<sub>max</sub> ):</html>", "", true, false)
            sd.addDialogDataInput("Minimum Number of neighbours in subsample.", "Number of Points:", "", true, false)
            sd.addDialogCheckBox("Use quadrant-based search", "Use quadrant-based search", false)
            
//            def listener = { evt -> if (evt.getPropertyName().equals("value")) { 
//            		String value = dfs.getValue()
//            		if (value != null && !value.isEmpty()) {
//            			value = value.trim()
//            			String[] strArray = dfs.getValue().split(";")
//            			String fileName = strArray[0]
//            			File file = new File(fileName)
//            			if (file.exists()) {
//	            			ShapeFile shapefile = new ShapeFile(fileName)
//		            		if (shapefile.getShapeType().getDimension() == ShapeTypeDimension.Z) {
//		            			dcb.setVisible(true)
//		            		} else {
//		            			dcb.setVisible(false)
//		            		}
//		            		double north, south, east, west
//		            		north = shapefile.getyMax()
//					        south = shapefile.getyMin()
//					        east = shapefile.getxMax()
//					        west = shapefile.getxMin()
//					        int numSamples = shapefile.getNumberOfRecords()
//					        double diagDist = Math.sqrt((north - south) * (north - south) + (east - west) * (east - west))
//							double R2 = (diagDist * diagDist) / (25 * numSamples)
//							DecimalFormat df = new DecimalFormat("0.000");
//							String output = df.format(R2)
//							ddi.setValue(output)
//		            	} else {
//		            		if (dcb.isVisible()) {
//		            			dcb.setVisible(false)
//		            		}
//		            	}
//            		}
//            	} 
//            } as PropertyChangeListener
//            dfs.addPropertyChangeListener(listener)
            
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
        	int progress, oldProgress, rows, cols, row, col
        	double x, y, z
        	double north, south, east, west
        	double cellSize = 1.0
        	double nodata = -32768
        	double[] weights
        	double minVal, maxVal, range
			List<KdTree.Entry<Integer>> results
            ArrayList<Double> xList = new ArrayList<>()
            ArrayList<Double> yList = new ArrayList<>()
            ArrayList<Double> zList = new ArrayList<>()
            
			if (args.length != 10) {
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
	        basisFunctionType = args[5].toLowerCase()
	        double R2 = 1.0;
			if (!args[6].toLowerCase().trim().equals("not specified") 
			  && !args[6].isEmpty()) {
				R2 = Double.parseDouble(args[6])
			}
			double maxSearchDist = Double.parseDouble(args[7])
			if (maxSearchDist < 0) {
				maxSearchDist = 0
			}
			int numNeighbours = Integer.parseInt(args[8])
			if (numNeighbours < 0) {
				numNeighbours = 0
			}
			boolean useQuadrandSearch = Boolean.parseBoolean(args[9])

			
            ShapeFile input = new ShapeFile(inputFile)
			AttributeTable table = new AttributeTable(inputFile.replace(".shp", ".dbf"))
			ShapeType shapeType = input.getShapeType()
            if (shapeType.getDimension() != ShapeTypeDimension.Z && useZValues) {
            	useZValues = false
            }
			int heightField = -1
            if (inputData.length == 2 && !inputData[1].trim().isEmpty()) {
            	String heightFieldName = inputData[1].trim()
            	String[] fieldNames = input.getAttributeTableFields()
            	for (int i = 0; i < fieldNames.length; i++) {
            		if (fieldNames[i].trim().equals(heightFieldName)) {
            			heightField = i
            			break
            		}
            	}
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
					z = (Double)(recData[heightField])
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
			minVal = Double.POSITIVE_INFINITY
			maxVal = Double.NEGATIVE_INFINITY
			KdTree<Integer> pointsTree = new KdTree.SqrEuclid<Integer>(2, new Integer(numSamples));
			for (i = 0; i < numSamples; i++) {
				double[] entry = new double[2]
				entry[0] = xList.get(i)
				entry[1] = yList.get(i);
				pointsTree.addPoint(entry, i);
				if (zList.get(i) < minVal) { minVal = zList.get(i) }
				if (zList.get(i) > maxVal) { maxVal = zList.get(i) }
				
			}
			range = maxVal - minVal
			
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

			double[] xArray = new double[numSamples]
			double[] yArray = new double[numSamples]
			double[] zArray = new double[numSamples]

			for (i = 0; i < numSamples; i++) {
				xArray[i] = xList.get(i)
				yArray[i] = yList.get(i)
				zArray[i] = zList.get(i) //- minVal) / range
			}

			// calculate the weights
			boolean globalWeights = true
			weights = new double[numSamples]
			if (!globalWeights) {
//			int targetNumNeighbours = Math.max(3, numNeighbours)
//			oldProgress = -1
//			for (i = 0; i < numSamples; i++) {
//				x = xArray[i]
//                y = yArray[i]
//                
//				double[] entry = new double[2]
//                entry[0] = x
//                entry[1] = y
//                
//                results = pointsTree.neighborsWithinRange(entry, maxSearchDist)
//				int numResults = results.size()
//				if (numResults < targetNumNeighbours) {
//                	results = pointsTree.nearestNeighbor(entry, targetNumNeighbours, true, useQuadrandSearch)
//                	numResults = results.size()
//                }
//                if (numResults >= targetNumNeighbours) {
//	                double[] xVals = new double[numResults]
//					double[] yVals = new double[numResults]
//					double[] zVals = new double[numResults]
//					int m = -1
//					
//					for (int j = 0; j < numResults; j++) {
//						KdTree.Entry result = results.get(j)
//						int k = result.value
//						xVals[j] = xArray[k]
//						yVals[j] = yArray[k]
//						zVals[j] = zArray[k]
//						if (result.distance == 0) {
//							m = j
//						}
//					}
//					if (m >= 0) {
//	                	double[] localWeights = getWeights(xVals, yVals, zVals, R2)
//	                	if (localWeights != null) {
//	                		weights[i] = localWeights[m]
//	                	}
//					}
//                }
//				
//				progress = (int)(100f * i / (numSamples - 1))
//				if (progress > oldProgress) {
//					pluginHost.updateProgress("Calculating Weights:", progress)
//					oldProgress = progress
//				}
//				// check to see if the user has requested a cancellation
//				if (pluginHost.isRequestForOperationCancelSet()) {
//					pluginHost.showFeedback("Operation cancelled")
//					return
//				}
//			}
			} else {
				pluginHost.updateProgress("Calculating weights:", 0)
				weights = getWeights(xArray, yArray, zArray, R2)
			}

//			for (i = 0; i < numSamples; i++) {
//				println weights[i]
//			}

			int[] neighbourIDs = new int[numNeighbours]
			double[] xVals = new double[numNeighbours]
			double[] yVals = new double[numNeighbours]
			double[] zVals = new double[numNeighbours]
			int sameWeightCount = 0		
			oldProgress = -1;
            for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                    x = output.getXCoordinateFromColumn(col);
                    y = output.getYCoordinateFromRow(row);
                    if (!globalWeights) {
	                    double[] entry = new double[2]
	                    entry[0] = x
	                    entry[1] = y
	                    
		                results = pointsTree.neighborsWithinRange(entry, maxSearchDist)
						int numResults = results.size()
//						if (numResults < numNeighbours) {
//		                	results = pointsTree.nearestNeighbor(entry, numNeighbours, true, useQuadrandSearch)
//		                	numResults = results.size()
//		                }
						if (numResults >= numNeighbours) {
			                // have the neighbours changed since last time?
			                boolean neighboursChanged = false
			                if (numResults != neighbourIDs.length) {
			                	neighboursChanged = true
			                	neighbourIDs = new int[numResults]
			                	for (int j = 0; j < numResults; j++) {
									neighbourIDs[j] = results.get(j).value
								}
								Arrays.sort(neighbourIDs)
			                } else {
			                	int[] newNeighbours = new int[numResults]
			                	for (int j = 0; j < numResults; j++) {
									newNeighbours[j] = results.get(j).value
								}
								Arrays.sort(newNeighbours)
								for (int j = 0; j < numResults; j++) {
									if (newNeighbours[j] != neighbourIDs[j]) {
										System.arraycopy(newNeighbours, 0, neighbourIDs, 0, newNeighbours.length)
										neighboursChanged = true
										break
									}
								}
			                }
		                
		                	if (neighboursChanged) {
				                xVals = new double[numResults]
								yVals = new double[numResults]
								zVals = new double[numResults]
								//double[] localWeights = new double[numResults]
								for (int j = 0; j < numResults; j++) {
									int k = neighbourIDs[j] //results.get(j).value
									//neighbourIDs[j] = k
									xVals[j] = xArray[k]
									yVals[j] = yArray[k]
									zVals[j] = zArray[k]
									//localWeights[j] = weights[k]
								}
								weights = getWeights(xVals, yVals, zVals, R2)
		                	} else {
		                		sameWeightCount++
		                	}
							z = interpolate(x, y, xVals, yVals, weights, R2)
		                } else {
		                	z = nodata
		                }
                    } else {
                    	z = interpolate(x, y, xArray, yArray, weights, R2)
                    }
                    output.setValue(row, col, z)
                }
                progress = (int)(100f * row / rows)
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

			println("number of same weights = " + sameWeightCount)
			
			output.addMetadataEntry("Created by the "
	                    + descriptiveName + " tool.")
	        output.addMetadataEntry("Created on " + new Date())
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

//	@CompileStatic
//    private int getQuadrant(double val1, double val2) {
//		if (val1 > 0.0 && val2 > 0.0) {
//			return 0
//		} else if (val1 > 0.0 && val2 <= 0.0) {
//			return 1
//		} else if (val1 <= 0.0 && val2 <= 0.0) {
//			return 2
//		} else {
//			return 3
//		}
//	}


	// based on algoirthm by Adrien Herubel on http://www.altdevblogaday.com/2011/09/04/interpolation-using-radial-basis-functions/
	@CompileStatic
    double RBFGauss(double d, double r) { return Math.exp(-(d / r)); }

	@CompileStatic
    double RBFInvMultiquadradic(double d, double r) { return 1 / Math.sqrt(d + r); }

//	@CompileStatic
//    double RBFMultiquadradic(double d, double r) { return Math.sqrt(d + r); }

	@CompileStatic
    double RBFMultiquadradic(double d, double r) { return Math.sqrt(1 + (d / r) * (d / r)); }

	final double log10 = Math.log(10)
	@CompileStatic
    double RBFThinPlateSpline(double d, double r) { return (d + r) * (Math.log(d + r) / log10); }

	@CompileStatic
    double RBFNaturalCubicSpline(double d, double r) { return Math.pow((d + r), 1.5); }

	@CompileStatic
    double[] getWeights(double[] x, double[] y, double[] v, double R2) {
	    int numSamples = x.length
	    if (y.length != numSamples || v.length != numSamples) {
	    	return null
	    }
	    RealMatrix m = new Array2DRowRealMatrix(numSamples, numSamples)
	    double norm
	    switch (basisFunctionType) {
	    	case ("multiquadric"):
			    for (int i = 0; i < numSamples; i++) {
			        for (int j = 0; j < numSamples; j++) {
			        	norm = ((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]))
			            m.addToEntry(i, j, RBFMultiquadradic(norm, R2));
			        }
			    }
			    break;
			case ("thin plate spline"):
			    for (int i = 0; i < numSamples; i++) {
			        for (int j = 0; j < numSamples; j++) {
			        	norm = ((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]))
			            m.addToEntry(i, j, RBFThinPlateSpline(norm, R2));
			        }
			    }
			    break;
			case ("inverse multiquadric"):
			    for (int i = 0; i < numSamples; i++) {
			        for (int j = 0; j < numSamples; j++) {
			        	norm = ((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]))
			            m.addToEntry(i, j, RBFInvMultiquadradic(norm, R2));
			        }
			    }
			    break;
			case ("gaussian"):
			    for (int i = 0; i < numSamples; i++) {
			        for (int j = 0; j < numSamples; j++) {
			        	norm = ((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]))
			            m.addToEntry(i, j, RBFGauss(norm, R2));
			        }
			    }
			    break;
			case ("natural cubic spline"):
			    for (int i = 0; i < numSamples; i++) {
			        for (int j = 0; j < numSamples; j++) {
			        	norm = ((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]))
			            m.addToEntry(i, j, RBFNaturalCubicSpline(norm, R2));
			        }
			    }
			    break;
	    }
	    DecompositionSolver solver = new SingularValueDecomposition(m).getSolver();
	    RealVector values = new ArrayRealVector(v)
	    RealVector solution = solver.solve(values);
        double[] weights = new double[numSamples];
        for (int i = 0; i < numSamples; i++) {
            weights[i] = solution.getEntry(i);
        }
	    
		return weights;
	}

	@CompileStatic
	double interpolate(double x0, double y0, double[] x, double[]y, double[] w, double R2) {
	    int numSamples = x.length
	    double res = 0.0;
	    double norm
	    switch (basisFunctionType) {
	    	case ("multiquadric"):
			    for (int i = 0; i < numSamples; i++) {
			    	norm = ((x0 - x[i]) * (x0 - x[i]) + (y0 - y[i]) * (y0 - y[i]))
			        res += w[i] * RBFMultiquadradic(norm, R2);
			    }
			    break;
			case ("thin plate spline"):
			    for (int i = 0; i < numSamples; i++) {
			    	norm = ((x0 - x[i]) * (x0 - x[i]) + (y0 - y[i]) * (y0 - y[i]))
			        res += w[i] * RBFThinPlateSpline(norm, R2);
			    }
			    break;
			case ("inverse multiquadric"):
			    for (int i = 0; i < numSamples; i++) {
			    	norm = ((x0 - x[i]) * (x0 - x[i]) + (y0 - y[i]) * (y0 - y[i]))
			        res += w[i] * RBFInvMultiquadradic(norm, R2);
			    }
			    break;
			case ("gaussian"):
			    for (int i = 0; i < numSamples; i++) {
			    	norm = ((x0 - x[i]) * (x0 - x[i]) + (y0 - y[i]) * (y0 - y[i]))
			        res += w[i] * RBFGauss(norm, R2);
			    }
			    break;
			case ("natural cubic spline"):
			    for (int i = 0; i < numSamples; i++) {
			        norm = ((x0 - x[i]) * (x0 - x[i]) + (y0 - y[i]) * (y0 - y[i]))
			        res += w[i] * RBFNaturalCubicSpline(norm, R2);
			    }
			    break;
	    }
	    
	    return res;
	}

}

if (args == null) {
    pluginHost.showFeedback("Plugin arguments not set.")
} else {
    def f = new RadialBasisFunctionInterpolation(pluginHost, args, descriptiveName)
}
