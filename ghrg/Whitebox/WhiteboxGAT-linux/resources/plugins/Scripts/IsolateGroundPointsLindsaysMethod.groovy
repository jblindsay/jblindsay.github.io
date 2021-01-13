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
 
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.util.concurrent.Future;
import java.util.concurrent.*;
import java.text.DecimalFormat;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.IntStream;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.geospatialfiles.LASReader;
import whitebox.geospatialfiles.LASReader.PointRecord;
import whitebox.structures.KdTree;
import whitebox.ui.plugin_dialog.ScriptDialog;
import whitebox.utilities.StringUtilities;
import whitebox.geospatialfiles.LASReader.VariableLengthRecord;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.*;
import whitebox.geospatialfiles.shapefile.attributes.*;
import whitebox.geospatialfiles.VectorLayerInfo;
import whitebox.utilities.Topology;
import groovy.transform.CompileStatic;
import whitebox.structures.BoundingBox;
import com.vividsolutions.jts.geom.PrecisionModel;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.LinearRing;
import com.vividsolutions.jts.geom.impl.CoordinateArraySequence;
import com.vividsolutions.jts.triangulate.VoronoiDiagramBuilder;
import whitebox.structures.BooleanBitArray1D;

/*
 * This tool can be used to identify points within a point cloud
 * contained within a LAS file that correspond with the ground 
 * surface. The points are then output into a MultiPoint shapefile.
 */
def name = "IsolateGroundPointsLindsaysMethod"
def descriptiveName = "Isolate Ground Points (Lindsay's Method)"
def description = "Isolates points associated with the ground surface in a LiDAR point cloud."
def toolboxes = ["LidarTools"]

public class IsolateGroundPointsLindsaysMethod {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	private InterpolationRecord[] data;
	private double threshold;
    private double searchDist;
    private long numClassifiedPoints = 0;
    private BooleanBitArray1D done;
    private KdTree<Integer> pointsTree;
	private int numPoints;
	private int progress, oldProgress = -1;
	
	public IsolateGroundPointsLindsaysMethod(WhiteboxPluginHost pluginHost, 
		String[] args, def name, def descriptiveName) {
		this.pluginHost = pluginHost
		this.descriptiveName = descriptiveName
			
		if (args.length > 0) {
			execute(args)
		} else {
			// create an ActionListener to handle the return from the dialog
            def ac = new ActionListener() {
                public void actionPerformed(ActionEvent event) {
                    if (event.getActionCommand().equals("ok")) {
                        args = sd.collectParameters()
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
        	};
			// Create a dialog for this tool to collect user-specified
			// tool parameters.
			sd = new ScriptDialog(pluginHost, descriptiveName, ac)	
		
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
			sd.addDialogFile("Input LAS file", "Input LAS File:", "open", "LAS Files (*.las), LAS", true, false)
            sd.addDialogFile("Output file", "Output Vector File:", "close", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogDataInput("Radius (m)", "Radius (m)", "", true, false)
			sd.addDialogDataInput("Pin height (m)", "Pin Height (m)", "", true, false)
			
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
	  		// make sure there are the appropriate number of arguments
		  	if (args.length != 4) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}

			// declare some variables
			double x, y, z, zN, xV, yV, xN, yN, dist, angle, dx, dy;
			int i, j, p, a, intensity;
			double angleN, startingAngle1, startingAngle2
			boolean pointWithinSquare1, pointWithinSquare2		
			PointRecord point;
	    	List<KdTree.Entry<Integer>> results
	    	InterpolationRecord value
	    	double[] entry;
	    	double halfPi = Math.PI / 2;
			double twoPi = Math.PI * 2;

			// Read the input parameters
			String inputFile = args[0];
            String outputFile = args[1];
			double r = Double.parseDouble(args[2]);
			double rSqrd = r * r;
		    double twoR = 2 * r;
		    double rSqrRoot2 = r * Math.sqrt(2);
		    searchDist = r;

		    double pinHeight = Double.parseDouble(args[3]);
			threshold = pinHeight * 1.5;
			
	    	// Create the LAS object
            LASReader las = new LASReader(inputFile);
			numPoints = (int) las.getNumPointRecords();

			// Read the points into the k-dimensional tree.
            pointsTree = new KdTree.SqrEuclid<Integer>(2, new Integer(numPoints));
            data = new InterpolationRecord[numPoints];
            done = new BooleanBitArray1D(numPoints);
			for (a = 0; a < numPoints; a++) {
                point = las.getPointRecord(a);
                if (!point.isPointWithheld()) {
                    x = point.getX();
                    y = point.getY();
                    z = point.getZ();
                    intensity = point.getIntensity();
                    entry = [x, y];
                    pointsTree.addPoint(entry, a);
                    
					data[a] = new InterpolationRecord(x, y, z, intensity, a);
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
	           points that can be touched by a corner of a square prism
	           when approached from the bottom. */
	        BooleanBitArray1D isOnTerrain = new BooleanBitArray1D(numPoints);
            for (a = 0; a < numPoints; a++) {
                x = data[a].x; 
                y = data[a].y; 
                z = data[a].z - pinHeight;
               
                entry = [x, y];
                results = pointsTree.neighborsWithinRange(entry, rSqrRoot2);
				int numResults = results.size()
				
				/* Create a list of the neighbouring points that 
				 *  are lower than the point of interest
				 */
				ArrayList<NeighbouringPoint> coordsList = new ArrayList<>();
				for (p = 0; p < numResults; p++) {
            		int pointNum = results.get(p).value;
        			if (z > data[pointNum].z && pointNum != a) { // only consider points that are lower than the centre point
        				coordsList.add(new NeighbouringPoint(data[pointNum].x, data[pointNum].y, 0, 0));
        			}
            	}

				if  (coordsList.size() < 3) {
					// the prism can touch the point
            		isOnTerrain.setValue(a, true);
            	} else {
					// calculate the dist and angle for the neighbouring points
					for (p = 0; p < coordsList.size(); p++) {
						xN = (coordsList.get(p)).x;
	                    yN = (coordsList.get(p)).y;
	                    dx = xN - x;
						dy = yN - y;
						(coordsList.get(p)).dist = Math.sqrt(dx * dx + dy * dy);
						angle = Math.atan2(dy, dx);
						angle = (angle + twoPi) % twoPi
						(coordsList.get(p)).angle = angle;
					}
					
					// see if there is an angle in which you can fit a 
					// square prism touching the centre point with a 
					// corner of the prism but no other neighbour
					for (p = 0; p < coordsList.size(); p++) {
						dist = coordsList.get(p).dist
						if (dist > 0) {
							if (dist < r) {
								dist = r
							}
							startingAngle1 = coordsList.get(p).angle - Math.acos(r / dist)
							startingAngle1 = (startingAngle1 + twoPi) % twoPi

							startingAngle2 = coordsList.get(p).angle + Math.acos(r / dist) - halfPi
							startingAngle2 = (startingAngle2 + twoPi) % twoPi

							pointWithinSquare1 = false
							pointWithinSquare2 = false
							for (int q = 0; q < coordsList.size(); q++) {
								if (q != p) { // of course the defining neighbouring point will be in its squares but no others are allowed to be
									if (!pointWithinSquare1) {
										angleN = ((coordsList.get(q).angle - startingAngle1) + twoPi) % twoPi;
										if (angleN < halfPi) {
											dist = coordsList.get(q).dist
											if (dist <= r) {
												// the neighbour is within the square
												pointWithinSquare1 = true;
											} else {
												// it could still be within the square; we have to check
												if (Math.sin(angleN) * dist <= r) {
													// the neighbour is within the square
													pointWithinSquare1 = true;
												}
											}
										}
									}
									if (!pointWithinSquare2) {
										angleN = ((coordsList.get(q).angle - startingAngle2) + twoPi) % twoPi;
										if (angleN < halfPi) {
											dist = coordsList.get(q).dist
											if (dist <= r) {
												// the neighbour is within the square
												pointWithinSquare2 = true;
											} else {
												// it could still be within the square; we have to check
												if (Math.sin(angleN) * dist <= r) {
													// the neighbour is within the square
													pointWithinSquare2 = true;
												}
											}
										}
									}

									if (pointWithinSquare1 && pointWithinSquare2) {
										break;
									}
								}
							}
							if (!pointWithinSquare1 || !pointWithinSquare2) {
								isOnTerrain.setValue(a, true);
								break
							}
						} else {
							// there is a point directly beneath our point. It's not on the terrain therefore.
							isOnTerrain.setValue(a, false);
							break
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


            

//            /* visit the neighbourhood around each point hunting for 
//	           points that can be touched by a cylinder of r radius when
//	           approached from the bottom. */
//	        BooleanBitArray1D isOnTerrain = new BooleanBitArray1D(numPoints);
//            for (a = 0; a < numPoints; a++) {
//                x = data[a].x; 
//                y = data[a].y; 
//                z = data[a].z - pinHeight;
//               
//                entry = [x, y];
//                results = pointsTree.neighborsWithinRange(entry, twoR);
//				int numResults = results.size()
//				
//				/* Create a list of the neighbouring points that 
//				 *  are lower than the point of interest
//				 */
//				ArrayList<Coordinate> coordsList = new ArrayList<>();
//				for (p = 0; p < numResults; p++) {
//            		int pointNum = results.get(p).value;
//        			if (z > data[pointNum].z) { // only consider points that are lower
//        				coordsList.add(new Coordinate(data[pointNum].x, data[pointNum].y));
//        			}
//            	}
//
//				if  (coordsList.size() < 2) {
//					// the cylinder can touch the point
//            		isOnTerrain.setValue(a, true);
//            	} else {
//					/* See whether the edge of the cyclinder
//					 *  can touch the point
//					 */
//					boolean[] degreesIntersected = new boolean[360];
//					int numIntersected = 0;
//					for (p = 0; p < coordsList.size(); p++) {
//						xN = coordsList.get(p).x;
//	                    yN = coordsList.get(p).y;
//	                    double dx = xN - x;
//						double dy = yN - y;
//						dist = Math.sqrt(dx * dx + dy * dy);
//						if (dist == 0) {
//							// there is a point directly beneath the point of interest
//							isOnTerrain.setValue(a, false);
//							numIntersected = 360;
//							break;
//						}
//						double alpha = -Math.atan2(dy, dx) + halfPi;
//						if (alpha < 0) alpha += twoPi;
//						alpha = Math.toDegrees(alpha);
//						double beta = Math.acos(dist / twoR);
//						beta = Math.toDegrees(beta);
//						int lower = (int)alpha - beta;
//						int upper = (int)alpha + beta;
//						if (lower >= 0 && upper < 360) {
//							for (i = lower; i <= upper; i++) {
//								if (!degreesIntersected[i]) {
//									degreesIntersected[i] = true;
//									numIntersected++;
//								}
//							}
//						} else if (lower < 0 && upper < 360) {
//							for (i = 0; i <= upper; i++) {
//								if (!degreesIntersected[i]) {
//									degreesIntersected[i] = true;
//									numIntersected++;
//								}
//							}
//							for (i = 360 + lower; i <= 359; i++) {
//								if (!degreesIntersected[i]) {
//									degreesIntersected[i] = true;
//									numIntersected++;
//								}
//							}
//						} else if (lower >= 0 && upper >= 360) {
//							for (i = lower; i < 360; i++) {
//								if (!degreesIntersected[i]) {
//									degreesIntersected[i] = true;
//									numIntersected++;
//								}
//							}
//							for (i = 0; i <= upper - 360; i++) {
//								if (!degreesIntersected[i]) {
//									degreesIntersected[i] = true;
//									numIntersected++;
//								}
//							}
//						}
//						if (numIntersected == 360) {
//							break;
//						}
//					}
//						
//					if (numIntersected < 360) {
//						isOnTerrain.setValue(a, true);
//						//pointCount++;
//					} else {
//						data[a].setClassValue(-99, true);
//					}
//            	}
//				progress = (int) (100f * (a + 1) / numPoints);
//                if (progress != oldProgress) {
//                    oldProgress = progress;
//                    pluginHost.updateProgress("Finding ground points:", progress);
//                    if (pluginHost.isRequestForOperationCancelSet()) {
//                        pluginHost.showFeedback("Operation cancelled")
//						return
//                    }
//                }
//            }

			
//            // perform the segmentation
//            int currentClass = 0;
//			long oldNumClassifiedPoints = 0;
//			List<Long> histo = new ArrayList<>();
//			boolean didSomething = true;
//            while (numClassifiedPoints < numPoints || didSomething) {
//            	didSomething = false;
//            	
//                // find the lowest unclassified point
//                int startingPoint = -1;
//                double lowestPointZ = Double.POSITIVE_INFINITY;
//                for (a = 0; a < numPoints; a++) {
//                    if (data[a].classValue == -1 && data[a].z < lowestPointZ) {
//                        lowestPointZ = data[a].z;
//                        startingPoint = a;
//                    }
//                }
//                if (startingPoint == -1) {
//                    break;
//                }
//
//                currentClass++;
//
//                List<Integer> seeds = new ArrayList<>();
//                seeds.add(startingPoint);
//                boolean flag = true;
//
//                while (flag) {
//                    flag = false;
//                    for (Integer s : seeds) {
//                        if (!done.getValue(s)) {
//                            data[s].setClassValue(currentClass);
//                            scanNeighbours(s);
//                        }
//                    }
//                    seeds.clear();
//                    if (seedPoints.size() > 0) {
//                        flag = true;
//                        seedPoints.forEach() { it ->
//                            if (!done.getValue((int)it)) {
//                                seeds.add((int)it);
//                            }
//                        };
//                        seedPoints.clear();
//                    }
//
//                    startingPoint = -1;
//
//                    if (pluginHost.isRequestForOperationCancelSet()) {
//                        pluginHost.showFeedback("Operation cancelled")
//						return
//                    }
//                }
//                //int numPointsInClass = (int)(numClassifiedPoints - oldNumClassifiedPoints);
//                //pluginHost.showFeedback("Class: $currentClass Num Points: $numPointsInClass");
//                //histo.add(numClassifiedPoints - oldNumClassifiedPoints);
//                if (oldNumClassifiedPoints < numClassifiedPoints) {
//                	didSomething = true;
//                }
//                oldNumClassifiedPoints = numClassifiedPoints;
//
//                if (pluginHost.isRequestForOperationCancelSet()) {
//                    pluginHost.showFeedback("Operation cancelled")
//					return
//                }
//            }

            // output the data
			DBFField[] fields = new DBFField[1];

            fields[0] = new DBFField();
            fields[0].setName("ID");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setFieldLength(10);
            fields[0].setDecimalCount(0);
            
			ShapeFile output = new ShapeFile(outputFile, ShapeType.MULTIPOINTZ, fields);

			// how many points are being output?
			int numOutputPoints = 0;
			for (a = 0; a < numPoints; a++) {
				if (isOnTerrain.getValue(a)) {// && histo.get(data[a].classValue - 1) > 1000) {
					numOutputPoints++;
				}
            }
            double[][] xyData = new double[numOutputPoints][2]
            double[] zData = new double[numOutputPoints]
			double[] mData = new double[numOutputPoints]
			int q = 0
			for (a = 0; a < numPoints; a++) {
				if (isOnTerrain.getValue(a)) {// && histo.get(data[a].classValue - 1) > 1000) {
					xyData[q][0] = data[a].x;
					xyData[q][1] = data[a].y;
					zData[q] = data[a].z;
					mData[q] = data[a].intensity;
					q++;
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
			
            output.write()

            //pluginHost.returnData(outputFile);

            // display the output image
			String paletteDirectory = pluginHost.getResourcesDirectory() + "palettes" + File.separator;
			VectorLayerInfo vli = new VectorLayerInfo(outputFile, paletteDirectory, 255i, -1);
			vli.setPaletteFile(paletteDirectory + "spectrum.pal");
			vli.setFilledWithOneColour(false);
			vli.setFillAttribute("Feature Z Value");
			vli.setPaletteScaled(true);
			vli.setDisplayMinValue(las.getMinZ());
			vli.setDisplayMaxValue(las.getMaxZ());
			vli.setMarkerSize(2.5f);
			vli.setRecordsColourData();
			pluginHost.returnData(vli);
      
	  	} catch (OutOfMemoryError oe) {
            pluginHost.showFeedback("An out-of-memory error has occurred during operation.")
	    } catch (Exception e) {
	        pluginHost.showFeedback("$e"); //"An error has occurred during operation. See log file for details.")
	        pluginHost.logException("Error in " + descriptiveName, e)
        } finally {
        	// reset the progress bar
        	pluginHost.updateProgress("", 0)
        }
	}

	private long depth = 0;
    private final long maxDepth = 1000;
    private final List<Integer> seedPoints = new ArrayList<>();
	
//	@CompileStatic
//    private void scanNeighbours(int refPointNum) {
//    	if (pluginHost.isRequestForOperationCancelSet()) {
//            return
//        }
//        depth++;
//        if (depth > maxDepth) {
//            if (seedPoints.size() < 80000) {
//                seedPoints.add(refPointNum);
//            }
//            depth--;
//            return;
//        }
//
//        if (done.getValue(refPointNum)) {
//            depth--;
//            return;
//        }
//
//        int classValue = data[refPointNum].classValue;
//		double x = data[refPointNum].x;
//		double y = data[refPointNum].y;
//        double[] entry = [x, y];
//        List<KdTree.Entry<Integer>> results = pointsTree.neighborsWithinRange(entry, searchDist);
//        for (int i = 0; i < results.size(); i++) {
//            int pointNum = results.get(i).value;
//            if (pointNum != data[refPointNum].index) {
//            	if (data[pointNum].classValue == -1) {
//                	if (Math.abs(data[pointNum].z - data[refPointNum].z) <= threshold) {
//                        data[pointNum].setClassValue(classValue);
//                        scanNeighbours(pointNum);
//                    }
//                }
//            }
//        }
//        done.setValue(refPointNum, true);
//        depth--;
//    }

	@CompileStatic
	class NeighbouringPoint implements Comparable<NeighbouringPoint> {
		double x;
        double y;
        double dist;
        double angle;

        public NeighbouringPoint(double x, double y, double dist, double angle) {
        	this.x = x
        	this.y = y
        	this.dist = dist
        	this.angle = angle
        }
        
        @Override
        public int compareTo(NeighbouringPoint other) {
            final int BEFORE = -1;
            final int EQUAL = 0;
            final int AFTER = 1;

            if (this.angle < other.angle) {
                return BEFORE;
            } else if (this.angle > other.angle) {
                return AFTER;
            }

            if (this.dist < other.dist) {
                return BEFORE;
            } else if (this.dist > other.dist) {
                return AFTER;
            }
            
            return EQUAL;
        }

        @Override
        public String toString() {
        	return "X: " + String.valueOf(this.x) + " Y: " + String.valueOf(this.y) + " Dist: " + String.valueOf(this.dist) + " Angle: " + String.valueOf(this.angle) + "\n";
    	}
	}

	@CompileStatic
	class InterpolationRecord {
        
        double x;
        double y;
        double z;
        int intensity;
        int index;
        int classValue = -1;
        
        InterpolationRecord(double x, double y, double z, int intensity, int index) {
            this.x = x;
            this.y = y;
            this.z = z;
            this.index = index;
            this.intensity = intensity;
        }

        public void setClassValue(int value) {
            if (classValue != value) {
                this.classValue = value;
                numClassifiedPoints++;

                progress = (int) (100f * numClassifiedPoints / numPoints);
                if (progress != oldProgress) {
                	oldProgress = progress;
                    pluginHost.updateProgress("Segmenting points...", progress);
                    if (pluginHost.isRequestForOperationCancelSet()) {
                        pluginHost.showFeedback("Operation cancelled")
						return
                    }
                }
            }
        }

        public void setClassValue(int value, boolean suppressProgressUpdate) {
            if (classValue != value) {
                this.classValue = value;
                numClassifiedPoints++;
                
				if (!suppressProgressUpdate) {
	                progress = (int) (100f * numClassifiedPoints / numPoints);
	                if (progress != oldProgress) {
	                	oldProgress = progress;
	                    pluginHost.updateProgress("Segmenting points...", progress);
	                    if (pluginHost.isRequestForOperationCancelSet()) {
	                        pluginHost.showFeedback("Operation cancelled")
							return
	                    }
	                }
				}
            }
        }
    }
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def myClass = new IsolateGroundPointsLindsaysMethod(pluginHost, args, name, descriptiveName)
}
