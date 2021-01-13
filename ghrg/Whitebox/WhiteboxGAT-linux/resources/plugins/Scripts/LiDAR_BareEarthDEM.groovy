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
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.LASReader
import whitebox.geospatialfiles.LASReader.PointRecord
import whitebox.geospatialfiles.LASReader.PointRecColours
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.WhiteboxRasterBase.DataScale
import whitebox.structures.BoundingBox
import whitebox.structures.KdTree
import whitebox.structures.BooleanBitArray1D
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import whitebox.utilities.FileUtilities;
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "LiDAR_BareEarthDEM"
def descriptiveName = "Bare-Earth DEM (LiDAR)"
def description = "Interpolates LAS files into a bare-Earth digital elevation model (DEM)."
def toolboxes = ["LidarTools"]

public class LiDAR_BareEarthDEM implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	private AtomicInteger numSolvedTiles = new AtomicInteger(0)
	
	public LiDAR_BareEarthDEM(WhiteboxPluginHost pluginHost, 
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
			sd.addDialogMultiFile("Select the input LAS files", "Input LAS Files:", "LAS Files (*.las), LAS")
			sd.addDialogDataInput("Output File Suffix (e.g. bare earth)", "Output File Suffix (e.g. bare earth)", "bare earth", false, false)
			sd.addDialogDataInput("IDW Exponent", "IDW Exponent", "2", true, false)
			sd.addDialogDataInput("Max Search Distance (m)", "Max Search Distance (m)", "", true, false)
			sd.addDialogDataInput("Grid Resolution (m)", "Grid Resolution:", "", true, false)
			sd.addDialogDataInput("Threshold in the slope between points to define an off-terrain point.", "Inter-point Slope Threshold:", "30.0", true, false)
			sd.addDialogDataInput("Max Scan Angle Deviation (optional)", "Max Scan Angle Deviation (optional)", "", true, true)
			sd.addDialogDataInput("Minimum elevation (used to exclude drop-out points; optional).", "Minimum Elevation:", "", true, true)
			
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
		long start = System.currentTimeMillis()  
	  try {
	  	if (args.length != 8) {
			pluginHost.showFeedback("Incorrect number of arguments given to tool.")
			return
		}
		
		// read the input parameters
		final String inputFileString = args[0]
		String suffix = ""
		if (args[1] != null) {
			suffix = args[1].trim();
		}
        double weight = Double.parseDouble(args[2]);
		double maxDist = Double.parseDouble(args[3]);
		double resolution = Double.parseDouble(args[4]);
        double maxSlope = Double.parseDouble(args[5]);
        double maxScanAngleDeviation = 1000.0d
        if (args[6] != null && !args[6].isEmpty()) {
        	if (!args[6].toLowerCase().equals("not specified")) {
        		maxScanAngleDeviation= Double.parseDouble(args[6])
        	}
        }
        if (maxScanAngleDeviation < 1.0) { maxScanAngleDeviation = 1.0; }
        double minHeight = Double.MIN_VALUE;
        if (args[7] != null && !args[7].isEmpty()) {
        	if (!args[7].toLowerCase().equals("not specified")) {
        		minHeight = Double.parseDouble(args[7]);
        	}
        }
        
		String[] inputFiles = inputFileString.split(";")

		// check for empty entries in the inputFiles array
		int i = 0;
		for (String inputFile : inputFiles) {
			if (!inputFile.isEmpty()) {
				i++
			}
		}

		if (i != inputFiles.length) {
			// there are empty entries in the inputFiles array
			// remove them.
			ArrayList<String> temp = new ArrayList<>();
			for (String inputFile : inputFiles) {
				if (!inputFile.isEmpty()) {
					temp.add(inputFile)
				}
			}
			inputFiles = new String[i];
    		temp.toArray(inputFiles);
		}
		
		//int rows, cols
		final int numFiles = inputFiles.length
		BoundingBox[] bb = new BoundingBox[numFiles]
		i = 0
		for (String inputFile : inputFiles) {
			LASReader las = new LASReader(inputFile)
			bb[i] = new BoundingBox(las.getMinX(), las.getMinY(), las.getMaxX(), las.getMaxY());
			i++
		}


		pluginHost.updateProgress("Please wait...", 0)
		ArrayList<DoWork> tasks = new ArrayList<>();
		for (i = 0; i < numFiles; i++) {
			tasks.add(new DoWork(i, inputFiles, suffix, 
	      		 bb, resolution, maxSlope, maxScanAngleDeviation,
	      		 weight, maxDist, minHeight))
		}

		ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
  	    // the only reason for the getExecutorResults method 
  	    // is that Groovy throws a compilation type mis-match
  	    // error when compiled statically. I think it's a bug.
  	    List<Future<Boolean>> results = getExecutorResults(executor, tasks); //executor.invokeAll(tasks);
        executor.shutdown();

        i = 0
        int progress = 0
        int oldProgress = -1
        int numSuccessfulInterpolations = 0
	    for (Future<Boolean> result : results) {
    		Boolean data = result.get()
    		if (data) { numSuccessfulInterpolations++ }
        	i++
			// update progress bar
			progress = (int)(100f * i / numFiles)
			if (progress > oldProgress) {
				pluginHost.updateProgress("Progress", progress)
				oldProgress = progress
			}
			// check to see if the user has requested a cancellation
			if (pluginHost.isRequestForOperationCancelSet()) {
				pluginHost.showFeedback("Operation cancelled")
				return
			}
	    }
	    
		String inputFileExtension = FileUtilities.getFileExtension(inputFiles[0])
        String outputHeader = inputFiles[0].replace(".${inputFileExtension}", suffix + ".dep");
        pluginHost.returnData(outputHeader);

		long end = System.currentTimeMillis()  
		double duration = (end - start) / 1000.0
		pluginHost.showFeedback("Interpolation completed in " + duration + " seconds.\n" + 
		  numSuccessfulInterpolations + " tiles were successfully interpolated.\n" + 
		  "One has been displayed on the map.")
		
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

    public List<Future<Boolean>> getExecutorResults(ExecutorService executor, ArrayList<DoWork> tasks) {
    	List<Future<Boolean>> results = executor.invokeAll(tasks);
		return results
    }

	class DoWork implements Callable<Boolean> {
		private int tileNum
		private boolean[] classValuesToExclude
		private BoundingBox[] bb
		private String[] inputFiles
		private String suffix
	    private double resolution
	    private double maxSlope
	    private double maxScanAngleDeviation
	    private double weight
		private double maxDist
		private double minHeight
		
	    DoWork(int tileNum, String[] inputFiles, String suffix, 
	      BoundingBox[] bb, double resolution, double maxSlope,
	      double maxScanAngleDeviation, double weight, 
	      double maxDist, double minHeight) {
			this.tileNum = tileNum;
			this.inputFiles = inputFiles.clone();
			this.suffix = suffix;
			this.bb = bb.clone();
			this.resolution = resolution;
			this.maxSlope = maxSlope;
			this.maxScanAngleDeviation = maxScanAngleDeviation
			this.weight = weight
			this.maxDist = maxDist
			this.minHeight = minHeight
       	}
        	
        @Override
        @CompileStatic
	    public Boolean call() {
	    	// check to see if the user has requested a cancellation
			if (pluginHost.isRequestForOperationCancelSet()) {
				return Boolean.FALSE
			}

			String inputFile = inputFiles[tileNum]
	    	int numFiles = inputFiles.length
	    	double x, y, z;
			double easting, northing;
			double scanAngle;
			final double noData = -32768;
			int i;
			List<KdTree.Entry<InterpolationRecord>> results;
		
			LASReader las;
			BoundingBox expandedBB = new BoundingBox(bb[tileNum].getMinX() - resolution, bb[tileNum].getMinY() - resolution, bb[tileNum].getMaxX() + resolution, bb[tileNum].getMaxY() + resolution);
			
			// count how many valid points there are
			int numPoints = 0;
			byte classVal
			ArrayList<PointRecord> recs = new ArrayList<>();
			for (int a = 0; a < numFiles; a++) {
				if (bb[a].entirelyContainedWithin(expandedBB) || 
				  bb[a].intersectsAnEdgeOf(expandedBB)) {
				 	las = new LASReader(inputFiles[a])
				 	ArrayList<PointRecord> points = las.getPointRecordsInBoundingBox(expandedBB)
			 		for (PointRecord point : points) {
				 		classVal = point.getClassification();
				 		if (!point.isPointWithheld()
				 		    && point.getReturnNumber() == point.getNumberOfReturns()
				 		    && classVal != 4 && classVal != 5 
				 		    && classVal != 6 && classVal != 7 && point.z >= minHeight) {
                            recs.add(point);
                        }
			 		}
				 	points.clear();
				}
			}
			
			numPoints = recs.size();
			KdTree<InterpolationRecord> pointsTree = new KdTree.SqrEuclid<InterpolationRecord>(2, new Integer(numPoints))

			BooleanBitArray1D nongroundBitArray = new BooleanBitArray1D(numPoints)
			
			double[] entry;
			PointRecColours pointColours;
			i = 0
			for (PointRecord point : recs) {
				entry = [point.getY(), point.getX()]
				pointsTree.addPoint(entry, new InterpolationRecord(point.getX(), point.getY(), point.getZ(), point.getScanAngle(), i));
				i++
			}
            
			recs.clear();
			
            // create the output grid
            String inputFileExtension = FileUtilities.getFileExtension(inputFile)
            String outputHeader = inputFile.replace(".${inputFileExtension}", suffix + ".dep");
            
	        // see if the output files already exist, and if so, delete them.
	        if ((new File(outputHeader)).exists()) {
	            (new File(outputHeader)).delete();
	            (new File(outputHeader.replace(".dep", ".tas"))).delete();
	        }
	
	        // What are north, south, east, and west and how many rows and 
	        // columns should there be?
	        double west = bb[tileNum].getMinX() - 0.5 * resolution;
	        double north = bb[tileNum].getMaxY() + 0.5 * resolution;
	        int nrows = (int) (Math.ceil((north - bb[tileNum].getMinY()) / resolution));
	        int ncols = (int) (Math.ceil((bb[tileNum].getMaxX() - west) / resolution));
	        double south = north - nrows * resolution;
	        double east = west + ncols * resolution;
	
	        try {
	            // Create the whitebox raster object.
                WhiteboxRaster image = new WhiteboxRaster(outputHeader, 
                	north, south, east, west, nrows, ncols, 
                	WhiteboxRasterBase.DataScale.CONTINUOUS,
                    WhiteboxRasterBase.DataType.FLOAT, noData, noData);
                image.setPreferredPalette("spectrum.plt")

                final double radToDeg = 180.0 / Math.PI
                double slopeThreshold = maxSlope / radToDeg
	            InterpolationRecord value;
                double dist, val, minVal, maxVal;
	            double minDist, minDistVal;
	            double halfResolution = resolution / 2;
	            int oldProgress = -1;
	            int progress;
                for (int row in 0..(nrows - 1)) {
                    for (int col in 0..(ncols - 1)) {
                        easting = image.getXCoordinateFromColumn(col);
                        northing = image.getYCoordinateFromRow(row);
                        entry = [northing, easting];
                        results = pointsTree.neighborsWithinRange(entry, maxDist);
                        if (results.size() > 1) {

							if (maxScanAngleDeviation < 90.0) {
								double minScanAngle = Double.POSITIVE_INFINITY;
		                        double maxScanAngle = Double.NEGATIVE_INFINITY;
		                        for (i = 0; i < results.size(); i++) {
		                        	value = (InterpolationRecord)results.get(i).value;
		                            scanAngle = value.scanAngle;
		                            if (scanAngle > maxScanAngle) { maxScanAngle = scanAngle; }
		                            if (scanAngle < minScanAngle) { minScanAngle = scanAngle; }
		                        }
		                        for (i = 0; i < results.size(); i++) {
		                        	value = (InterpolationRecord)results.get(i).value;
		                            scanAngle = value.scanAngle;
		                            if ((scanAngle - minScanAngle) > maxScanAngleDeviation) {
		                                nongroundBitArray.setValue(value.getIndex(), true)
		                            }
		                        }
							}
							
							// eliminate all non-ground points based on slope
							int n = results.size()
	                        double slope
	                        //int higherPoint
	                        int higherPointIndex
	                        double higherVal, lowerVal
	                        InterpolationRecord rec1, rec2
							for (i = 0; i < n - 1; i++) {
								rec1 = (InterpolationRecord)(results.get(i).value)
								if (!nongroundBitArray.getValue(rec1.getIndex())) {
									for (int j = i + 1; j < n; j++) {
										rec2 = (InterpolationRecord)(results.get(j).value)
										if (!nongroundBitArray.getValue(rec2.getIndex())) {
											dist = Math.sqrt((rec1.getX() - rec2.getX()) * (rec1.getX() - rec2.getX()) + (rec1.getY() - rec2.getY()) * (rec1.getY() - rec2.getY()))
											
											if (rec1.getValue() > rec2.getValue()) {
												higherVal = rec1.getValue()
												lowerVal = rec2.getValue()
												//higherPoint = i
												higherPointIndex = rec1.getIndex()
											} else {
												higherVal = rec2.getValue()
												lowerVal = rec1.getValue()
												//higherPoint = j
												higherPointIndex = rec2.getIndex()
											}
											slope = Math.atan((higherVal - lowerVal) / dist)
											if (slope > slopeThreshold) {
												nongroundBitArray.setValue(higherPointIndex, true)
											}
										}
									}
								}
							}

//							// now find the nearest ground point and assign it as the z
//	                        minDist = Double.POSITIVE_INFINITY
//	                        z = noData
//	                        n = 0
//	                        for (i = 0; i < results.size(); i++) {
//	                        	rec1 = (InterpolationRecord)results.get(i).value
//	                        	if (!nongroundBitArray.getValue(rec1.getIndex())) {
//	                        		dist = results.get(i).distance
//		                        	val = rec1.value
//		                            if (dist < minDist) { 
//		                            	minDist = dist
//		                            	z = val
//		                            }
//		                            n++
//	                        	}
//	                        }

							// now perform an IDW interpolation
	                        z = noData
	                        n = 0
	                        double sumWeights = 0
//	                        double weight = 2.0
	                        ArrayList<Double> weights = new ArrayList<>()
	                        ArrayList<Double> vals = new ArrayList<>()
	                        for (i = 0; i < results.size(); i++) {
	                        	rec1 = (InterpolationRecord)results.get(i).value
	                        	if (!nongroundBitArray.getValue(rec1.getIndex())) {
	                        		if (results.get(i).distance > 0) {
		                        		dist = 1 / Math.pow(Math.sqrt(results.get(i).distance), weight)
			                        	weights.add(dist)
			                            sumWeights += dist
			                            vals.add(rec1.value)
			                            n++
	                        		} else {
	                        			weights = new ArrayList<>()
	                        			vals = new ArrayList<>()
	                        			weights.add(1.0d)
			                            sumWeights = 1.0
			                            vals.add(rec1.value)
			                            n = 1
	                        			break
	                        		}
	                        	}
	                        }
	                        if (n > 0) {
	                        	z = 0
	                        	for (int s = 0; s < n; s++) {
	                        		z += (weights.get(s) * vals.get(s)) / sumWeights
	                        	}
	                        }
	                        
	                        image.setValue(row, col, z);
                        } else if (results.size() == 1) {
                        	value = (InterpolationRecord)results.get(0).value
                        	image.setValue(row, col, value.value);
                        } else {
                        	image.setValue(row, col, noData);
                        }
                    }

                    if (numFiles == 1) {
                    	progress = (int)(100f * row / (nrows - 1))
                    	if (progress > oldProgress) {
                    		oldProgress = progress
                    		pluginHost.updateProgress("Interpolation Progress:", progress)
                    	}
                    }
                    
                    // check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return Boolean.FALSE
					}
                }

                image.addMetadataEntry("Created by the " + descriptiveName + " tool.");
                image.addMetadataEntry("Created on " + new Date());
                image.addMetadataEntry("IDW weight: $weight");
                image.addMetadataEntry("Slope threshold: $maxSlope");
                image.addMetadataEntry("Max search dist: $maxDist");
                image.close();

	        } catch (Exception e) {
	            pluginHost.showFeedback(e.getMessage());
	            return Boolean.FALSE;
	        }

			int solved = numSolvedTiles.incrementAndGet() //++;
			int progress = (int) (100f * solved / numFiles)
			pluginHost.updateProgress("Interpolated $solved tiles:", progress)
	        return Boolean.TRUE;
        }                
    }

	@CompileStatic
	class InterpolationRecord {
        
        double value;
        byte scanAngle;
        double x;
        double y;
        int index;
        
        InterpolationRecord(double x, double y, double value, byte scanAngle, int index) {
            this.value = value;
            this.scanAngle = (byte)Math.abs(scanAngle);
            this.x = x;
            this.y = y;
            this.index = index;
        }

        double getValue() {
            return value;
        }
        
        byte getScanAngle() {
            return scanAngle;
        }

        double getX() {
        	return x;
        }

        double getY() {
        	return y;
        }     

        int getIndex() {
        	return index;
        }
    }
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def myClass = new LiDAR_BareEarthDEM(pluginHost, args, name, descriptiveName)
}
