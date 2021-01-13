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
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
//def name = "LiDAR_BareEarthDEM"
def descriptiveName = "Bare-Earth DEM (LiDAR)"
def description = "Interpolates LAS files into a bare-Earth digital elevation model (DEM)."
def toolboxes = ["LidarTools"]

public class LiDAR_BareEarthDEM implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	private AtomicInteger numSolvedRows = new AtomicInteger(0)
	
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
			sd.addDialogFile("Output file", "Output Raster File:", "saveAs", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogDataInput("IDW Exponent", "IDW Exponent", "2", true, false)
			sd.addDialogDataInput("Max Search Distance (m)", "Max Search Distance (m)", "", true, false)
			sd.addDialogDataInput("Grid Resolution (m)", "Grid Resolution:", "", true, false)
			sd.addDialogDataInput("Threshold in the slope between points to define an off-terrain point.", "Inter-point Slope Threshold:", "30.0", true, false)
			sd.addDialogDataInput("Max Scan Angle Deviation (optional)", "Max Scan Angle Deviation (optional)", "", true, true)
			
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
	  	if (args.length != 7) {
			pluginHost.showFeedback("Incorrect number of arguments given to tool.")
			return
		}
		
		// read the input parameters
		final String inputFileString = args[0]
		String outputFile = args[1]
        double weight = Double.parseDouble(args[2])
		double maxDist = Double.parseDouble(args[3])
		double resolution = Double.parseDouble(args[4])
        double maxSlope = Double.parseDouble(args[5])
        double maxScanAngleDeviation = 1000.0d
        if (args[6] != null && !args[6].isEmpty()) {
        	if (!args[6].toLowerCase().equals("not specified")) {
        		maxScanAngleDeviation= Double.parseDouble(args[6])
        	}
        }
        if (maxScanAngleDeviation < 1.0) { maxScanAngleDeviation = 1.0; }
        
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
		
		final int numFiles = inputFiles.length
		BoundingBox[] bb = new BoundingBox[numFiles]
		LASReader[] lasReaders = new LASReader[numFiles]
		i = 0
		double minX = Double.POSITIVE_INFINITY
		double minY = Double.POSITIVE_INFINITY
		double maxX = Double.NEGATIVE_INFINITY
		double maxY = Double.NEGATIVE_INFINITY
		for (String inputFile : inputFiles) {
			LASReader las = new LASReader(inputFile)
			lasReaders[i] = las
			bb[i] = new BoundingBox(las.getMinX(), las.getMinY(), las.getMaxX(), las.getMaxY());
			if (las.getMinX() < minX) { minX = las.getMinX(); }
			if (las.getMaxX() > maxX) { maxX = las.getMaxX(); }
			if (las.getMinY() < minY) { minY = las.getMinY(); }
			if (las.getMaxY() > maxY) { maxY = las.getMaxY(); }
			i++
		}


		// What are north, south, east, and west and how many rows and 
        // columns should there be?
        double west = minX - 0.5 * resolution;
        double north = maxY + 0.5 * resolution;
        int nrows = (int) (Math.ceil((north - minY) / resolution));
        int ncols = (int) (Math.ceil((maxX - west) / resolution));
        double south = north - nrows * resolution;
        double east = west + ncols * resolution;
        final double noData = -32768;

        ArrayList<DoWork> tasks = new ArrayList<>();

        // Create the whitebox raster object.
        WhiteboxRaster image = new WhiteboxRaster(outputFile, 
        	north, south, east, west, nrows, ncols, 
        	WhiteboxRasterBase.DataScale.CONTINUOUS,
            WhiteboxRasterBase.DataType.FLOAT, noData, noData);
        image.setPreferredPalette("spectrum.plt")
        image.setForceAllDataInMemory(true)

		ArrayList<PointRecord>[] dataArray = new ArrayList<PointRecord>[nrows];
		int offset = 0 //(int)(maxDist / resolution)
		//if (offset < 1) { offset = 1 }
		int progress
		int oldProgress = -1
		for (int a = 0; a < numFiles; a++) {
			long numPoints = lasReaders[a].getNumPointRecords()
			int n = 0
			for (int p = 0; p < numPoints; p++) {
				PointRecord point = lasReaders[a].getPointRecord(p)
				int row = image.getRowFromYCoordinate(point.getY())
				for (int k = (row - offset); k <= (row + offset); k++) {
					if (k >= 0 && k < nrows) {
						if (dataArray[k] == null) { dataArray[k] = new ArrayList<PointRecord>() }
						dataArray[k].add(point)
					}
				}
				n++
				progress = (int)(100f * n / numPoints)
				if (progress > oldProgress) {
					oldProgress = progress
					pluginHost.updateProgress("Getting data (${a + 1} of ${numFiles}):", progress)
	
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}
		}

		progress = 0
        oldProgress = -1
        for (int row in 0..(nrows - 1)) {
        	if (dataArray[row] != null && dataArray[row].size() > 0) {
        		tasks.add(new DoWork(row, dataArray[row], resolution, maxSlope, maxScanAngleDeviation, weight, maxDist, image))
        	}
        	progress = (int)(100f * row / (nrows - 1))
			if (progress > oldProgress) {
				oldProgress = progress
				pluginHost.updateProgress("Getting data:", progress)

				// check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
			}
        }

//		east += maxDist; // = image.getXCoordinateFromColumn(ncols - 1) + maxDist;
//        west -= maxDist; // = image.getXCoordinateFromColumn(0) - maxDist;
//        int progress = 0
//        int oldProgress = -1
//        long n = 0        
//        for (int row in 0..(nrows - 1)) {
//            // get all of the points that are within this row
//            north = image.getYCoordinateFromRow(row) + maxDist;
//            south = north - 2 * maxDist; //image.getYCoordinateFromRow(row) - maxDist;
//            BoundingBox rowBB = new BoundingBox(west, south, east, north);
//			
//			// count how many valid points there are
//			int numPoints = 0;
//			ArrayList<PointRecord> recs = new ArrayList<>();
//			for (int a = 0; a < numFiles; a++) {
				//if (bb[a].entirelyContainedWithin(rowBB) || 
				  //bb[a].intersectsAnEdgeOf(rowBB)) {
				 	//LASReader las = new LASReader(inputFiles[a])
//				 	ArrayList<PointRecord> points = lasReaders[a].getPointRecordsInBoundingBox(rowBB)
				 	//n++
				 	//pluginHost.showFeedback("${north} ${south} ${east} ${west}")
//				 	for (PointRecord point : points) {
//				 		if (!point.isPointWithheld()) {
//                            recs.add(point);
//                        }
//			 		}
//			 		n += points.size()
//				 	points.clear();
				//}
//			}
//			if (recs.size() > 0) {
//				tasks.add(new DoWork(row, recs, resolution, maxSlope, maxScanAngleDeviation, weight, maxDist, image))
//				//n++
//			}
//			progress = (int)(100f * row / (nrows - 1))
//			if (progress > oldProgress) {
//				oldProgress = progress
//				pluginHost.updateProgress("Getting data:", progress)
//
//				// check to see if the user has requested a cancellation
//				if (pluginHost.isRequestForOperationCancelSet()) {
//					pluginHost.showFeedback("Operation cancelled")
//					return
//				}
//			}
//        }
		//pluginHost.showFeedback("${n} points")


		
		pluginHost.updateProgress("Please wait...", 0)
		
        ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
  	    // the only reason for the getExecutorResults method 
  	    // is that Groovy throws a compilation type mis-match
  	    // error when compiled statically. I think it's a bug.
  	    List<Future<Boolean>> results = getExecutorResults(executor, tasks);
        executor.shutdown();

//        i = 0
        for (Future<RowNumberAndData> result : results) {
    		RowNumberAndData data = result.get()
    		if (data != null) {
    			int row = data.getRow()
    			//image.setRowValues(row, data.getData())
    			double[] outputData = data.getData()
    			for (int col = 0; col < ncols; col++) {
    				image.setValue(row, col, outputData[col])
    			}
    			dataArray[row].clear()
    		}
//        	i++
//			// update progress bar
//			progress = (int)(100f * i / numFiles)
//			if (progress > oldProgress) {
//				pluginHost.updateProgress("Progress", progress)
//				oldProgress = progress
//			}
//			// check to see if the user has requested a cancellation
//			if (pluginHost.isRequestForOperationCancelSet()) {
//				pluginHost.showFeedback("Operation cancelled")
//				return
//			}
	    }

	    image.addMetadataEntry("Created by the " + descriptiveName + " tool.");
        image.addMetadataEntry("Created on " + new Date());
        image.addMetadataEntry("Max. Search Distance:\t${maxDist}");
        image.addMetadataEntry("Resolution:\t${resolution}");
        image.addMetadataEntry("Max. Slope:\t${maxSlope}");
        image.close();
	    
		pluginHost.returnData(outputFile);

		long end = System.currentTimeMillis()  
		double duration = (end - start) / 1000.0
		pluginHost.showFeedback("Interpolation completed in " + duration + " seconds.")
		
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

	class DoWork implements Callable<RowNumberAndData> {
		private int row
		private double resolution
	    private double maxSlope
	    private double maxScanAngleDeviation
	    private double weight
		private double maxDist
		private ArrayList<PointRecord> recs
		private WhiteboxRaster image
		
	    DoWork(int row, ArrayList<PointRecord> recs, double resolution, double maxSlope,
	      double maxScanAngleDeviation, double weight, double maxDist, WhiteboxRaster image) {
	      	this.row = row;
	      	this.recs = recs;
			this.resolution = resolution;
			this.maxSlope = maxSlope;
			this.maxScanAngleDeviation = maxScanAngleDeviation
			this.weight = weight
			this.maxDist = maxDist
			this.image = image
       	}
        	
        @Override
        @CompileStatic
	    public RowNumberAndData call() {
	    	// check to see if the user has requested a cancellation
			if (pluginHost.isRequestForOperationCancelSet()) {
				return null
			}
			final double noData = -32768;
			double z
			int nrows = image.getNumberRows()
			int ncols = image.getNumberColumns()
			int numPoints = recs.size();
			KdTree<InterpolationRecord> pointsTree = new KdTree.SqrEuclid<InterpolationRecord>(2, new Integer(numPoints))

			BooleanBitArray1D nongroundBitArray = new BooleanBitArray1D(numPoints)
			
			double[] entry;
			List<KdTree.Entry<InterpolationRecord>> results;
			int i = 0
			for (PointRecord point : recs) {
				entry = [point.getY(), point.getX()]
				pointsTree.addPoint(entry, new InterpolationRecord(point.getX(), point.getY(), point.getZ(), point.getScanAngle(), i));
				i++
			}
            
			recs.clear();
			
            double[] retData = new double[ncols]
                
	        try {
	            final double radToDeg = 180.0 / Math.PI
                double slopeThreshold = maxSlope / radToDeg
	            InterpolationRecord value;
                double dist, val, minVal, maxVal;
	            double minDist, minDistVal;
//	            double halfResolution = resolution / 2;
	            
                for (int col in 0..(ncols - 1)) {
                    double easting = image.getXCoordinateFromColumn(col);
                    double northing = image.getYCoordinateFromRow(row);
                    entry = [northing, easting];
                    results = pointsTree.neighborsWithinRange(entry, maxDist);
                    if (results.size() > 1) {

						if (maxScanAngleDeviation < 90.0) {
							double minScanAngle = Double.POSITIVE_INFINITY;
	                        double maxScanAngle = Double.NEGATIVE_INFINITY;
	                        for (i = 0; i < results.size(); i++) {
	                        	value = (InterpolationRecord)results.get(i).value;
	                            double scanAngle = value.scanAngle;
	                            if (scanAngle > maxScanAngle) { maxScanAngle = scanAngle; }
	                            if (scanAngle < minScanAngle) { minScanAngle = scanAngle; }
	                        }
	                        for (i = 0; i < results.size(); i++) {
	                        	value = (InterpolationRecord)results.get(i).value;
	                            double scanAngle = value.scanAngle;
	                            if ((scanAngle - minScanAngle) > maxScanAngleDeviation) {
	                                nongroundBitArray.setValue(value.getIndex(), true)
	                            }
	                        }
						}
						
						// eliminate all non-ground points based on slope
						int n = results.size()
                        double slope
                        int higherPoint, higherPointIndex
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
											higherPoint = i
											higherPointIndex = rec1.getIndex()
										} else {
											higherVal = rec2.getValue()
											lowerVal = rec1.getValue()
											higherPoint = j
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

						// now find the nearest ground point and assign it as the z
                        minDist = Double.POSITIVE_INFINITY
                        z = noData
                        n = 0
                        for (i = 0; i < results.size(); i++) {
                        	rec1 = (InterpolationRecord)results.get(i).value
                        	if (!nongroundBitArray.getValue(rec1.getIndex())) {
                        		dist = results.get(i).distance
	                        	val = rec1.value
	                            if (dist < minDist) { 
	                            	minDist = dist
	                            	z = val
	                            }
	                            n++
                        	}
                        }
	                        
//						// now perform an IDW interpolation
//                        z = noData
//                        n = 0
//                        double sumWeights = 0
//                        ArrayList<Double> weights = new ArrayList<>()
//                        ArrayList<Double> vals = new ArrayList<>()
//                        for (i = 0; i < results.size(); i++) {
//                        	rec1 = (InterpolationRecord)results.get(i).value
//                        	if (!nongroundBitArray.getValue(rec1.getIndex())) {
//                        		if (results.get(i).distance > 0) {
//	                        		dist = 1 / Math.pow(Math.sqrt(results.get(i).distance), weight)
//		                        	weights.add(dist)
//		                            sumWeights += dist
//		                            vals.add(rec1.value)
//		                            n++
//                        		} else {
//                        			weights = new ArrayList<>()
//                        			vals = new ArrayList<>()
//                        			weights.add(1.0d)
//		                            sumWeights += 1.0
//		                            vals.add(rec1.value)
//		                            n = 1
//                        			break
//                        		}
//                        	}
//                        }
//                        if (n > 0) {
//                        	z = 0
//                        	for (int s = 0; s < n; s++) {
//                        		z += (weights.get(s) * vals.get(s)) / sumWeights
//                        	}
//                        }

                        retData[col] = z;
                    } else if (results.size() == 1) {
                    	value = (InterpolationRecord)results.get(0).value
                    	retData[col] = value.value;
                    } else {
                    	retData[col] = noData;
                    }
                }

	        } catch (Exception e) {
	            pluginHost.showFeedback(e.getMessage());
	            return null;
	        }

			int solved = numSolvedRows.incrementAndGet()
			int progress = (int) (100f * solved / (nrows - 1))
			pluginHost.updateProgress("Progress:", progress)
			
			def ret = new RowNumberAndData(row, retData)
	        return ret;
        }                
    }

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

    @CompileStatic
    class RowNumberAndData {
		private int row
		private double[] data
		private int numDataEntries = 0
		
    	RowNumberAndData(int row, double[] data) {
    		this.row = row
    		this.data = data
    		this.numDataEntries = data.length
    	}

    	public int getRow() {
    		return row
    	}

    	public void setRow(int row) {
    		this.row = row
    	}

    	public double[] getData() {
    		return data
    	}

		public void setData(double[] data) {
			this.data = data
		}
		
    	public double getDataAt(int n) {
    		if (n < numDataEntries) {
    			return data[n]
    		} else {
    			return null
    		}
    	}
    }
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def myClass = new LiDAR_BareEarthDEM(pluginHost, args, name, descriptiveName)
}
