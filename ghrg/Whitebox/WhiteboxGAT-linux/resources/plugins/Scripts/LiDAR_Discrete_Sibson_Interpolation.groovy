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
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.structures.BoundingBox
import whitebox.structures.KdTree
import whitebox.structures.SimpleGridCell
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import whitebox.utilities.FileUtilities;
import whitebox.structures.DoubleArray2D
import whitebox.structures.IntArray2D
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "LiDAR_Discrete_Sibson_Interpolation"
def descriptiveName = "Discrete Sibson Interpolation (LiDAR)"
def description = "Interpolates LAS files using an discrete Sibson (Natural Neighbour) scheme."
def toolboxes = ["LidarTools"]

public class LiDAR_Discrete_Sibson_Interpolation implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	private AtomicInteger numSolvedTiles = new AtomicInteger(0)
	
	public LiDAR_Discrete_Sibson_Interpolation(WhiteboxPluginHost pluginHost, 
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
			sd.addDialogDataInput("Output File Suffix (e.g. _lastReturn)", "Output File Suffix (e.g. _lastReturn)", "", false, false)
			sd.addDialogComboBox("Interpolation Paramter:", "Interpolation Parameter:", ["Z (elevation)", "Intensity", "Scan Angle"], 0)
			sd.addDialogComboBox("Point Return:", "Point Return:", ["All Points", "First Return", "Last Return"], 0)
			sd.addDialogDataInput("Grid Resolution (m)", "Grid Resolution (m)", "", true, false)
			sd.addDialogDataInput("Discretization factor (5-20)", "Discretization factor (5-20)", "10", true, false)
			sd.addDialogLabel("Exclude points with the following classification")
			sd.addDialogLabel("values from the interpolation:")
			sd.addDialogCheckBox("Never Classified", "Never Classified", false)
			sd.addDialogCheckBox("Unclassified", "Unclassified", false)
			sd.addDialogCheckBox("Bare Ground", "Bare Ground", false)
			sd.addDialogCheckBox("Low Vegetation", "Low Vegetation", false)
			sd.addDialogCheckBox("Medium Vegetation", "Medium Vegetation", false)
			sd.addDialogCheckBox("High Vegetation", "High Vegetation", false)
			sd.addDialogCheckBox("Building", "Building", false)
			sd.addDialogCheckBox("Low Points", "Low Points", false)
			sd.addDialogCheckBox("Model Key Points", "Model Key Points", false)
			sd.addDialogCheckBox("Water", "Water", false)
			sd.addDialogDataInput("Minimum height (Optional: used to exclude points).", "Minimum Height (optional):", "", true, true)
			sd.addDialogDataInput("Maximum height (Optional: used to exclude points).", "Maximum Height (optional):", "", true, true)
			
			
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
	  	if (args.length < 16) {
			pluginHost.showFeedback("Incorrect number of arguments given to tool.")
			return
		}
		
		// read the input parameters
		final String inputFileString = args[0]
		String suffix = ""
		if (args[1] != null) {
			suffix = args[1].trim();
		}
        String whatToInterpolate = args[2].toLowerCase();
        String returnNumberToInterpolate = args[3].toLowerCase();
        double resolution = Double.parseDouble(args[4]);
        int discretizationFactor = (int)(Double.parseDouble(args[5]));
		final boolean excludeNeverClassified = Boolean.parseBoolean(args[6]);
        final boolean excludeUnclassified = Boolean.parseBoolean(args[7]);
        final boolean excludeBareGround = Boolean.parseBoolean(args[8]);
        final boolean excludeLowVegetation = Boolean.parseBoolean(args[9]);
        final boolean excludeMediumVegetation = Boolean.parseBoolean(args[10]);
        final boolean excludeHighVegetation = Boolean.parseBoolean(args[11]);
        final boolean excludeBuilding = Boolean.parseBoolean(args[12]);
        final boolean excludeLowPoint = Boolean.parseBoolean(args[13]);
        final boolean excludeModelKeyPoint = Boolean.parseBoolean(args[14]);
        final boolean excludeWater = Boolean.parseBoolean(args[15]);

        double minHeight = Double.NEGATIVE_INFINITY 
        double maxHeight = Double.POSITIVE_INFINITY
        if (args.length == 18) { // min and max heights have been specified
	        if (!args[16].toLowerCase().contains("not")) {
	        	minHeight = Double.parseDouble(args[16]);
	        }
	        if (!args[17].toLowerCase().contains("not")) {
	        	maxHeight = Double.parseDouble(args[17]);
	        }
        }

        boolean[] classValuesToExclude = new boolean[32]; // there can be up to 32 different classes in future versions

        if (excludeNeverClassified) {
            classValuesToExclude[0] = true;
        }
        if (excludeUnclassified) {
            classValuesToExclude[1] = true;
        }
        if (excludeBareGround) {
            classValuesToExclude[2] = true;
        }
        if (excludeLowVegetation) {
            classValuesToExclude[3] = true;
        }
        if (excludeMediumVegetation) {
            classValuesToExclude[4] = true;
        }
        if (excludeHighVegetation) {
            classValuesToExclude[5] = true;
        }
        if (excludeBuilding) {
            classValuesToExclude[6] = true;
        }
        if (excludeLowPoint) {
            classValuesToExclude[7] = true;
        }
        if (excludeModelKeyPoint) {
            classValuesToExclude[8] = true;
        }
        if (excludeWater) {
            classValuesToExclude[9] = true;
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
	      		 bb, whatToInterpolate, returnNumberToInterpolate, 
	      		 resolution, discretizationFactor, classValuesToExclude, 
	      		 minHeight, maxHeight))
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
	    private String whatToInterpolate
	    private String returnNumberToInterpolate
	    private int discretizationFactor
		private double resolution
		private double minVal = Double.NEGATIVE_INFINITY;
		private double maxVal = Double.POSITIVE_INFINITY;
		
	    DoWork(int tileNum, String[] inputFiles, String suffix, 
	      BoundingBox[] bb, String whatToInterpolate, 
	      String returnNumberToInterpolate, double resolution, int discretizationFactor, 
	      boolean[] classValuesToExclude, double minHeight, double maxHeight) {
			this.tileNum = tileNum;
			this.inputFiles = inputFiles.clone();
			this.suffix = suffix;
			this.bb = bb.clone();
			this.whatToInterpolate = whatToInterpolate;
			this.returnNumberToInterpolate = returnNumberToInterpolate;
			this.discretizationFactor = discretizationFactor;
			this.resolution = resolution;
			this.classValuesToExclude = classValuesToExclude.clone();
			this.minVal = minHeight;
			this.maxVal = maxHeight;
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
			int row, col, i;
			int oldProgress = -1;
	        int progress;    
			KdTree.Entry<Double> result;
			double value;
			double distance;
			double maxDist2 = resolution * 2;
                
			LASReader las;
			BoundingBox expandedBB = new BoundingBox(bb[tileNum].getMinX() - maxDist2, bb[tileNum].getMinY() - maxDist2, bb[tileNum].getMaxX() + maxDist2, bb[tileNum].getMaxY() + maxDist2);
			
			// count how many valid points there are
			int numPoints = 0;
			ArrayList<PointRecord> recs = new ArrayList<>();
			for (int a = 0; a < numFiles; a++) {
				if (bb[a].entirelyContainedWithin(expandedBB) || 
				  bb[a].intersectsAnEdgeOf(expandedBB)) {
				 	las = new LASReader(inputFiles[a])
				 	ArrayList<PointRecord> points = las.getPointRecordsInBoundingBox(expandedBB)
				 	int numLasPoints = points.size();
				 	int lasPoint = 0
				 	if (returnNumberToInterpolate.equals("all points")) {
				 		for (PointRecord point : points) {
				 			z = point.getZ()
					 		if (!point.isPointWithheld() &&
	                          !(classValuesToExclude[point.getClassification()]) &&
	                          z > minVal && z < maxVal) {
	                            recs.add(point);
	                        }

	                        if (numFiles == 1) {
	                        	lasPoint++
			                	progress = (int)(100f * lasPoint / numLasPoints)
			                	if (progress > oldProgress) {
			                		oldProgress = progress
			                		pluginHost.updateProgress("Reading Data:", progress)

					                // check to see if the user has requested a cancellation
									if (pluginHost.isRequestForOperationCancelSet()) {
										pluginHost.showFeedback("Operation cancelled")
										return Boolean.FALSE
									}
			                	}
			                }
				 		}
				 	} else if (returnNumberToInterpolate.equals("first return")) {
				 		for (PointRecord point : points) {
					 		z = point.getZ()
					 		if (!point.isPointWithheld() &&
	                          !(classValuesToExclude[point.getClassification()])
                              && point.getReturnNumber() == 1 && z > minVal && z < maxVal) {
	                            recs.add(point);
	                        }
	                        
	                        if (numFiles == 1) {
	                        	lasPoint++
			                	progress = (int)(100f * lasPoint / numLasPoints)
			                	if (progress > oldProgress) {
			                		oldProgress = progress
			                		pluginHost.updateProgress("Reading Data:", progress)

					                // check to see if the user has requested a cancellation
									if (pluginHost.isRequestForOperationCancelSet()) {
										pluginHost.showFeedback("Operation cancelled")
										return Boolean.FALSE
									}
			                	}
			                }
				 		}
				 	} else { // if (returnNumberToInterpolate.equals("last return")) {
				 		for (PointRecord point : points) {
					 		z = point.getZ()
					 		if (!point.isPointWithheld() &&
	                          !(classValuesToExclude[point.getClassification()])
                              && point.getReturnNumber() == point.getNumberOfReturns()
                              && z > minVal && z < maxVal) {
	                            recs.add(point);
	                        }

	                        if (numFiles == 1) {
	                        	lasPoint++
			                	progress = (int)(100f * lasPoint / numLasPoints)
			                	if (progress > oldProgress) {
			                		oldProgress = progress
			                		pluginHost.updateProgress("Reading Data:", progress)

					                // check to see if the user has requested a cancellation
									if (pluginHost.isRequestForOperationCancelSet()) {
										pluginHost.showFeedback("Operation cancelled")
										return Boolean.FALSE
									}
			                	}
			                }
				 		}
				 	}
				 	points.clear();
				}
			}
			
			numPoints = recs.size();
			KdTree<Double> pointsTree = new KdTree.SqrEuclid<Double>(2, new Integer(numPoints))

			double[] entry;
			PointRecColours pointColours;
			if (whatToInterpolate.equals("z (elevation)")) {
				for (PointRecord point : recs) {
					entry = [point.getY(), point.getX()]
					pointsTree.addPoint(entry, new Double(point.getZ()));
				}
            } else if (whatToInterpolate.equals("intensity")) {
                for (PointRecord point : recs) {
					entry = [point.getY(), point.getX()]
					pointsTree.addPoint(entry, new Double(point.getIntensity()));
				}
            } else if (whatToInterpolate.equals("classification")) {
                for (PointRecord point : recs) {
					entry = [point.getY(), point.getX()]
					pointsTree.addPoint(entry, new Double(point.getClassification()));
				}
            } else if (whatToInterpolate.equals("scan angle")) {
                for (PointRecord point : recs) {
					entry = [point.getY(), point.getX()]
					pointsTree.addPoint(entry, new Double(point.getScanAngle()));
				}
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
	        int rowsLessOne = nrows - 1;
	        int ncols = (int) (Math.ceil((bb[tileNum].getMaxX() - west) / resolution));
	        double south = north - nrows * resolution;
	        double east = west + ncols * resolution;
	        double halfResolution = resolution / 2.0;
			IntArray2D nArray = new IntArray2D(nrows, ncols, 0);
			DoubleArray2D sumArray = new DoubleArray2D(nrows, ncols, noData);

			int numCells = nrows * ncols;
			List<KdTree.Entry<SimpleGridCell>> results;
			KdTree<SimpleGridCell> gridTree = new KdTree.SqrEuclid<SimpleGridCell>(2, new Integer(numCells))
			oldProgress = -1
			for (row = 0; row < nrows; row++) {
                northing = (north - halfResolution) - (row * resolution);
                for (col = 0; col < ncols; col++) {
					easting = (col * resolution) + (west + halfResolution);
                    entry = [northing, easting]
					gridTree.addPoint(entry, new SimpleGridCell(row, col));    
                }

                if (numFiles == 1) {
                	progress = (int)(100f * row / rowsLessOne)
                	if (progress > oldProgress) {
                		oldProgress = progress
                		pluginHost.updateProgress("Initializing:", progress)
                	}
                }
                
                // check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return Boolean.FALSE
				}
			}

	        int fineRows = nrows * discretizationFactor;
	        int fineRowsLessOne = fineRows - 1;
			int fineCols = ncols * discretizationFactor;
			double fineResolution = resolution / discretizationFactor;
			double fineHalfResolution = fineResolution / 2.0;
			int r, c;
			SimpleGridCell cellValue;
			oldProgress = -1
			for (row = 0; row < fineRows; row++) {
                northing = (north - fineHalfResolution) - (row * fineResolution);
            	for (col = 0; col < fineCols; col++) {
                	easting = (col * fineResolution) + (west + fineHalfResolution);
                    entry = [northing, easting];
                    result = pointsTree.nearestNeighbor(entry);
                    value = (double)result.value;
                    distance = Math.sqrt(result.distance)

                    // see if there are any grid cells within this distance of the point
					results = gridTree.neighborsWithinRange(entry, distance);
					if (results.size() > 0) {
						for (i = 0; i < results.size(); i++) {
                        	cellValue = (SimpleGridCell)results.get(i).value;
                            r = cellValue.row
                            c = cellValue.column
                            if (sumArray.getValue(r, c) != noData) {
                            	sumArray.incrementValue(r, c, value)
                            } else {
                            	sumArray.setValue(r, c, value)
                            }
                            nArray.incrementValue(r, c)
                        }
					}
                }
                if (numFiles == 1) {
                	progress = (int)(100f * row / fineRowsLessOne)
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
			
	
	        try {
	            // create the whitebox header file.
	            FileWriter fw = new FileWriter(outputHeader, false);
	            BufferedWriter bw = new BufferedWriter(fw);
	            PrintWriter out = new PrintWriter(bw, true);
	
	            String str1 = "Min:\t" + Double.toString(Integer.MAX_VALUE);
	            out.println(str1);
	            str1 = "Max:\t" + Double.toString(Integer.MIN_VALUE);
	            out.println(str1);
	            str1 = "North:\t" + Double.toString(north);
	            out.println(str1);
	            str1 = "South:\t" + Double.toString(south);
	            out.println(str1);
	            str1 = "East:\t" + Double.toString(east);
	            out.println(str1);
	            str1 = "West:\t" + Double.toString(west);
	            out.println(str1);
	            str1 = "Cols:\t" + Integer.toString(ncols);
	            out.println(str1);
	            str1 = "Rows:\t" + Integer.toString(nrows);
	            out.println(str1);
	            str1 = "Data Type:\t" + "float";
	            out.println(str1);
	            str1 = "Z Units:\t" + "not specified";
	            out.println(str1);
	            str1 = "XY Units:\t" + "not specified";
	            out.println(str1);
	            str1 = "Projection:\t" + "not specified";
	            out.println(str1);
	            str1 = "Data Scale:\tcontinuous";
	            out.println(str1);
	            if (whatToInterpolate.equals("intensity")) {
	                str1 = "Preferred Palette:\t" + "grey.pal";
	            } else {
	                str1 = "Preferred Palette:\t" + "spectrum.pal";
	            }
	            out.println(str1);
	            str1 = "NoData:\t" + noData;
	            out.println(str1);
	            if (java.nio.ByteOrder.nativeOrder() == java.nio.ByteOrder.LITTLE_ENDIAN) {
	                str1 = "Byte Order:\t" + "LITTLE_ENDIAN";
	            } else {
	                str1 = "Byte Order:\t" + "BIG_ENDIAN";
	            }
	            out.println(str1);
	
	            out.close();

	            // Create the whitebox raster object.
	        	WhiteboxRaster image = new WhiteboxRaster(outputHeader, "rw");
	            oldProgress = -1
	            for (row = 0; row < nrows; row++) {
                    for (col = 0; col < ncols; col++) {
                    	if (nArray.getValue(row, col) > 0) {
                    		image.setValue(row, col, sumArray.getValue(row, col) / nArray.getValue(row, col));
                    	} else {
                    		// use the nearest neighbour
                    		easting = (col * resolution) + (west + halfResolution);
                    		northing = (north - halfResolution) - (row * resolution);
                    		entry = [northing, easting]
							result = pointsTree.nearestNeighbor(entry);
                    		value = (double)result.value;
                    
                    		image.setValue(row, col, value);
                    	}
                    }

                    if (numFiles == 1) {
                    	progress = (int)(100f * row / rowsLessOne)
                    	if (progress > oldProgress) {
                    		oldProgress = progress
                    		pluginHost.updateProgress("Saving data:", progress)
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
                image.addMetadataEntry("Discretization factor: $discretizationFactor");
                if (minVal > Double.NEGATIVE_INFINITY) {
                	image.addMetadataEntry("Min Height: $minVal");
                	image.addMetadataEntry("Max Height: $maxVal");                
                }
                
                image.close();
	
	        } catch (Exception e) {
	            pluginHost.showFeedback(e.getMessage());
	            return Boolean.FALSE;
	        }

			int solved = numSolvedTiles.incrementAndGet() //++;
			progress = (int) (100f * solved / numFiles)
			pluginHost.updateProgress("Interpolated ${solved} tiles:", progress)
	        return Boolean.TRUE;
        }                
    }
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def myClass = new LiDAR_Discrete_Sibson_Interpolation(pluginHost, args, name, descriptiveName)
}
