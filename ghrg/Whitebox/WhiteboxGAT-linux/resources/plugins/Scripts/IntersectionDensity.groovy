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
import java.util.concurrent.Future
import java.util.concurrent.*
import java.util.Date
import java.util.ArrayList
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import whitebox.structures.KdTree
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "IntersectionDensity"
def descriptiveName = "Intersection Density"
def description = "Calculates the spatial pattern of polygon intersection with a buffer distance"
def toolboxes = ["VectorTools"]

public class IntersectionDensity implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
	
    public IntersectionDensity(WhiteboxPluginHost pluginHost, 
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
            def helpFile = "IntersectionDensity"
            sd.setHelpFile(helpFile)
		
            // Specifying the source file allows the 'view code' 
            // button on the tool dialog to be displayed.
            def pathSep = File.separator
            def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "IntersectionDensity.groovy"
            sd.setSourceFile(scriptFile)
			
            // add some components to the dialog
            sd.addDialogFile("Input file", "Input Vector Polygon File:", "open", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogFile("Output file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
            sd.addDialogDataInput("Buffer Distance:", "Enter a buffer distance", "", true, false)
            sd.addDialogDataInput("Output Grid Resolution:", "Enter a grid resolution", "", true, false)
			
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
            int i, f, progress, oldProgress
	  		double x1, x2, y1, y2, dist1, dist2
            if (args.length != 4) {
                pluginHost.showFeedback("Incorrect number of arguments given to tool.")
                return
            }
            // read the input parameters
            String inputFile = args[0]
            String outputFile = args[1]
            double neighbourhoodRadius = Double.parseDouble(args[2])
            double radiusSquared = neighbourhoodRadius * neighbourhoodRadius
            double cellSize = Double.parseDouble(args[3])

            def input = new ShapeFile(inputFile)

            // make sure that input is of a POLYLINE or POLYGON base shapetype
            if (input.getShapeType().getBaseType() != ShapeType.POLYLINE
                && input.getShapeType().getBaseType() != ShapeType.POLYGON) {
                pluginHost.showFeedback("Input shapefile must be of a POLYLINE or POLYGON base shapetype.")
                return
            }

            double[][] points
           
            double x
            double y
            double z = 0.0
            List<KdTree.Entry<Double>> results
            int recNum = 0
            int numFeatures = input.getNumberOfRecords()
            
			// find out how many points there are to add to the kd-tree
			def numPoints = 0
			progress = 0
	    	oldProgress = -1
	    	for (ShapeFileRecord record : input.records) {
                recNum++
                f = record.getRecordNumber() - 1
                
                points = record.getGeometry().getPoints()

                numPoints += record.getGeometry().getPoints().length
                
                // update progress
                progress = (int)(100.0 * recNum / numFeatures)
				if (progress != oldProgress) {
                    oldProgress = progress
                    pluginHost.updateProgress("Loop 1 of 3:", progress)
                }
                // check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
            }
            
            KdTree<Double> pointsTree = new KdTree.SqrEuclid<Double>(2, new Integer(numPoints))
            
            progress = 0
	    	oldProgress = -1
	    	for (ShapeFileRecord record : input.records) {
                recNum++
                f = record.getRecordNumber() - 1
                
                points = record.getGeometry().getPoints()
                
                for (i = 0; i < points.length; i++) {
                	def entry = new double[2]
                	entry[1] = points[i][0]
                    entry[0] = points[i][1]
                    
					pointsTree.addPoint(entry, (double)f)
                }

                // update progress
                progress = (int)(100.0 * recNum / numFeatures)
				if (progress != oldProgress) {
                    oldProgress = progress
                    pluginHost.updateProgress("Loop 2 of 3:", progress)
                }
                // check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
            }

            // create the output raster
            def rows = (int) (Math.ceil((input.getyMax() - input.getyMin()) / cellSize))
            def cols = (int) (Math.ceil((input.getxMax() - input.getxMin()) / cellSize))
            def north = input.getyMax()
            def south = input.getyMax() - rows * cellSize
            def east = input.getxMin() + cols * cellSize
            def west = input.getxMin()
            double noData = -32768.0
            def dataType = WhiteboxRasterBase.DataType.FLOAT
            
            def output = new WhiteboxRaster(outputFile, north, south, east, west,
                rows, cols, WhiteboxRasterBase.DataScale.CONTINUOUS,
                dataType, 0.0, noData);
            output.setPreferredPalette("spectrum.plt")





//			double[] xCoords = new double[cols]
//			for (int col = 0; col < cols; col++) {
//               xCoords[col] = output.getXCoordinateFromColumn(col)
//			}
//			
//			// this just primes the raster for reading/writing
//			output.getValue(0, 0)
//			
//			pluginHost.updateProgress("Please wait...", 0)
//			ArrayList<DoWorkOnRow> tasks = new ArrayList<>();
//			for (row in 0..(rows - 1)) {
//				y = output.getYCoordinateFromRow(row)
//				tasks.add(new DoWorkOnRow(input, pointsTree, 
//    	    	y, xCoords, row, neighbourhoodRadius))
//			}
//        
//			ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
//	  	    // the only reason for the getExecutorResults method 
//	  	    // is that Groovy throws a compilation type mis-match
//	  	    // error when compiled statically. I think it's a bug.
//	  	    List<Future<RowNumberAndData>> exResults = getExecutorResults(executor, tasks);
//	        executor.shutdown();
//
//			i = 0
//    	    for (Future<RowNumberAndData> result : exResults) {
//    	    	RowNumberAndData data = result.get()
//    	    	if (data != null) {
//        		def row = data.getRow()
//        		output.setRowValues(row, data.getData())
//    	    	}
//            	i++
//				// update progress bar
//				progress = (int)(100f * i / rows)
//				if (progress > oldProgress) {
//					pluginHost.updateProgress("Progress", progress)
//					oldProgress = progress
//				}
//				// check to see if the user has requested a cancellation
//				if (pluginHost.isRequestForOperationCancelSet()) {
//					pluginHost.showFeedback("Operation cancelled")
//					return
//				}
//    	    }





    	    
            // interpolate into the output raster
            def numNeighbouringFeatures = 0
            oldProgress = -1
            for (int row = 0; row < rows; row++) {
                for (int col = 0; col < cols; col++) {
                    x = output.getXCoordinateFromColumn(col)
                    y = output.getYCoordinateFromRow(row)
                    double[] entry = new double[2]
                    entry[0] = y
                    entry[1] = x
                    results = pointsTree.neighborsWithinRange(entry, neighbourhoodRadius)
					if (results.size() > 0) {
	                    def neighbouringFeatures = new boolean[numFeatures]
    	                for (KdTree.Entry entry2 : results) {
        	                i = (int)entry2.value
            	            neighbouringFeatures[i] = true
                	    }
                    
	                    int intersections = 0
    	                for (f = 0; f < numFeatures; f++) {
        	                if (neighbouringFeatures[f]) {
            	            	ShapeFileRecord record = input.getRecord(f)
								points = record.getGeometry().getPoints()
                	
                				// scan each line segment for an intersection
            	    			for (i = 1; i < points.length; i++) {
	        	        			x1 = points[i - 1][0]
    		                		y1 = points[i - 1][1]
	    	                		x2 = points[i][0]
            	        			y2 = points[i][1]
            	        			if (lineSegmentInstersectsCircle(x1, 
            	        			     y1, x2, y2, x, y, neighbourhoodRadius)) {
            	        				intersections++
            	        			}
		                		}
        	                }
            	        }
                	    
                    	output.setValue(row, col, intersections)
                    
					}
					
                }
                // update progress
                progress = (int) (100f * row / (rows - 1))
                if (progress != oldProgress) {
                    oldProgress = progress
                    pluginHost.updateProgress("Loop 3 of 3:", progress)
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
            output.close()

            // display the output image
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

	
    public List<Future<RowNumberAndData>> getExecutorResults(ExecutorService executor, ArrayList<DoWorkOnRow> tasks) {
    	List<Future<RowNumberAndData>> results = executor.invokeAll(tasks);
		return results
    }

	@CompileStatic
    class DoWorkOnRow implements Callable<RowNumberAndData> {
		private final int row
	    private final WhiteboxRaster raster
	    private final KdTree<Double> pointsTree
	    private final ShapeFile input
	    private double radius
	    private double y
	    private double[] x
	    
    	DoWorkOnRow(ShapeFile input, KdTree<Double> pointsTree, 
    	    double y, double[] x,  int row, double radius) {
        	this.raster = raster
        	this.row = row
        	this.input = input
        	this.pointsTree = pointsTree
        	this.radius = radius
			this.y = y
			this.x = x
       	}
        	
        @Override
	    public RowNumberAndData call() {
	    	try {
	       	int cols = x.length
//	       	double nodata = raster.getNoDataValue()
            double[] retData = new double[cols]
   	        double z, x1, y1, x2, y2
   	        int i, f
			double[][] points
   	        int numFeatures = input.getNumberOfRecords()
   	        
//   	        List<KdTree.Entry<Double>> results

	        for (int col = 0; col < cols; col++) {
	        	double[] entry = new double[2]
                entry[0] = y
                entry[1] = x[col]
                List<KdTree.Entry<Double>> results = pointsTree.neighborsWithinRange(entry, radius)
                retData[col] = results.size()
//				if (results.size() > 0) {
//	                def neighbouringFeatures = new boolean[numFeatures]
//    	            for (KdTree.Entry entry2 : results) {
//      	                i = (int)entry2.value
//           	            neighbouringFeatures[i] = true
//               	    }
//                    
//	                int intersections = 0
//    	            for (f = 0; f < numFeatures; f++) {
//       	                if (neighbouringFeatures[f]) {
//           	            	ShapeFileRecord record = input.getRecord(f)
//							points = record.getGeometry().getPoints()
//                	
//              				// scan each line segment for an intersection
//           	    			for (i = 1; i < points.length; i++) {
//	       	        			x1 = points[i - 1][0]
//    		                	y1 = points[i - 1][1]
//	    	                	x2 = points[i][0]
//            	        		y2 = points[i][1]
//            	       			if (lineSegmentInstersectsCircle(x1, 
//            	        		    y1, x2, y2, x[col], y, radius)) {
//            	        			intersections++
//            	        		}
//		               		}
//        	            }
//            	    }
//                	    
//                    retData[col] = intersections
//                    
//				}
					
            }

	        def ret = new RowNumberAndData(row, retData)
	        return ret
	    	} catch (Exception e) {
//	    		println row;
	    		//println e.printStackTrace()
	    		return null
	    	}
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
    
    @CompileStatic
    public static boolean lineSegmentInstersectsCircle(double x1, double y1, 
    	double x2, double y2, double cX, double cY, double r) {
    	double a = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y1 - y1);
		double b  = 2.0 * ((x2 - x1) * (x1 - cX) +(y2 - y1) * (y1 - cY));
		double cc = cX * cX + cY * cY + x1 * x1 + y1 * y1 - 2.0 * (cX * x1 + cY * y1) - r * r;
		double deter = b * b - 4.0 * a * cc;
		if (deter <= 0 ) {
			return false;
		} else {
			double e = Math.sqrt(deter);
			double u1 = ( - b + e ) / (2 * a );
			double u2 = ( - b - e ) / (2 * a );
			if ((u1 < 0 || u1 > 1) && (u2 < 0 || u2 > 1)) {
				return false;
			} else {
				return true;
			}
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
    def f = new IntersectionDensity(pluginHost, args, descriptiveName)
}
