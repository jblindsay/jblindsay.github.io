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
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "CorrectVignetting"
def descriptiveName = "Correct Vignetting"
def description = "Corrects the darkening of images towards corners."
def toolboxes = ["Photogrammetry", "ImageEnhancement"]

public class CorrectVignetting implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd
	private String descriptiveName
	
	public CorrectVignetting(WhiteboxPluginHost pluginHost, 
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
			def helpFile = "CorrectVignetting"
			sd.setHelpFile(helpFile)
		
			// Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.
			def pathSep = File.separator
			def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "CorrectVignetting.groovy"
			sd.setSourceFile(scriptFile)
			
			// add some components to the dialog
			sd.addDialogFile("Input image file", "Input Image:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Input principal point file", "Input Principal Point Vector:", "open", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogFile("Output file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogDataInput("Enter the camera focal length in millimetres.", "Camera Focal Length (mm):", "304.8", true, false)
            sd.addDialogDataInput("Enter the distance between the left and right edges in millimetres.", "Distance Between Left-Right Edges (mm):", "228.6", true, false)
            sd.addDialogDataInput("<html>The <i>n</i> parameter controls the effect gradient.</html>", "<html><i>n</i> Parameter:</html>", "4.0", true, false)
            
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
	  	
			if (args.length != 6) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String inputFile = args[0]
			String inputPPFile = args[1]
			String outputFile = args[2]
			double focalLength = Double.parseDouble(args[3])
			double distBetweenEdges = Double.parseDouble(args[4])
			double n = Double.parseDouble(args[5])

			// read the input image and PP vector files
			WhiteboxRasterInfo image = new WhiteboxRasterInfo(inputFile)
			double nodata = image.getNoDataValue()
			int rows = image.getNumberRows()
			int cols = image.getNumberColumns()
			double scaleFactor = distBetweenEdges / cols
			
			def inputPP = new ShapeFile(inputPPFile)
            ShapeType shapeType = inputPP.getShapeType()
			if (shapeType.getBaseType() != ShapeType.POINT) {
				pluginHost.showFeedback("The principal point file does not appear to be of a POINT shape type.")
				return;
			}

			// get the x and y coordinates of the PP in pixel values.
			ShapeFileRecord record = inputPP.getRecord(0)
			double[][] point = record.getGeometry().getPoints()
			double PPx = image.getColumnFromXCoordinate(point[0][0])
			double PPy = image.getRowFromYCoordinate(point[0][1])
			inputPP = null
			
			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", 
  		  	  inputFile, DataType.FLOAT, nodata)


  		  	int progress = 0
			int oldProgress = -1
				
			pluginHost.updateProgress("Loop 1 of 2, Please wait...", 0)
			ArrayList<DoWorkOnRow> tasks = new ArrayList<>();
			for (row in 0..(rows - 1)) {
				double[] data = image.getRowValues(row)
				tasks.add(new DoWorkOnRow(data, nodata, row, scaleFactor,
	              n, PPx, PPy, focalLength))
			}
	    
			ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
	  	    // the only reason for the getExecutorResults method 
	  	    // is that Groovy throws a compilation type mis-match
	  	    // error when compiled statically. I think it's a bug.
	  	    List<Future<RowNumberAndData>> results = getExecutorResults(executor, tasks); 
	        executor.shutdown();
	
	        int i = 0
		    for (Future<RowNumberAndData> result : results) {
	    		RowNumberAndData data = result.get()
	    		def row = data.getRow()
	    		output.setRowValues(row, data.getData())
	        	i++
				// update progress bar
				progress = (int)(100f * i / rows)
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

			output.findMinAndMaxVals()
			
		    double minValOut = output.getMinimumValue()
		    double maxValOut = output.getMaximumValue()
		    double minValIn = image.getMinimumValue()
		    double maxValIn = image.getMaximumValue()

			pluginHost.updateProgress("Loop 2 of 2, Please wait...", 0)
			ArrayList<DoWorkOnRow2> tasks2 = new ArrayList<>();
			for (row in 0..(rows - 1)) {
				double[] data = output.getRowValues(row)
				tasks2.add(new DoWorkOnRow2(data, nodata, row, minValOut,
	              maxValOut, minValIn, maxValIn))
			}
	    
			ExecutorService executor2 = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
	  	    // the only reason for the getExecutorResults method 
	  	    // is that Groovy throws a compilation type mis-match
	  	    // error when compiled statically. I think it's a bug.
	  	    List<Future<RowNumberAndData>> results2 = getExecutorResults2(executor2, tasks2); 
	        executor2.shutdown();
	
	        i = 0
	        oldProgress = -1
		    for (Future<RowNumberAndData> result : results2) {
	    		RowNumberAndData data = result.get()
	    		def row = data.getRow()
	    		output.setRowValues(row, data.getData())
	        	i++
				// update progress bar
				progress = (int)(100f * i / rows)
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

		    output.setDisplayMinimum(image.getDisplayMinimum())
		    output.setDisplayMaximum(image.getDisplayMaximum())

			image.close()

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

    public List<Future<RowNumberAndData>> getExecutorResults(ExecutorService executor, ArrayList<DoWorkOnRow> tasks) {
    	List<Future<RowNumberAndData>> results = executor.invokeAll(tasks);
		return results
    }

	@CompileStatic
    class DoWorkOnRow implements Callable<RowNumberAndData> {
		private double[] data
		private int row
		private double nodata
		private double scaleFactor
		private double n
		private double PPx
		private double PPy
		private double focalLength
		
	    DoWorkOnRow(double[] data, double nodata, int row, double scaleFactor,
	      double n, double PPx, double PPy, double focalLength) {
        	this.data = data
        	this.row = row
            this.nodata = nodata
            this.scaleFactor = scaleFactor
            this.n = n
            this.PPx = PPx
            this.PPy = PPy
            this.focalLength = focalLength
       	}
        	
        @Override
	    public RowNumberAndData call() {
	    	// check to see if the user has requested a cancellation
			if (pluginHost.isRequestForOperationCancelSet()) {
				pluginHost.showFeedback("Operation cancelled")
				return
			}
			
	    	int cols = data.length
	       	double[] retData = new double[cols]
			double dist
			double theta
			double value
			
			for (int col = 0; col < cols; col++) {
  				if (data[col] != nodata) {
  					dist = Math.sqrt((row - PPy) * (row - PPy) + (col - PPx) * (col - PPx))
  					theta = Math.atan((dist * scaleFactor / focalLength))
 					value = Math.pow(Math.cos(theta), n)
  					retData[col] = data[col] / value
				} else {
  					retData[col] = nodata
  				}
  			}
  			
	        def ret = new RowNumberAndData(row, retData)
	        
    	    return ret
        }                
    }

	public List<Future<RowNumberAndData>> getExecutorResults2(ExecutorService executor, ArrayList<DoWorkOnRow2> tasks) {
    	List<Future<RowNumberAndData>> results = executor.invokeAll(tasks);
		return results
    }
    
	@CompileStatic
    class DoWorkOnRow2 implements Callable<RowNumberAndData> {
		private double[] data
		private int row
		private double nodata
		private double minValOut
		private double maxValOut
		private double minValIn
		private double maxValIn
		private double rangeOut
		private double rangeIn
		
	    DoWorkOnRow2(double[] data, double nodata, int row, double minValOut,
	      double maxValOut, double minValIn, double maxValIn) {
        	this.data = data
        	this.row = row
            this.nodata = nodata
            this.minValOut = minValOut
            this.maxValOut = maxValOut
            this.minValIn = minValIn
            this.maxValIn = maxValIn
            this.rangeOut = maxValOut - minValOut
            this.rangeIn = maxValIn - minValIn
       	}
        	
        @Override
	    public RowNumberAndData call() {
	    	// check to see if the user has requested a cancellation
			if (pluginHost.isRequestForOperationCancelSet()) {
				pluginHost.showFeedback("Operation cancelled")
				return
			}
			
	    	int cols = data.length
	       	double[] retData = new double[cols]
			double value
			
			for (int col = 0; col < cols; col++) {
  				if (data[col] != nodata) {
  					value = minValIn + (data[col] - minValOut) / rangeOut * rangeIn
  					retData[col] =  value
				} else {
  					retData[col] = nodata
  				}
  			}
  			
	        def ret = new RowNumberAndData(row, retData)
	        
    	    return ret
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
	def f = new CorrectVignetting(pluginHost, args, descriptiveName)
}
