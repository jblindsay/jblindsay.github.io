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
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "AverageOverlay"
def descriptiveName = "Average Overlay"
def description = "Calculates the average for each grid cell from a group of raster images."
def toolboxes = ["OverlayTools", "StatisticalTools"]

public class AverageOverlay implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public AverageOverlay(WhiteboxPluginHost pluginHost, 
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
			sd.addDialogMultiFile("Select some input raster files", "Input Raster Files:", "Raster Files (*.dep), DEP")
			sd.addDialogFile("Output file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
			
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
	  	
		if (args.length != 2) {
			pluginHost.showFeedback("Incorrect number of arguments given to tool.")
			return
		}
		// read the input parameters
		String inputFileString = args[0]
		String outputFile = args[1]
		
		String[] inputFiles = inputFileString.split(";")
		
		ArrayList<WhiteboxRaster> rasters = new ArrayList<>()
		int rows, cols
		int numFiles = 0
		double[] nodataVals = new double[inputFiles.length]

		
		for (String inputFile : inputFiles) {
			if (!inputFile.isEmpty()) {
				rasters.add(new WhiteboxRaster(inputFile, "r"))
				if (numFiles == 0) {
					rows = rasters.get(0).getNumberRows()
					cols = rasters.get(0).getNumberColumns()
				} else {
					if (rasters.get(numFiles).getNumberRows() != rows ||
					 rasters.get(numFiles).getNumberColumns() != cols) {
					 	pluginHost.showFeedback("Each input raster must have the same dimensions.")
					 	return
					 }
				}
				nodataVals[numFiles] = rasters.get(numFiles).getNoDataValue()
				numFiles++
			}
		}

		WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", 
		rasters.get(0).getHeaderFile(), DataType.FLOAT, nodataVals[0])

		int progress = 0
		int oldProgress = -1
			
		pluginHost.updateProgress("Please wait...", 0)
		ArrayList<DoWorkOnRow> tasks = new ArrayList<>();
		for (row in 0..(rows - 1)) {
			double[][] data = new double[numFiles][]
			for (int j = 0; j < numFiles; j++) {
				data[j] = rasters.get(j).getRowValues(row)
			}
			tasks.add(new DoWorkOnRow(data, nodataVals, row))
		}
    
		ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
  	    // the only reason for the getExecutorResults method 
  	    // is that Groovy throws a compilation type mis-match
  	    // error when compiled statically. I think it's a bug.
  	    List<Future<RowNumberAndData>> results = getExecutorResults(executor, tasks); //executor.invokeAll(tasks);
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
        
		for (int j = 0; j < numFiles; j++) {
			rasters.get(j).close()
		}

		output.addMetadataEntry("Created by the "
                    + descriptiveName + " tool.")
        output.addMetadataEntry("Created on " + new Date())
		output.close()

		// display the output image
		pluginHost.returnData(outputFile)

		// reset the progress bar
		pluginHost.updateProgress(0)
	  } catch (OutOfMemoryError oe) {
            pluginHost.showFeedback("An out-of-memory error has occurred during operation.")
      } catch (Exception e) {
            pluginHost.showFeedback("An error has occurred during operation. See log file for details.")
            pluginHost.logException("Error in " + descriptiveName, e)
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
		private double[][] data
		private int row
		private double[] nodata
		
	    DoWorkOnRow(double[][] data, double[] nodata, int row) {
        	this.data = data
        	this.row = row
            this.nodata = nodata
       	}
        	
        @Override
	    public RowNumberAndData call() {
	    	int numImages = data.length
	       	int cols = data[0].length
	       	double[] retData = new double[cols]
   	        double total, n
       	    for (int col = 0; col < cols; col++) {
       	    	n = 0
       	    	total = 0
       	    	for (int i = 0; i < numImages; i++) {
       	    		if (data[i][col] != nodata[i]) {
       	    			n++
       	    			total += data[i][col]
       	    		}
       	    	}
       	    	if (n > 0) {
           	    	retData[col] = total / n
       	    	} else {
       	    		retData[col] = nodata[0]
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
	def f = new AverageOverlay(pluginHost, args, name, descriptiveName)
}
