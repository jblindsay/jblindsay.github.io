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
def name = "FilterUserDefinedWeights"
def descriptiveName = "User-Defined Weights Filter"
def description = "Performs a convolution filter on a raster using a kernel of custom weights"
def toolboxes = ["Filters"]

public class FilterUserDefinedWeights implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public FilterUserDefinedWeights(WhiteboxPluginHost pluginHost, 
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
			def helpFile = "FilterUserDefinedWeights"
			sd.setHelpFile(helpFile)
		
			// Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.
			def pathSep = File.separator
			def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "FilterUserDefinedWeights.groovy"
			sd.setSourceFile(scriptFile)
			
			// add some components to the dialog
			sd.addDialogFile("Input file", "Input Raster File:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Output file", "Output Raster File:", "close", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Kernel file", "Input Kernel File (.txt):", "open", "Text Files (*.txt), txt", true, false)
			def centeringOptions = ["kernel center", "upper-left corner", "upper-right corner", "lower-left corner", "lower-right corner"]
			sd.addDialogComboBox("Center kernel on:", "Center kernel on:", centeringOptions, 0)
			sd.addDialogCheckBox("Reflect values at image edge?", "Reflect values at image edge?", true)
			sd.addDialogCheckBox("Normalize kernel weights?", "Normalize kernel weights?", true)
			sd.addDialogCheckBox("Use parallel processing?", "Use parallel processing?", true)

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
		if (args.length != 7) {
			pluginHost.showFeedback("Incorrect number of arguments given to tool.")
			return
		}
		// read the input parameters
		String inputFile = args[0]
		String outputFile = args[1]
		String kernelfile = args[2]
		String centerOn = args[3]
		boolean relectAtBorders = Boolean.parseBoolean(args[4])
		boolean normalize = Boolean.parseBoolean(args[5])
		boolean processParallel = Boolean.parseBoolean(args[6])

		// I'm sure there's some fancy way of doing this
		// same thing with Groovy closures, but I don't 
		// know how.
		String delimiter = "\t"
		String str = new File(kernelfile).getText()
		String[] str2 = str.split("\n")
		int kernelRows = str2.length
		int kernelCols = str2[0].split(delimiter).length
		if (kernelCols == 1) {
			delimiter = ","
			kernelCols = str2[0].split(delimiter).length
			if (kernelCols == 1) {
				delimiter = " "
				kernelCols = str2[0].split(delimiter).length
				if (kernelCols == 1) {
//					pluginHost.showFeedback("Unrecognized delimeter in kernel file.")
//					return
				}
			}
		}
		
		int kernelCenterX, kernelCenterY
		switch (centerOn) {
			case "kernel center":
				kernelCenterX = (int)(Math.floor(kernelCols / 2.0))
				kernelCenterY = (int)(Math.floor(kernelRows / 2.0))
				break;
			case "upper-left corner":
				kernelCenterX = 0
				kernelCenterY = 0
				break
			case "upper-right corner":
				kernelCenterX = kernelCols
				kernelCenterY = 0
				break
			case "lower-left corner":
				kernelCenterX = 0
				kernelCenterY = kernelRows
				break
			default: // "lower-right corner":
				kernelCenterX = kernelCols
				kernelCenterY = kernelRows
				break
		}

		ArrayList<Double> weights = new ArrayList<>()
		ArrayList<Integer> dX = new ArrayList<>()
		ArrayList<Integer> dY = new ArrayList<>()
		int filterSize = 0
		for (int i = 0; i < str2.length; i++) {
			String[] str3 = str2[i].split(delimiter)
			for (int j = 0; j < kernelCols; j++) {
				if (StringUtilities.isNumeric(str3[j].trim())) {
					weights.add(Double.parseDouble(str3[j].trim()))
					dX.add(j - kernelCenterX)
					dY.add(i - kernelCenterY)
					filterSize++
				}
			}
		}
	
		// set up the input and output rasters
		WhiteboxRaster input = new WhiteboxRaster(inputFile, "r")
		int rows = input.getNumberRows()
		int cols = input.getNumberColumns()
		double nodata = input.getNoDataValue()
		if (relectAtBorders) {
			input.isReflectedAtEdges = true
		}
		
		WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", 
		inputFile, DataType.FLOAT, nodata)

		// filter the image
		int progress = 0
		int oldProgress = -1

		if (processParallel) {
			input.setForceAllDataInMemory(true)

			// this just primes the raster for reading/writing
			input.getValue(0, 0)
			
			pluginHost.updateProgress("Please wait...", 0)
			ArrayList<FilterRow> tasks = new ArrayList<>();
			for (row in 0..(rows - 1)) {
				tasks.add(new FilterRow(input, row, weights, 
				 dX, dY, normalize))
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
        
		} else {
		
			double z, zN, total, totalWeight
			if (!normalize) {
				for (int row = 0; row < rows; row++) {
					for (int col = 0; col < cols; col++) {
						z = input.getValue(row, col)
						if (z != nodata) {
							total = 0
							for (int i = 0; i < filterSize; i++) {
								zN = input.getValue(row + dY.get(i), col + dX.get(i))
								if (zN != nodata) {
									total += zN * weights.get(i)
								}
							}
							output.setValue(row, col, total)
						}
					}
					// update progress bar
					progress = (int)(100f * row / (rows - 1.0))
					if (progress > oldProgress) {
						pluginHost.updateProgress(progress)
						oldProgress = progress
					}
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			} else {
				for (int row = 0; row < rows; row++) {
					for (int col = 0; col < cols; col++) {
						z = input.getValue(row, col)
						if (z != nodata) {
							total = 0
							totalWeight = 0
							for (int i = 0; i < filterSize; i++) {
								zN = input.getValue(row + dY.get(i), col + dX.get(i))
								if (zN != nodata) {
									total += zN * weights.get(i)
									totalWeight += weights.get(i)
								}
							}
							output.setValue(row, col, total / totalWeight)
						}
					}
					// update progress bar
					progress = (int)(100f * row / (rows - 1.0))
					if (progress > oldProgress) {
						pluginHost.updateProgress(progress)
						oldProgress = progress
					}
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}			
			}
		
		}
		
		input.close()

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

    public List<Future<RowNumberAndData>> getExecutorResults(ExecutorService executor, ArrayList<FilterRow> tasks) {
    	List<Future<RowNumberAndData>> results = executor.invokeAll(tasks);
		return results
    }

	@CompileStatic
    class FilterRow implements Callable<RowNumberAndData> {
		private final int row
	    private final ArrayList<Double> weights
	    private final ArrayList<Integer> dX
	    private final ArrayList<Integer> dY
	    private final WhiteboxRaster raster
	    private final boolean normalize
    	FilterRow(WhiteboxRaster raster, int row, 
    	  ArrayList<Double> weights,  ArrayList<Integer> dX, 
    	  ArrayList<Integer> dY, boolean normalize) {
        	this.raster = raster
        	this.row = row
            this.weights = weights
            this.dX = dX
       	    this.dY = dY
       	    this.normalize = normalize
       	}
        	
        @Override
	    public RowNumberAndData call() {
	       	int cols = raster.getNumberColumns()
	       	int filterSize = dX.size()
	       	double nodata = raster.getNoDataValue()
            double[] retData = raster.getRowValues(row)
   	        double total, zN, z, totalWeight
       	    for (int col = 0; col < cols; col++) {
           	    z = retData[col]
           	    if (z != nodata) {
					total = 0
					totalWeight = 0
					for (int i = 0; i < filterSize; i++) {
						zN = raster.getValue(row + dY.get(i), col + dX.get(i))
						if (zN != nodata) {
							total += zN * weights.get(i)
							totalWeight += weights.get(i)
						}
					}
					if (normalize) {
						retData[col] = total / totalWeight
					} else {
						retData[col] = total
					}
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
	def f = new FilterUserDefinedWeights(pluginHost, args, descriptiveName)
}
