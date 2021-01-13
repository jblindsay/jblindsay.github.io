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
import groovy.time.TimeDuration
import groovy.time.TimeCategory

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "ElevationPercentileImageHisto"
def descriptiveName = "Elevation Percentile (Image Histogram)"
def description = "Calculates elevation percentile for a DEM."
def toolboxes = ["ElevResiduals"]

public class ElevationPercentileImageHisto implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public ElevationPercentileImageHisto(WhiteboxPluginHost pluginHost, 
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
			def helpFile = "ElevationPercentileImageHisto"
			sd.setHelpFile(helpFile)
		
			// Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.
			def pathSep = File.separator
			def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "ElevationPercentileImageHisto.groovy"
			sd.setSourceFile(scriptFile)
			
			// add some components to the dialog
			sd.addDialogFile("Input DEM file", "Input DEM:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Output file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogDataInput("Search neighbourhood radius (cells)", "Search Neighbourhood Radius (cells):", "", true, false)
            sd.addDialogDataInput("Number of histogram bins", "Number of Histogram Bins:", "16", true, false)
            
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
			Date start = new Date()
			
			int row, col, bin, progress, oldProgress, N, numLess, numInBin
			int x1, x2, y1, y2
			double z, percentile
			if (args.length != 4) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String inputFile = args[0]
			String outputFile = args[1]
			int neighbourhoodSize = Integer.parseInt(args[2])
			if (neighbourhoodSize < 1) { neighbourhoodSize = 1 }
			int numBins = Integer.parseInt(args[3])
			if (numBins < 2) { numBins = 2 }
			if (numBins > 1024) { numBins = 1024 }
			int numCells
			
			// read the input image
			WhiteboxRaster image = new WhiteboxRaster(inputFile, "r")
			double nodata = image.getNoDataValue()
			int rows = image.getNumberRows()
			int cols = image.getNumberColumns()
			double minValue = image.getMinimumValue()
			double maxValue = image.getMaximumValue()
			double range = maxValue - minValue
			double binSize = Math.ceil(range / numBins)

			double[] binLowerValue = new double[numBins]
			for (int i = 0; i < numBins; i++) {
				binLowerValue[i] = minValue + i * binSize
			}
			
			Histogram[][] histoImage = new Histogram[rows][cols]
			oldProgress = -1
			for (row = 0; row < rows; row++) {
				Histogram rowSum = new Histogram(numBins)
				
				for (col = 0; col < cols; col++) {
					histoImage[row][col] = new Histogram(numBins)
  					
					z = image.getValue(row, col)
	  				if (z != nodata) {
  						bin = (int)(Math.floor((z - minValue) / binSize))
						
  						rowSum.incrementBin(bin)
  					}

  					if (row > 0) {
  						histoImage[row][col].setHisto(addHistograms(rowSum.histo, (histoImage[row - 1][col]).histo))
  					} else {
  						histoImage[row][col].setHisto(rowSum.histo)
  					}
				}
				progress = (int)(100f * row / (rows - 1))
				if (progress != oldProgress) {
					pluginHost.updateProgress("Loop 1 of 2:", progress)
					oldProgress = progress
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}
			
			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", 
  		  	  inputFile, DataType.FLOAT, nodata)
  		  	output.setNoDataValue(nodata)
			output.setPreferredPalette("blue_white_red.plt")

			oldProgress = -1
			for (row = 0; row < rows; row++) {
				y1 = row - neighbourhoodSize
				if (y1 < 0) { y1 = 0 }
				if (y1 >= rows) { y1 = rows - 1 }

				y2 = row + neighbourhoodSize
				if (y2 < 0) { y2 = 0 }
				if (y2 >= rows) { y2 = rows - 1 }
				
				for (col = 0; col < cols; col++) {
	  				z = image.getValue(row, col)
	  				if (z != nodata) {
  						bin = (int)(Math.floor((z - minValue) / binSize))
						
	  					x1 = col - neighbourhoodSize
						if (x1 < 0) { x1 = 0 }
						if (x1 >= cols) { x1 = cols - 1 }
	
						x2 = col + neighbourhoodSize
						if (x2 < 0) { x2 = 0 }
						if (x2 >= cols) { x2 = cols - 1 }

						int[] a = (histoImage[y2][x2]).histo
						int[] b = (histoImage[y1][x1]).histo
						int[] c = (histoImage[y1][x2]).histo
						int[] d = (histoImage[y2][x1]).histo

						int[] e = addHistograms(a, b)
						int[] f = subtractHistograms(e, c)
						int[] g = subtractHistograms(f, d)

						N = 0
						numLess = 0
						for (int i = 0; i < numBins; i++) {
							N += g[i]
							if (i < bin) {
								numLess += g[i]
							}
						}

						if (N > 0) {
							percentile = 100.0 * (numLess + (z - binLowerValue[bin]) / binSize * g[bin]) / N
							output.setValue(row, col, percentile)
						}
  					}
  				}
  				progress = (int)(100f * row / (rows - 1))
				if (progress > oldProgress) {
					pluginHost.updateProgress("Loop 2 of 2:", progress)
					oldProgress = progress
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}

			image.close()

			output.flush()
			output.findMinAndMaxVals()
			output.setDisplayMaximum(0)
			output.setDisplayMinimum(100)
			output.addMetadataEntry("Created by the " + descriptiveName + " tool.")
	        output.addMetadataEntry("Created on " + new Date())
	        output.addMetadataEntry("Window size (D): ${neighbourhoodSize * 2 + 1}")
			output.addMetadataEntry("Histogram bins: $numBins")
			Date stop = new Date()
			TimeDuration td = TimeCategory.minus(stop, start)
			output.addMetadataEntry("Elapsed time: $td")
			
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

	@CompileStatic
    class Histogram {
    	int[] histo 
//    	int numBins

    	public Histogram() {
    		
    	}

    	@CompileStatic
    	public Histogram(int numBins) {
    		this.histo = new int[numBins]
//    		this.numBins = numBins
    	}
//
//		@CompileStatic
//		public void setNumBins(int numBins) {
//			this.numBins = numBins
//			this.histo = new int[numBins]
//		}

		@CompileStatic
    	public int[] getHisto() {
    		return histo
    	}

		@CompileStatic
    	public void setHisto(int[] data) {
//    		this.numBins = data.length
    		this.histo = new int[data.length]
    		System.arraycopy(data, 0, this.histo, 0, data.length);
    	}

		@CompileStatic
    	public void incrementBin(int bin) {
    		this.histo[bin]++
    	}

		@CompileStatic
    	public void decrementBin(int bin) {
    		this.histo[bin]--
    	}
    }

	@CompileStatic
    private int[] subtractHistograms(int[] a, int[] b) {
//    	int numBins = a.length
// 		if (b.length != numBins) {
// 			pluginHost.showFeedback("Histogram lengths are not identical")
// 		}
    	int[] result = new int[a.length]
    	for (int i = 0; i < a.length; i++) {
    		result[i] = a[i] - b[i]
    	}
    	return result
    }

	@CompileStatic
    private int[] addHistograms(int[] a, int[] b) {
//    	int numBins = a.length
// 		if (b.length != numBins) {
// 			pluginHost.showFeedback("Histogram lengths are not identical")
// 		}
    	int[] result = new int[a.length]
    	for (int i = 0; i < a.length; i++) {
    		result[i] = a[i] + b[i]
    	}
    	return result
    }
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def f = new ElevationPercentileImageHisto(pluginHost, args, descriptiveName)
}
