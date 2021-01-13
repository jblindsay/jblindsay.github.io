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
import java.io.File
import java.util.Date
import java.util.ArrayList
import java.util.Arrays
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.FileUtilities;
import whitebox.utilities.Topology;
import whitebox.interfaces.InteropPlugin.InteropPluginType
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "ImportArcAsciiGrid"
def descriptiveName = "Import ArcGIS ASCII Grid"
def description = "Imports an ArcGIS ASCII grid file (.txt, .asc)."
def toolboxes = ["IOTools"]

// The following variables tell the plugin host that 
// this tool should be used as a supported geospatial
// data format.
def extensions = ["txt", "asc"]
def fileTypeName = "ArcGIS ASCII Grid"
def isRasterFormat = true
def interopPluginType = InteropPluginType.importPlugin;

public class ImportArcAsciiGrid implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
	
    public ImportArcAsciiGrid(WhiteboxPluginHost pluginHost, 
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
            sd.addDialogMultiFile("Select the input ArcGIS ASCII grid files", "Input ArcGIS ASCII Raster Grid Files:", "ArcGIS ASCII Grid (*.txt; *.asc), TXT, ASC")
			
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
			int i = 0;
        	double x, y
        	int progress = 0;
        	int oldProgress = -1;

        	if (args.length != 1) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String inputFileString = args[0]
			if (inputFileString.isEmpty()) {
	            pluginHost.showFeedback("One or more of the input parameters have not been set properly.");
	            return;
	        }
			String[] inputFiles = inputFileString.split(";")
			int numFiles = inputFiles.length;

			
            for (i = 0; i < numFiles; i++) {
                progress = (int) (100f * i / (numFiles - 1));
                if (numFiles > 1) {
                	pluginHost.updateProgress("Loop " + (i + 1) + " of " + numFiles + ":", progress);
                }
                String arcFileName = inputFiles[i];
                File arcFile = new File(arcFileName)
                if (!arcFile.exists()) {
                    pluginHost.showFeedback("ArcGIS ASCII file does not exist.");
                    break;
                }				
                
				// set up the output file
				String inputFileExtension = FileUtilities.getFileExtension(arcFileName)
				String outputFile = arcFileName.replace(".${inputFileExtension}", ".dep")
				
				int rows, cols
				double north, south, east, west
				double xllcorner, yllcorner, xllcenter, yllcenter
				double cellsize, nodata, z
				String delimiter = " "
                int row = 0
                int col = 0
                oldProgress = -1
				WhiteboxRaster output;
				arcFile.eachLine { line -> 
					if (pluginHost.isRequestForOperationCancelSet()) {
						return
					}
					String str = ((String)(line)).trim().toLowerCase()
					if (str.startsWith("ncols")) {
						cols = Integer.parseInt(str.replace("ncols", "").trim())
					} else if (str.startsWith("nrows")) {
						rows = Integer.parseInt(str.replace("nrows", "").trim())
					} else if (str.startsWith("xllcorner")) {
						xllcorner = Double.parseDouble(str.replace("xllcorner", "").trim())
					} else if (str.startsWith("yllcorner")) {
						yllcorner = Double.parseDouble(str.replace("yllcorner", "").trim())
					} else if (str.startsWith("xllcenter")) {
						xllcenter = Double.parseDouble(str.replace("xllcenter", "").trim())
					} else if (str.startsWith("yllcenter")) {
						yllcenter = Double.parseDouble(str.replace("yllcenter", "").trim())
					} else if (str.startsWith("cellsize")) {
						cellsize = Double.parseDouble(str.replace("cellsize", "").trim())
						if (xllcorner != null) {
                            east = xllcorner + cols * cellsize;
                            west = xllcorner;
                            south = yllcorner;
                            north = yllcorner + rows * cellsize;
                        } else {
                            east = xllcenter - (0.5 * cellsize) + cols * cellsize;
                            west = xllcenter - (0.5 * cellsize);
                            south = yllcenter - (0.5 * cellsize);
                            north = yllcenter - (0.5 * cellsize) + rows * cellsize;
                        }
					} else if (str.startsWith("nodata_value")) {
						nodata = Double.parseDouble(str.replace("nodata_value", "").trim())
					} else {
						// it's a data line.
						if (output == null) {
							output = new WhiteboxRaster(outputFile, north, south, east, west,
                        		rows, cols, WhiteboxRasterBase.DataScale.CONTINUOUS,
                        		WhiteboxRasterBase.DataType.FLOAT, 0.0d, nodata);
                        	output.setPreferredPalette("grey.pal")
						}
						if (!str.isEmpty()) {
							String[] data = str.split(delimiter);
	                        if (data.length <= 1) {
	                            delimiter = "\t";
	                            data = str.split(delimiter);
	                            if (data.length <= 1) {
	                                delimiter = " ";
	                                data = str.split(delimiter);
	                                if (data.length <= 1) {
	                                    delimiter = ",";
	                                    data = str.split(delimiter);
	                                }
	                            }
	                        }
	                        for (int d = 0; d < data.length; d++) {
	                        	String str2 = data[d].trim()
	                        	if (!str2.isEmpty()) {
		                        	z = Double.parseDouble(data[d].trim())
		                        	output.setValue(row, col, z)
		                        	col++
		                        	if (col == cols) {
		                        		col = 0
		                        		row++
		                        	}
	                        	}
	                        }
						}

                        progress = (int)(100f * row / (rows - 1))
						if (progress > oldProgress) {
							pluginHost.updateProgress(progress)
							oldProgress = progress
							// check to see if the user has requested a cancellation
							if (pluginHost.isRequestForOperationCancelSet()) {
								pluginHost.showFeedback("Operation cancelled")
								return
							}
						}
					}
				}

				output.addMetadataEntry("Created by the "
	                    + descriptiveName + " tool.")
		        output.addMetadataEntry("Created on " + new Date())
				output.close()
		
				// display the output image
				pluginHost.returnData(outputFile)
				
            	// check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
                
            }

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
}

if (args == null) {
    pluginHost.showFeedback("Plugin arguments not set.")
} else {
    def f = new ImportArcAsciiGrid(pluginHost, args, name, descriptiveName)
}
