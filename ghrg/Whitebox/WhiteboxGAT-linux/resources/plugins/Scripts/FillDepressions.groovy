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
import java.util.PriorityQueue
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import groovy.transform.CompileStatic
import groovy.time.*

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
//def name = "FillDepressions"
//def descriptiveName = "Fill Depressions"
//def description = "This tool fills all of the depressions in a DEM using the Wang and Liu (2006) algorithm."
//def toolboxes = ["DEMPreprocessing"]

public class FillDepressions implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public FillDepressions(WhiteboxPluginHost pluginHost, 
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
			sd.addDialogFile("Input digital elevation model (DEM) file", "Input Digital Elevation Model (DEM) File:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Output file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogDataInput("Flat increment value", "Flat Increment Value:", "0.001", true, false)
            
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
	  	
		if (args.length != 3) {
			pluginHost.showFeedback("Incorrect number of arguments given to tool.")
			return
		}
		
		def timeStart = new Date()
		
		int progress
		int oldProgress = -1
		
		// read the input parameters
		String inputHeader = args[0]
		String outputHeader = args[1]
		double SMALL_NUM = Double.parseDouble(args[2])
		
        pluginHost.updateProgress("Initializing: ", -1);
        int row_n, col_n;
        int row, col;
        double z_n;
        long k = 0;
        GridCell gc = null;
        double z;
        int[] Dy = [-1, 0, 1, 1, 1, 0, -1, -1];
        int[] Dx = [1, 1, 1, 0, -1, -1, -1, 0];
//        boolean flag = false;
        
        WhiteboxRaster image = new WhiteboxRaster(inputHeader, "r");
        int rows = image.getNumberRows();
        int rowsLessOne = rows - 1
        int cols = image.getNumberColumns();
        int numCells = 0;
        String preferredPalette = image.getPreferredPalette();
        
        double noData = image.getNoDataValue();

        double[][] output = new double[rows][cols];
        double[][] input = new double[rows + 2][cols + 2];
        for (row = 0; row < rows + 2; row++) {
            input[row][0] = noData;
            input[row][cols + 1] = noData;
        }
        
        for (col = 0; col < cols + 2; col++) {
            input[0][col] = noData;
            input[rows + 1][col] = noData;
        }
        
        double[] data;
        for (row = 0; row < rows; row++) {
            data = image.getRowValues(row);
            for (col = 0; col < cols; col++) {
                output[row][col] = -999;
                input[row + 1][col + 1] = data[col];
            }
        }
        image.close();

        //TimeDuration duration1 = TimeCategory.minus(new Date(), timeStart)
		
        
        // initialize and fill the priority queue.
        pluginHost.updateProgress("Loop 1: ", -1);

        PriorityQueue<GridCell> queue = new PriorityQueue<GridCell>((2 * rows + 2 * cols) * 2);
		oldProgress = 0
        for (row = 0; row < rows; row++) {
            for (col = 0; col < cols; col++) {
                z = input[row + 1][col + 1];
                if (z != noData) {
                    numCells++;
                    for (int i = 0; i < 8; i++) {
                        row_n = row + Dy[i];
                        col_n = col + Dx[i];
                        z_n = input[row_n + 1][col_n + 1];
                        if (z_n == noData) {
                            // it's an edge cell.
                            queue.add(new GridCell(row, col, z));
                        	output[row][col] = z;
                            break
                        }
                    }
                } else {
//                    k++;
                    output[row][col] = noData;
                }

            }
            progress = (int)(100f * row / rowsLessOne)
			if ((progress - oldProgress) == 2) {
				pluginHost.updateProgress(progress)
				oldProgress = progress
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
			}
        }
		
		//TimeDuration duration2 = TimeCategory.minus(new Date(), timeStart)
		
        // now fill!
        pluginHost.updateProgress("Loop 2: ", 0);

        oldProgress = 0
        while (!queue.isEmpty()) {
            gc = queue.poll();
            row = gc.row;
            col = gc.col;
            z = gc.z;
            for (int i = 0; i < 8; i++) {
                row_n = row + Dy[i];
                col_n = col + Dx[i];
                z_n = input[row_n + 1][col_n + 1];
                if ((output[row_n][col_n] == -999) && (z_n != noData)) {
                    if (z_n <= z) {
                        z_n = z + SMALL_NUM;
                    }
                    output[row_n][col_n] = z_n;
                    queue.add(new GridCell(row_n, col_n, z_n));
                }
            }
            k++;
            progress = (int)(k * 100f / numCells);
			if ((progress - oldProgress) == 5) {
				pluginHost.updateProgress(progress)
				oldProgress = progress
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
			}
        }

        //TimeDuration duration3 = TimeCategory.minus(new Date(), timeStart)
		

        pluginHost.updateProgress("Saving Data: ", 0);
        WhiteboxRaster outputFile = new WhiteboxRaster(outputHeader, "rw", 
                inputHeader, DataType.DOUBLE, -999);
        outputFile.setPreferredPalette(preferredPalette);
        
        oldProgress = 0
        for (row = 0; row < rows; row++) {
            outputFile.setRowValues(row, output[row]);
            
			progress = (int)(100f * row / rowsLessOne)
			if ((progress - oldProgress) == 2) {
				pluginHost.updateProgress(progress)
				oldProgress = progress
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
			}
        }

        //TimeDuration duration4 = TimeCategory.minus(new Date(), timeStart)
		
            
		outputFile.addMetadataEntry("Created by the "
                    + descriptiveName + " tool.")
        outputFile.addMetadataEntry("Created on " + new Date())
		outputFile.close()

		//TimeDuration duration5 = TimeCategory.minus(new Date(), timeStart)
		//pluginHost.showFeedback("${duration1} ${duration2} ${duration3} ${duration4}")

		TimeDuration duration = TimeCategory.minus(new Date(), timeStart)
		pluginHost.showFeedback("Operation completed in ${duration} seconds.")
	

		// display the output image
		pluginHost.returnData(outputHeader)

		// reset the progress bar
		pluginHost.updateProgress(0)
	  } catch (OutOfMemoryError oe) {
            pluginHost.showFeedback("An out-of-memory error has occurred during operation.")
      } catch (Exception e) {
            pluginHost.showFeedback("An error has occurred during operation. See log file for details.")
            pluginHost.logException("Error in " + descriptiveName, e)
      } finally {
      	pluginHost.updateProgress("Progress:", 0)
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
	class GridCell implements Comparable<GridCell> {

        public int row;
        public int col;
        public double z;

        public GridCell(int row, int col, double z) {
            this.row = row;
            this.col = col;
            this.z = z;
        }

        @Override
        public int compareTo(GridCell cell) {
//        	return (int)((this.z - cell.z) * 10000)
            if (this.z < cell.z) {
                return -1;
            } else if (this.z > cell.z) {
                return 1;
            }

            if (this.row < cell.row) {
                return -1;
            } else if (this.row > cell.row) {
                return 1;
            }

            if (this.col < cell.col) {
                return -1;
            } else if (this.col > cell.col) {
                return 1;
            }

            return 0;
        }
    }
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def f = new FillDepressions(pluginHost, args, name, descriptiveName)
}
