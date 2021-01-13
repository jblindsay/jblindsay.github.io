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
import java.awt.Dimension
import java.awt.Color
import java.util.Date
import java.util.ArrayList
import java.util.PriorityQueue
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.WhiteboxRasterBase.DataScale
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "SimpleRegionGrow"
def descriptiveName = "Simple Region Grow"
def description = "Performs a simple region-growing operation."
def toolboxes = ["ImageProc"]

public class SimpleRegionGrow {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public SimpleRegionGrow(WhiteboxPluginHost pluginHost, 
		String[] args, String name, String descriptiveName) {
		this.pluginHost = pluginHost
		this.descriptiveName = descriptiveName
			
		if (args.length > 0) {
			execute(args)
		} else {
			// create an ActionListener to handle the return from the dialog
			def ac = new ActionListener() {
	            public void actionPerformed(ActionEvent event) {
	                if (event.getActionCommand().equals("ok")) {
			    		args = sd.collectParameters()
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
	        };
			// Create a dialog for this tool to collect user-specified
			// tool parameters.
		 	sd = new ScriptDialog(pluginHost, descriptiveName, ac)	
		
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
			sd.addDialogMultiFile("Select the input raster files", "Input Raster Files:", "Raster Files (*.dep), DEP")
			//sd.addDialogFile("Input digital elevation model (DEM) file", "Input Digital Elevation Model (DEM):", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Input seed points file", "Input Seed Points File:", "open", "Shapefiles (*.shp), SHP", true, false)
            def btn = sd.addDialogButton("Create a new seed point vector now...", "centre")
			btn.addActionListener(new ActionListener() {
 	            public void actionPerformed(ActionEvent e) {
 	            	// find an appropriate file name for it
					def outputFile = pluginHost.getWorkingDirectory() + "Points.shp";
					def file = new File(outputFile);
					if (file.exists()) {
						for (int i = 1; i < 101; i++) {
							outputFile = pluginHost.getWorkingDirectory() + "Points${i}.shp";
							file = new File(outputFile);
							if (!file.exists()) {
								break;
							}
						}
					}
					DBFField[] fields = new DBFField[1];
		            
		            fields[0] = new DBFField();
		            fields[0].setName("FID");
		            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
		            fields[0].setFieldLength(10);
		            fields[0].setDecimalCount(0);
		            
		            ShapeFile output = new ShapeFile(outputFile, ShapeType.POINT, fields);
					output.write();
		            
		            pluginHost.returnData(outputFile);
		            
		            pluginHost.editVector();

		            pluginHost.showFeedback("Press the 'Digitize New Feature' icon on the toolbar \n" +
		            "to add a point. Then toggle the 'Edit Vector' icon when you \n" + 
		            "are done digitizing");
 	            }
 	        });  
			//sd.addDialogFile("Output text file (optional)", "Output Text File (blank for none):", "save", "Text Files (*.csv), CSV", true, true)
            sd.addDialogFile("Output File", "Output File:", "save", "Raster Files (*.dep), DEP", true, false)
            sd.addDialogDataInput("Similarity Threshold", "Similarity Threshold:", "", true, false)
            sd.addDialogComboBox("Neighbourhood type", "Neighbourhood type:", ["von Neumann (4 cell)", "Moore (8 cell)"], 1)
	
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
	  		int progress, oldProgress, dir
	  		int[] dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ]
			int[] dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ]
			double x, y, z, zN, diff
			int rows, cols, rowsLessOne
        	int row, col, rowN, colN, i;
        	boolean flag
        	int numPoints = 1;
        	GridCell gc
	        
			if (args.length < 4) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String inputFileString = args[0];
			String[] inputFiles = inputFileString.split(";");
			int numRastersInList = inputFiles.length;
			WhiteboxRaster[] rastersInList = new WhiteboxRaster[numRastersInList];
			double[] nodataValues = new double[numRastersInList];
			for (i = 0; i < numRastersInList; i++) {
				rastersInList[i] = new WhiteboxRaster(inputFiles[i], "r");
			    nodataValues[i] = rastersInList[i].getNoDataValue();
			    if (i == 0) {
			    	rows = rastersInList[i].getNumberRows()
			    	rowsLessOne = rows -1
					cols = rastersInList[i].getNumberColumns()
			    }
			}
			String inputSeedFile = args[1]
			String outputFile = args[2]
			double allowableDiff = 0.0;
			if (!args[3].trim().isEmpty() && !(args[3].toLowerCase().equals("not specified"))) {
				allowableDiff = Double.parseDouble(args[3])
				allowableDiff = allowableDiff * allowableDiff
			}
			if (args[4].contains("4") || args[4].toLowerCase().contains("von")) {
				dX = [ 1, 0, -1, 0 ]
				dY = [ 0, 1, 0, -1 ]
			}

           	int dxSize = dX.size()
           	
            def input = new ShapeFile(inputSeedFile)
            ShapeType shapeType = input.getShapeType()
			if (shapeType.getBaseType() != ShapeType.POINT) {
				pluginHost.showFeedback("Error: The seed file must be of a POINT shape type.")
				return;
			}

			double backgroundVal = -999999
			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", 
  		  	     inputFiles[0], DataType.FLOAT, backgroundVal)
  		  	output.setDataScale(DataScale.CATEGORICAL);
            output.setPreferredPalette("qual.pal");
  		  	output.setNoDataValue(nodataValues[0])
			
			int featureNum = 1;
			numPoints = input.getNumberOfRecords();
        	double[][] point

			double[] seedVals = new double[numRastersInList]
        	double[] vals = new double[numRastersInList]
			PriorityQueue<GridCell> queue = new PriorityQueue<GridCell>((2 * rows + 2 * cols) * 2);
			int numSolvedCells = 0
			int numCells = rows * cols
			progress = 1
			for (ShapeFileRecord record : input.records) {
        	 	point = record.getGeometry().getPoints()
				col = rastersInList[0].getColumnFromXCoordinate(point[0][0])
				row = rastersInList[0].getRowFromYCoordinate(point[0][1])
				//z = rastersInList[0].getValue(row, col);
				boolean isValid = true
				for (i = 0; i < numRastersInList; i++) {
					seedVals[i] = rastersInList[i].getValue(row, col);
					if (seedVals[i] == nodataValues[i]) { isValid = false }
				}
				
				if (isValid) {
					numSolvedCells = 0
					output.setValue(row, col, featureNum)
					queue.add(new GridCell(row, col, 0d))
					while (!queue.isEmpty()) {
						gc = queue.poll();
		                row = gc.row;
		                col = gc.col;
		                for (i = 0; i < dxSize; i++) {
		                	rowN = row + dY[i];
		                    colN = col + dX[i];
		                    if (output.getValue(rowN, colN) == backgroundVal) {
		                        // is it similar enough to the seed cells
		                        diff = 0
		                        isValid = true
		                        for (int j = 0; j < numRastersInList; j++) {
									vals[j] = rastersInList[j].getValue(rowN, colN);
									if (vals[j] == nodataValues[j]) { 
										isValid = false 
										break
									}
									diff += (seedVals[j] - vals[j]) * (seedVals[j] - vals[j])
								}
								if (isValid && (diff <= allowableDiff)) {
									
		                        	output.setValue(rowN, colN, featureNum)
		                        	gc = new GridCell(rowN, colN, diff);
		                        	queue.add(gc);
		                        	numSolvedCells++;
		                        	pluginHost.updateProgress("Seed $featureNum ($numSolvedCells cells found)", progress)
								} else {
									output.setValue(rowN, colN, nodataValues[0])
								}
		                    }
		                }
					}
				} else {
					output.setValue(row, col, nodataValues[0])
				}
				
				progress = (int)(100f * featureNum / numPoints)
        		if (progress != oldProgress) {
					pluginHost.updateProgress("Seed $featureNum ($numSolvedCells cells found)", progress)
        			oldProgress = progress
        			// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
        		}
        		
				featureNum++
        	}
        	
        	for (i = 0; i < numRastersInList; i++) {
				rastersInList[i].close()
        	}

        	oldProgress = -1
  		  	for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					z = output.getValue(row, col)
					if (z == backgroundVal) {
						output.setValue(row, col, nodataValues[0])
					}
				}
  		  		progress = (int)(100f * row / rowsLessOne)
				if (progress != oldProgress) {
					pluginHost.updateProgress("Saving Data:", progress)
					oldProgress = progress

					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}
        	output.close()

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

    @CompileStatic
    class GridCell implements Comparable<GridCell> {

        public int row;
        public int col;
        public double z;

        public GridCell(int Row, int Col, double Z) {
            row = Row;
            col = Col;
            z = Z;
        }

        @Override
        public int compareTo(GridCell other) {
            final int BEFORE = -1;
            final int EQUAL = 0;
            final int AFTER = 1;

            if (this.z < other.z) {
                return BEFORE;
            } else if (this.z > other.z) {
                return AFTER;
            }

            if (this.row < other.row) {
                return BEFORE;
            } else if (this.row > other.row) {
                return AFTER;
            }

            if (this.col < other.col) {
                return BEFORE;
            } else if (this.col > other.col) {
                return AFTER;
            }

            return EQUAL;
        }
    }
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def tdf = new SimpleRegionGrow(pluginHost, args, name, descriptiveName)
}
