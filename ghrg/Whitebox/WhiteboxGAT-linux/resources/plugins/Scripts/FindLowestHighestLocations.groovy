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
import java.util.Collections
import java.beans.PropertyChangeEvent
import java.beans.PropertyChangeListener
import java.awt.Color;
import javax.swing.JPanel;
import java.awt.Dimension
import java.nio.file.Paths
import java.nio.file.Files
import com.vividsolutions.jts.geom.*
import whitebox.plugins.PluginInfo
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.*
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.geospatialfiles.shapefile.ShapeFileRecord
import whitebox.utilities.Topology
import whitebox.structures.BoundingBox
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "FindLowestHighestLocations"
def descriptiveName = "Find Lowest/Highest Locations"
def description = "Locates the lowest and/or highest cells in a raster"
def toolboxes = ["GISTools"]

public class FindLowestHighestLocations implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
	
    public FindLowestHighestLocations(WhiteboxPluginHost pluginHost, 
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
        	DialogFile dfIn1 = sd.addDialogFile("Input raster file", "Input Raster File:", "open", "Raster Files (*.dep), DEP", true, false)
            DialogFile dfIn2 = sd.addDialogFile("Output vector points file.", "Output Vector Points File:", "save", "Vector Files (*.shp), SHP", true, false)
			sd.addDialogComboBox("Which point to find?", "Which point(s) should be located?", ["lowest", "highest", "both"], 0)
            
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
        	int progress, oldProgress
            double nodata
	        int cols, rows 
	        Object[] rowData
	        whitebox.geospatialfiles.shapefile.Point wbGeometry
        	
            if (args.length < 2) {
                pluginHost.showFeedback("Incorrect number of arguments given to tool.")
                return
            }
            
            // read the input parameters
            String inputFile = args[0];
            WhiteboxRaster input = new WhiteboxRaster(inputFile, "r")
			cols = input.getNumberColumns()
			rows = input.getNumberRows()
			nodata = input.getNoDataValue()
			
            String outputFile = args[1];
			DBFField[] fields = new DBFField[4];

            fields[0] = new DBFField();
            fields[0].setName("FID");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setFieldLength(10);
            fields[0].setDecimalCount(0);

            fields[1] = new DBFField();
            fields[1].setName("XCOORD");
            fields[1].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[1].setFieldLength(14);
            fields[1].setDecimalCount(4);

            fields[2] = new DBFField();
            fields[2].setName("YCOORD");
            fields[2].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[2].setFieldLength(14);
            fields[2].setDecimalCount(4);

			fields[3] = new DBFField();
            fields[3].setName("VALUE");
            fields[3].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[3].setFieldLength(12);
            fields[3].setDecimalCount(4);
            
            ShapeFile output = new ShapeFile(outputFile, ShapeType.POINT, fields);

			String whichPoints = "both";
			if (args.length == 3 && args[2] != null && !(args[2].trim().isEmpty())) {
            	whichPoints = args[2].toLowerCase().trim();
			}
            
            if (!(whichPoints.equals("lowest")) && !(whichPoints.equals("highest"))) {
            	whichPoints = "both";
            }

			double minValue = Double.POSITIVE_INFINITY
			double minX = -1, minY = -1
			double maxValue = Double.NEGATIVE_INFINITY
			double maxX = -1, maxY = -1
			boolean flag = false
            for (int row = 0; row < rows; row++) {
				double[] data = input.getRowValues(row)
				for (int col = 0; col < cols; col++) {
					if (data[col] != nodata) {
						if (data[col] < minValue) {
							minValue = data[col];
							minX = input.getXCoordinateFromColumn(col);
							minY = input.getYCoordinateFromRow(row);
						}
						if (data[col] > maxValue) {
							maxValue = data[col];
							maxX = input.getXCoordinateFromColumn(col);
							maxY = input.getYCoordinateFromRow(row);
						}
					}
				}
				progress = (int)(100f * row / rows)
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

			int FID = 1
			if (whichPoints.equals("lowest")) {
				// output the point
				wbGeometry = new whitebox.geospatialfiles.shapefile.Point(minX, minY);                  
                rowData = new Object[4]
                rowData[0] = new Double(FID)
                rowData[1] = minX
                rowData[2] = minY
                rowData[3] = minValue
                output.addRecord(wbGeometry, rowData);
                FID++
			}
			if (whichPoints.equals("highest")) {
				// output the point
				wbGeometry = new whitebox.geospatialfiles.shapefile.Point(maxX, maxY);                  
                rowData = new Object[4]
                rowData[0] = new Double(FID)
                rowData[1] = maxX
                rowData[2] = maxY
                rowData[3] = maxValue
                output.addRecord(wbGeometry, rowData);
			}

			if (whichPoints.equals("both")) {
				// output the points
				wbGeometry = new whitebox.geospatialfiles.shapefile.Point(minX, minY);                  
                rowData = new Object[4]
                rowData[0] = new Double(FID)
                rowData[1] = minX
                rowData[2] = minY
                rowData[3] = minValue
                output.addRecord(wbGeometry, rowData);
                FID++
                
				wbGeometry = new whitebox.geospatialfiles.shapefile.Point(maxX, maxY);                  
                rowData = new Object[4]
                rowData[0] = new Double(FID)
                rowData[1] = maxX
                rowData[2] = maxY
                rowData[3] = maxValue
                output.addRecord(wbGeometry, rowData);
			}

            output.write()
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
}

if (args == null) {
    pluginHost.showFeedback("Plugin arguments not set.")
} else {
    def f = new FindLowestHighestLocations(pluginHost, args, name, descriptiveName)
}
