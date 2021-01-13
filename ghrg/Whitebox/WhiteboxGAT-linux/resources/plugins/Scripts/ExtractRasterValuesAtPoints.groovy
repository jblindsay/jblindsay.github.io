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
import java.awt.Frame
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
import whitebox.interfaces.MapLayer
import whitebox.interfaces.MapLayer.MapLayerType
import whitebox.geospatialfiles.VectorLayerInfo
import whitebox.utilities.Topology
import whitebox.structures.BoundingBox
import whiteboxgis.user_interfaces.AttributesFileViewer
import whiteboxgis.WhiteboxGui
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "ExtractRasterValuesAtPoints"
def descriptiveName = "Extract Raster Values At Points"
def description = "Extracts the values of raster(s) at vector point locations"
def toolboxes = ["GISTools", "StatisticalTools"]

public class ExtractRasterValuesAtPoints implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
	
    public ExtractRasterValuesAtPoints(WhiteboxPluginHost pluginHost, 
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
        	sd.addDialogMultiFile("Select the input raster files", "Input Raster Files:", "Raster Files (*.dep), DEP")
			sd.addDialogFile("Input vector points file.", "Input Vector Points File:", "open", "Vector Files (*.shp), SHP", true, false)
			
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
        	int progress, oldProgress;
            double nodata, x, y, value;
	        int cols, rows, row, col;
	        Object[] rowData;
	        double[][] pointData;
	        double[] data;
	        whitebox.geospatialfiles.shapefile.Point wbGeometry;
        	
            if (args.length < 2) {
                pluginHost.showFeedback("Incorrect number of arguments given to tool.");
                return;
            }
            
            // read the input parameters
			String inputFileString = args[0];
			String[] inputFiles = inputFileString.split(";");
			int numRastersInList = inputFiles.length;
			WhiteboxRasterBase[] rastersInList = new WhiteboxRasterBase[numRastersInList];
			double[] nodataValues = new double[numRastersInList];
			for (int i = 0; i < numRastersInList; i++) {
				rastersInList[i] = new WhiteboxRaster(inputFiles[i], "r");
			    nodataValues[i] = rastersInList[i].getNoDataValue();
			}
            
            String pointsFile = args[1];
			ShapeFile points = new ShapeFile(pointsFile);

			// make sure the vector file is of POINT ShapeType
			if (points.getShapeType() != ShapeType.POINT) {
				pluginHost.showFeedback("The input vector must be of a POINT ShapeType");
				return;
			}
			
			// add a new field to the points attribute table for each raster input
			AttributeTable table = points.getAttributeTable();
			int[] attributeIndex = new int[numRastersInList];
			for (int i = 0; i < numRastersInList; i++) {
				 DBFField field = new DBFField();

            	field.setName("VALUE${i + 1}");
            	field.setDataType(DBFField.DBFDataType.NUMERIC);
            	field.setFieldLength(12);
            	field.setDecimalCount(4);
            	table.addField(field);
            	attributeIndex[i] = table.getFieldCount() - 1;
			}

			// read the data
			int numPoints = points.getNumberOfRecords();
			double[][] outData = new double[numPoints][numRastersInList];
			for (int i = 0; i < numRastersInList; i++) {
				for (ShapeFileRecord record : points.records) {
					int r = record.getRecordNumber() - 1;
	            	pointData = record.getGeometry().getPoints();
	            	x = pointData[0][0];
	            	y = pointData[0][1];
	            	row = rastersInList[i].getRowFromYCoordinate(y);
	            	data = rastersInList[i].getRowValues(row);
					col = rastersInList[i].getColumnFromXCoordinate(x);
					value = data[col];
					outData[r][i] = value;

					progress = (int) (100f * (r + 1) / numPoints);
	                if (progress != oldProgress) {
	                    oldProgress = progress;
	                    pluginHost.updateProgress("Loop ${i + 1} of $numRastersInList:", progress);
	                    if (pluginHost.isRequestForOperationCancelSet()) {
	                        pluginHost.showFeedback("Operation cancelled")
							return
	                    }
	                }
				}
			}

			// output the data
			for (ShapeFileRecord record : points.records) {
            	int r = record.getRecordNumber() - 1;
            	pointData = record.getGeometry().getPoints();
            	x = pointData[0][0];
            	y = pointData[0][1];
                rowData = table.getRecord(r)
                for (int i = 0; i < numRastersInList; i++) {
					value = outData[r][i];
					if (value != nodataValues[i]) {
						rowData[attributeIndex[i]] = value;
					}
                }
                table.updateRecord(r, rowData);

                progress = (int) (100f * (r + 1) / numPoints);
                if (progress != oldProgress) {
                    oldProgress = progress;
                    pluginHost.updateProgress("Saving the data:", progress);
                    if (pluginHost.isRequestForOperationCancelSet()) {
                        pluginHost.showFeedback("Operation cancelled")
						return
                    }
                }
            }
			
			
//            for (ShapeFileRecord record : points.records) {
//            	int r = record.getRecordNumber() - 1;
//            	pointData = record.getGeometry().getPoints();
//            	x = pointData[0][0];
//            	y = pointData[0][1];
//                rowData = table.getRecord(r)
//                for (int i = 0; i < numRastersInList; i++) {
//					col = rastersInList[i].getColumnFromXCoordinate(x);
//					row = rastersInList[i].getRowFromYCoordinate(y);
//					value = rastersInList[i].getValue(row, col);
//					if (value != nodataValues[i]) {
//						rowData[attributeIndex[i]] = value;
//					}
//                }
//                table.updateRecord(r, rowData);
//
//                progress = (int) (100f * (r + 1) / numPoints);
//                if (progress != oldProgress) {
//                    oldProgress = progress;
//                    pluginHost.updateProgress("Saving the data:", progress);
//                    if (pluginHost.isRequestForOperationCancelSet()) {
//                        pluginHost.showFeedback("Operation cancelled")
//						return
//                    }
//                }
//            }

            points.write();
			
			//pluginHost.returnData(pointsFile);
			if (pluginHost.showFeedback("The raster values have been output to the points file's attribute table.\nWould you like to display the table?", 0, 1) == 0) {
				//displayTable(pointsFile)
				MapLayer[] a = pluginHost.getAllMapLayers();
				
				boolean found = false;
				for (MapLayer k : a) {
					if (k.getLayerType() == MapLayerType.VECTOR) {
						VectorLayerInfo vli = (VectorLayerInfo)k
						if (vli.getFileName().equals(pointsFile)) {
							found = true;
							break;
						}
					}
				}
				
				if (!found) {
					pluginHost.returnData(pointsFile);
				}
				
				for (MapLayer k : a) {
					if (k.getLayerType() == MapLayerType.VECTOR) {
						VectorLayerInfo vli = (VectorLayerInfo)k
						if (vli.getFileName().equals(pointsFile)) {
							AttributesFileViewer afv = new AttributesFileViewer((Frame)((WhiteboxGui)(pluginHost)), false, vli);
							afv.setActiveTab(0);
							afv.setSize((int)(500 * 1.61803399), 500);
							afv.setVisible(true);
							break;
						}
					}
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
    def f = new ExtractRasterValuesAtPoints(pluginHost, args, name, descriptiveName)
}
