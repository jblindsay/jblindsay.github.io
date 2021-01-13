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
 
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Point
import com.vividsolutions.jts.geom.Geometry
import com.vividsolutions.jts.geom.GeometryCollection
import com.vividsolutions.jts.geom.GeometryFactory
import java.awt.event.ActionListener
import java.awt.event.ActionEvent
import java.io.File
import java.util.Date
import java.util.ArrayList
import java.util.Arrays
import java.beans.PropertyChangeEvent
import java.beans.PropertyChangeListener
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.ui.plugin_dialog.*
import whitebox.utilities.FileUtilities;
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "CentroidVector"
def descriptiveName = "Centroid (Vector)"
def description = "Identifies the centroid point of a group of input vector features."
def toolboxes = ["VectorTools"]

public class CentroidVector implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
	
    public CentroidVector(WhiteboxPluginHost pluginHost, 
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
            DialogFile df = sd.addDialogFile("Input polygon vector file", "Input Polygon Vector File:", "open", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogFile("Output points file", "Output Points Vector File:", "save", "Vector Files (*.shp), SHP", true, false)
            DialogCheckBox cb = sd.addDialogCheckBox("Create one point for each part in a multipart feature?", "Create multiple points for multipart features?", false)
			cb.setVisible(false)
			def lstrDF = { evt -> if (evt.getPropertyName().equals("value")) { 
            		String fileName = df.getValue()
            		if (fileName != null && !fileName.isEmpty()) { 
            			def file = new File(fileName)
            			if (file.exists()) {
            				ShapeFile sf = new ShapeFile(fileName)
            				if (sf.getShapeType().getBaseType() == ShapeType.POLYGON) {
            					cb.setVisible(true)
            				} else {
            					cb.setVisible(false)
            				}
            			}
            		}
            	}
            } as PropertyChangeListener
            df.addPropertyChangeListener(lstrDF)

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
            int i, n, progress, oldProgress
            
            if (args.length < 2) {
                pluginHost.showFeedback("Incorrect number of arguments given to tool.")
                return
            }
            
            // read the input parameters
            String inputFile = args[0]
            String outputFile = args[1]
            boolean breakApartMultipart = false
            if (args.length >= 3) {
            	breakApartMultipart = Boolean.parseBoolean(args[2])
            }
            
           	def input = new ShapeFile(inputFile);
			def shapeType = input.getShapeType()
			
			def numRecs = input.getNumberOfRecords();
			AttributeTable table = input.getAttributeTable()
			DBFField[] fields = table.getAllFields()
			
			def output = new ShapeFile(outputFile, ShapeType.POINT, fields);
			
			com.vividsolutions.jts.geom.Point p = null;
			com.vividsolutions.jts.geom.Geometry[] jtsGeometries = null;
			GeometryCollection jtsCollection = null;
			GeometryFactory factory = new GeometryFactory()
			Coordinate pCoord
			oldProgress = -1
			progress = 0;
			if (shapeType.getBaseType() != ShapeType.POINT) {
				for (i = 0; i < numRecs; i++) {
					ShapeFileRecord record = input.getRecord(i)
				    if (record.getShapeType() != ShapeType.NULLSHAPE) {
				        jtsGeometries = record.getGeometry().getJTSGeometries();
				        if (breakApartMultipart) {
				            for (n = 0; n < jtsGeometries.length; n++) {
					            p = jtsGeometries[n].getCentroid();
					            pCoord = p.getCoordinate();
					            whitebox.geospatialfiles.shapefile.Point wbGeometry = new whitebox.geospatialfiles.shapefile.Point(pCoord.x, pCoord.y);
					    		Object[] rowData = table.getRecord(i);
					    		output.addRecord(wbGeometry, rowData);
					        }
				        } else {
				        	jtsCollection = new GeometryCollection(jtsGeometries, factory)
					        p = jtsCollection.getCentroid()
					        pCoord = p.getCoordinate();
					        whitebox.geospatialfiles.shapefile.Point wbGeometry = new whitebox.geospatialfiles.shapefile.Point(pCoord.x, pCoord.y);
						    Object[] rowData = table.getRecord(i);
						    output.addRecord(wbGeometry, rowData);
				        }
				    }
				    progress = (int)(100f * i / (numRecs - 1))
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
			} else {
				com.vividsolutions.jts.geom.Geometry[] geomList = new com.vividsolutions.jts.geom.Geometry[numRecs]
				for (i = 0; i < numRecs; i++) {
					ShapeFileRecord record = input.getRecord(i)
				    if (record.getShapeType() != ShapeType.NULLSHAPE) {
				        jtsGeometries = record.getGeometry().getJTSGeometries();
				        geomList[i] = jtsGeometries[0]
				    }
				    progress = (int)(100f * i / (numRecs - 1))
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
				jtsCollection = new GeometryCollection(geomList, factory)
			    p = jtsCollection.getCentroid()
		        pCoord = p.getCoordinate();
		        whitebox.geospatialfiles.shapefile.Point wbGeometry = new whitebox.geospatialfiles.shapefile.Point(pCoord.x, pCoord.y);
			    fields = new DBFField[1];

	            fields[0] = new DBFField();
	            fields[0].setName("FID");
	            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
	            fields[0].setFieldLength(10);
	            fields[0].setDecimalCount(0);
			    output = new ShapeFile(outputFile, ShapeType.POINT, fields);
				
				Object[] rowData = new Object[1]
			    rowData[0] = new Double(1)
			    output.addRecord(wbGeometry, rowData);   
			}
			
			output.write()
			
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
}

if (args == null) {
    pluginHost.showFeedback("Plugin arguments not set.")
} else {
    def f = new CentroidVector(pluginHost, args, name, descriptiveName)
}
