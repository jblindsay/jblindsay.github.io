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
import java.util.Collections
import java.beans.PropertyChangeEvent
import java.beans.PropertyChangeListener
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.geospatialfiles.shapefile.attributes.AttributeTable
import whitebox.geospatialfiles.shapefile.attributes.DBFField
import whitebox.geospatialfiles.shapefile.attributes.DBFField.DBFDataType
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.ui.plugin_dialog.*
import whitebox.utilities.StringUtilities
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "AddPointCoordinatesToTable"
def descriptiveName = "Add Point Coordinates To Table"
def description = "Adds the xy coordinates of a point file to its attribute table"
def toolboxes = ["DatabaseTools"]

public class AddPointCoordinatesToTable implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public AddPointCoordinatesToTable(WhiteboxPluginHost pluginHost, 
		String[] args, String name, String descriptiveName) {
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
			sd.addDialogFile("Input file", "Input Vector Points File:", "open", "Vector Files (*.shp), SHP", true, false)
            
			// resize the dialog to the standard size and display it
			sd.setSize(800, 400)
			sd.visible = true
		}
	}

	@CompileStatic
	private void execute(String[] args) {
		try {
	  		int progress, oldProgress, i
	  		double x, y
	  		Object[] rec;
	  		int numFeatures
			if (args.length < 1) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			
			// read the input parameters
			String inputFile = args[0]

			ShapeFile input = new ShapeFile(inputFile)
			if (input.getShapeType() != ShapeType.POINT) {
            	pluginHost.showFeedback("The input file must be of a POINT ShapeType.")
            	return
            }

			AttributeTable table = input.getAttributeTable()

			DBFField field = new DBFField();
			field.setName("XCOORD");
            field.setDataType(DBFField.DBFDataType.NUMERIC);
            field.setFieldLength(14);
            field.setDecimalCount(5);
            table.addField(field);

            field = new DBFField();
            field.setName("YCOORD");
            field.setDataType(DBFField.DBFDataType.NUMERIC);
            field.setFieldLength(14);
            field.setDecimalCount(5);
            table.addField(field);
			
			DBFField[] fields = table.getAllFields()

			numFeatures = input.getNumberOfRecords()
			double[][] point
			int numFields = table.getFieldCount();
			i = 0
			for (ShapeFileRecord record : input.records) {
				point = record.getGeometry().getPoints()
				x = point[0][0];
				y = point[0][1];

				Object[] rowData = table.getRecord(i);
				rowData[numFields - 2] = new Double(x);
				rowData[numFields - 1] = new Double(y);
				table.updateRecord(i, rowData);
				
				i++;
                progress = (int)(100f * i / numFeatures);
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Reading Points:", progress)
            		oldProgress = progress
            		// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}
			}

			pluginHost.showFeedback("Operation complete.")

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
	def tdf = new AddPointCoordinatesToTable(pluginHost, args, name, descriptiveName)
}
