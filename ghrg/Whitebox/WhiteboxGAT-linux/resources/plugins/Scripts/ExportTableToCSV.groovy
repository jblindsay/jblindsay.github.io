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
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.ui.plugin_dialog.*
import whitebox.utilities.StringUtilities
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "ExportTableToCSV"
def descriptiveName = "Export Table to CSV"
def description = "Exports an attribute table to a CSV text file."
def toolboxes = ["DatabaseTools"]

public class ExportTableToCSV implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public ExportTableToCSV(WhiteboxPluginHost pluginHost, 
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
			sd.addDialogFile("Input file", "Input Vector File:", "open", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogFile("Output comma-separated values (CSV) file", "Output Comma-Separated Values (CSV) File:", "saveAs", "CSV Files (*.csv), CSV", true, false)
            sd.addDialogCheckBox("Export field names as header", "Export field names as header", true)
			
			// resize the dialog to the standard size and display it
			sd.setSize(800, 400)
			sd.visible = true
		}
	}

	private void execute(String[] args) {
		try {
	  		int progress, oldProgress, line
	  		String delimiter = ","
	  		
			if (args.length != 3) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			
			// read the input parameters
			String inputFile = args[0]
			String databaseFile = inputFile.trim().replace(".shp", ".dbf")
			String csvFile = args[1]
			boolean exportFieldNames = Boolean.parseBoolean(args[2])

			if (!(new File(databaseFile)).exists()) {
				pluginHost.showFeedback("The attribute file could not be located")
				return
			}
			
			File outFile = new File(csvFile)
			if (outFile.exists()) {
				outFile.delete()
			}
			AttributeTable attributes = new AttributeTable(databaseFile)
            int numRecords = attributes.getNumberOfRecords()
            int fieldCount = attributes.getFieldCount()
            outFile.withWriter { out ->
            	if (exportFieldNames) {
            		StringBuilder sb = new StringBuilder()
            		(attributes.getAttributeTableFieldNames()).each() { str ->
            			sb.append(str + ",")
            		}
            		out.writeLine(sb.toString())
            	}
		        for (int r = 0; r < numRecords; r++) {
		        	Object[] rec = attributes.getRecord(r)
		        	StringBuilder sb = new StringBuilder()
		        	for (int c = 0; c < fieldCount; c++) {
						sb.append(rec[c].toString().replace(",", ";")).append(",")
		        	}
		        	out.writeLine(sb.toString())
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
	def tdf = new ExportTableToCSV(pluginHost, args, name, descriptiveName)
}
