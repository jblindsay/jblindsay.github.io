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
import java.util.HashSet;
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
def name = "MergeShapefiles"
def descriptiveName = "Merge Shapefiles"
def description = "Merges two shapefiles of the same shape type"
def toolboxes = ["VectorTools"]

public class MergeShapefiles implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public MergeShapefiles(WhiteboxPluginHost pluginHost, 
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
			//sd.addDialogFile("Input file", "Input Vector File:", "open", "Vector Files (*.shp), SHP", true, false)
            //sd.addDialogFile("Input file", "Input Vector File:", "open", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogMultiFile("Select input vector files", "Input Vector Files:", "Vector Files (*.shp), SHP")
			sd.addDialogFile("Output file", "Output Vector File:", "save", "Vector Files (*.shp), SHP", true, false)
            
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
	  		if (args.length < 2) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}

			String inputFileString = args[0]
			String[] inputFiles = inputFileString.split(";")
			if (inputFiles.length < 2) {
				pluginHost.showFeedback("At least two input vector files are needed to run this tool.")
				return
			}
			String outputFile = args[1]

			ShapeFile[] inputs = new ShapeFile[inputFiles.length]
			ShapeType shapetype;
			ArrayList<DBFField> fields = new ArrayList<DBFField>()
			int[][] fieldsMap = new int[inputFiles.length][]
			for (int f = 0; f < inputFiles.length; f++) {
				String inputFile = inputFiles[f].trim()
				if (!inputFile.isEmpty()) {
					inputs[f] = new ShapeFile(inputFile)
					if (f == 0) {
						shapetype = inputs[f].getShapeType()
					} else {
						if (inputs[f].getShapeType() != shapetype) {
							pluginHost.showFeedback("ERROR: All of the input vectors must be of the same ShapeType.")
							return
						}
					}
					AttributeTable table = new AttributeTable(inputFile.replace(".shp", ".dbf"))
					DBFField[] myFields = table.getAllFields()
					fieldsMap[f] = new int[myFields.length]
					for (i = 0; i < myFields.length; i++) {
						DBFField field = myFields[i]
						if (!fields.contains(field)) {
							// add it
							fields.add(field)
							fieldsMap[f][i] = fields.size() - 1;
						} else {
							// find it
							fieldsMap[f][i] = fields.indexOf(field);
						}
					}
				}
			}

			int numOutFields = fields.size();
			DBFField[] outFields = new DBFField[numOutFields]
			for (i = 0; i < fields.size(); i++) {
				outFields[i] = fields.get(i)
			}
			ShapeFile output = new ShapeFile(outputFile, shapetype, outFields);

			Object[] rowData1;
			Object[] rowData2;
			for (int f = 0; f < inputFiles.length; f++) {
				AttributeTable table = inputs[f].getAttributeTable()
				int r = 0;
				int numFeatures = inputs[f].getNumberOfRecords();
				for (ShapeFileRecord record : inputs[f].records) {
					rowData1 = table.getRecord(record.getRecordNumber() - 1); // table records are base zero
					rowData2 = new Object[numOutFields]
					for (i = 0; i < rowData1.length; i++) {
						rowData2[fieldsMap[f][i]] = rowData1[i]
					}
					output.addRecord(record.getGeometry(), rowData2);

					r++;
	                progress = (int)(100f * r / numFeatures);
	            	if (progress != oldProgress) {
						pluginHost.updateProgress("Merging file ${f + 1} of ${inputFiles.length}:", progress)
	            		oldProgress = progress
	            		// check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
	            	}
				}
			}

			output.write();

			pluginHost.returnData(outputFile);
			
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
	def tdf = new MergeShapefiles(pluginHost, args, name, descriptiveName)
}
