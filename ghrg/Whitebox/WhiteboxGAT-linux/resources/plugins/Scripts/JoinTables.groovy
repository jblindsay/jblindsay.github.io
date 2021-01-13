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
def name = "JoinTables"
def descriptiveName = "Join Tables"
def description = "Joins (merges) two attribute tables."
def toolboxes = ["DatabaseTools"]

public class JoinTables implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public JoinTables(WhiteboxPluginHost pluginHost, 
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
			DialogFieldSelector dfs = sd.addDialogFieldSelector("Select the field that is the 'Primary Key' or unique ID field.", "Select the field that is the 'Primary Key' or ID field.", false)
            DialogFile df = sd.addDialogFile("File from which to derive the appended data", "Input Shapefile:", "open", "Shapefiles (*.shp), SHP", true, false)
            String[] listItems = []
            DialogComboBox foreignKey = sd.addDialogComboBox("Select the field that is the 'Foreign Key' or unique ID field.", "'Foreign Key' or ID field", listItems, 0)
			DialogList includedFields = sd.addDialogList("Select the fields to include in the merge.", "Select the fields to include in the merge.", listItems, true)

            def lstrDF = { evt -> if (evt.getPropertyName().equals("value")) { 
            		if (df.getValue() != null && !(df.getValue().isEmpty())) {
	            		String fileName = (df.getValue()).replace(".shp", ".dbf")
	            		if (fileName != null && !fileName.isEmpty()) { 
	            			def file = new File(fileName)
	            			if (file.exists()) {
		            			AttributeTable table = new AttributeTable(fileName)
		            			String[] fields = table.getAttributeTableFieldNames()
		            			foreignKey.setListItems(fields)
		            			includedFields.setListItems(fields)
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

	@CompileStatic
	private void execute(String[] args) {
		try {
	  		int progress, oldProgress, line
	  		Object[] rec;
	  		int numFeatures
			if (args.length != 4) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			
			// read the input parameters
			String[] inputData = args[0].split(";")
			if (inputData.length < 2) {
				pluginHost.showFeedback("Error reading one of the input parameters.")
				return
			}
			String databaseFile1 = inputData[0].trim().replace(".shp", ".dbf")
			String primaryKeyString = inputData[1].trim()
			String databaseFile2 = args[1].replace(".shp", ".dbf")
			String foreignKeyString = args[2]
			String[] includeFieldString = args[3].split(";")
			
			int numIncludedFields = includeFieldString.length
			DBFField[] includedDBFFields = new DBFField[numIncludedFields]
			
			AttributeTable table2 = new AttributeTable(databaseFile2)
			
			// get the field numbers of the foreign key and include fields
			int foreignKey = table2.getFieldColumnNumberFromName(foreignKeyString)
			int[] includedFields = new int[numIncludedFields]
			String[] fields = table2.getAttributeTableFieldNames()
			for (int a in 0..<fields.length) {
				for (int b in 0..<numIncludedFields) {
					if (fields[a].equals(includeFieldString[b])) {
						includedFields[b] = a
						includedDBFFields[b] = table2.getField(a)
					}
				}
			}
				
//			if (foreignKey == -1) {
//				pluginHost.showFeedback("Could not locate the foreign key (unique ID field).")
//				return
//			}
			if (foreignKey == null || foreignKey < 0) {
        		pluginHost.showFeedback("Could not locate the foreign key (unique ID field). Check your spelling.")
        		return
        	}
        	DBFField[] dbfFields = table2.getAllFields()
			if (dbfFields[foreignKey].getDataType() == DBFDataType.NUMERIC || 
        	     dbfFields[foreignKey].getDataType() == DBFDataType.FLOAT) {
        	    pluginHost.showFeedback("The primary and foreign keys must be text strings and not numeric.\nYou can convert the keys using the Field Calculator in the Attribute Table Viewer.")
        		return
        	}

			HashMap<String, Object[]> hm = new HashMap<String, Object[]>()

			numFeatures = table2.getNumberOfRecords()
			oldProgress = -1
			for (int i in 0..<numFeatures) {
				rec = table2.getRecord(i);
				String keyValue = rec[foreignKey].toString().trim()
				Object[] values = new Object[numIncludedFields]
				for (int k in 0..<numIncludedFields) {
					values[k] = rec[includedFields[k]]
				}
				hm.put(keyValue, values)
			}
			
			AttributeTable table = new AttributeTable(databaseFile1)
			
			int initialFieldCount = table.getFieldCount()

			// find the primary key
			int primaryKey = table.getFieldColumnNumberFromName(primaryKeyString)
        	
        	if (primaryKey == null || primaryKey < 0) {
        		pluginHost.showFeedback("Could not locate the primary key (unique ID field). Check your spelling.")
        		return
        	}
        	dbfFields = table.getAllFields()
			if (dbfFields[primaryKey].getDataType() == DBFDataType.NUMERIC || 
        	     dbfFields[primaryKey].getDataType() == DBFDataType.FLOAT) {
        	    pluginHost.showFeedback("The primary and foreign keys must be text strings and not numeric.\nYou can convert the primary key using the Field Calculator in the Attribute Table Viewer.")
        		return
        	}
//			int primaryKey = -1
//
//			String[] fieldName = table.getAttributeTableFieldNames()
//			int f = 0
//			fieldName.each() {
//				String str = (String)it
//				if (((String)it).equals(primaryKeyString)) {
//					primaryKey = f
//				}
//				f++
//			}
//			
//			if (primaryKey == -1) {
//				pluginHost.showFeedback("Could not locate the primary key (unique ID field).")
//				return
//			}

			// append the include fields to the table
			int[] outputFieldNums = new int[numIncludedFields]
			for (int i in 0..<numIncludedFields) {
				table.addField(includedDBFFields[i])
				outputFieldNums[i] = initialFieldCount + i 
			}
			
			numFeatures = table.getNumberOfRecords()
			oldProgress = -1
			for (int i in 0..<numFeatures) {
				rec = table.getRecord(i);
				String keyValue = rec[primaryKey].toString().trim()
				if (hm.containsKey(keyValue)) {
					Object[] values = hm.get(keyValue)
					for (int k in 0..<numIncludedFields) {
						rec[initialFieldCount + k] = values[k]
					}
					table.updateRecord(i, rec)
				}

				progress = (int)(100f * i / numFeatures)
            	if (progress != oldProgress) {
					pluginHost.updateProgress(progress)
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
	def tdf = new JoinTables(pluginHost, args, name, descriptiveName)
}
