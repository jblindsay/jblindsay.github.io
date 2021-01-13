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
def name = "MergeTableWithCSV"
def descriptiveName = "Merge Table With CSV"
def description = "Merges an attribute table with a CSV text file."
def toolboxes = ["DatabaseTools"]

public class MergeTableWithCSV implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public MergeTableWithCSV(WhiteboxPluginHost pluginHost, 
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
            DialogFile df = sd.addDialogFile("Input comma separated values (CSV) file", "Input Comma Separated Values (CSV) File:", "open", "Comma-separated Values Files (*.csv; *.txt), CSV, TXT", true, false)
            DialogCheckBox headerCB = sd.addDialogCheckBox("Does the first line contain column header data?", "Does the first line contain column header data?", false)
			String[] listItems = []
            DialogComboBox foreignKey = sd.addDialogComboBox("Select the field that is the 'Foreign Key' or unique ID field.", "'Foreign Key' or CSV ID field", listItems, 0)
			DialogList includedFields = sd.addDialogList("Select the fields to include in the merge.", "Select the CSV's fields to include in the merge.", listItems, true)

            def lstrCSVFile = { evt -> if (evt.getPropertyName().equals("value")) { 
            		String csvFile = df.getValue()
            		if (csvFile != null && !csvFile.isEmpty()) { 
            			def file = new File(csvFile)
            			if (file.exists()) {
	            			boolean headerValue = Boolean.parseBoolean(headerCB.getValue())
	            			boolean firstLineRead = false
	            			String delimiter = ","
	            			file.eachLine {
								if (headerValue && !firstLineRead) {
									String[] entry = ((String)(it)).split(delimiter)
									if (entry.length == 1) {
										delimiter = "\t"
										entry = ((String)(it)).split(delimiter)
										if (entry.length == 1) {
											delimiter = " "
											entry = ((String)(it)).split(delimiter)
											
										}
									}
									foreignKey.setListItems(entry)
									includedFields.setListItems(entry)
									firstLineRead = true
								} else if (!firstLineRead) {
									String[] entry = ((String)(it)).split(delimiter)
									if (entry.length == 1) {
										delimiter = "\t"
										entry = ((String)(it)).split(delimiter)
										if (entry.length == 1) {
											delimiter = " "
											entry = ((String)(it)).split(delimiter)
											
										}
									}
									int numEntries = entry.length
									String[] colHeaders = new String[numEntries]
									int i = 0
									entry.each() { 
										colHeaders[i] = "COLUMN${i + 1}" 
										i++
									}
									foreignKey.setListItems(colHeaders)
									includedFields.setListItems(colHeaders)
									firstLineRead = true
								}
	            			}
            			}
            		}
            	}
            } as PropertyChangeListener
            df.addPropertyChangeListener(lstrCSVFile)


            def lstrHeader = { evt -> if (evt.getPropertyName().equals("value")) { 
            		String value = headerCB.getValue()
            		if (value != null && !value.isEmpty()) { 
            			String csvFile = df.getValue()
	            		if (csvFile != null && !csvFile.isEmpty()) { 
	            			def file = new File(csvFile)
	            			if (file.exists()) {
		            			boolean headerValue = Boolean.parseBoolean(headerCB.getValue())
		            			boolean firstLineRead = false
		            			String delimiter = ","
		            			file.eachLine {
									if (headerValue && !firstLineRead) {
										String[] entry = ((String)(it)).split(delimiter)
										if (entry.length == 1) {
											delimiter = "\t"
											entry = ((String)(it)).split(delimiter)
											if (entry.length == 1) {
												delimiter = " "
												entry = ((String)(it)).split(delimiter)
												if (entry.length == 1) {
													delimiter = ";"
													entry = ((String)(it)).split(delimiter)
												
												}
											}
										}
										foreignKey.setListItems(entry)
										includedFields.setListItems(entry)
										firstLineRead = true
									} else if (!firstLineRead) {
										String[] entry = ((String)(it)).split(delimiter)
										if (entry.length == 1) {
											delimiter = "\t"
											entry = ((String)(it)).split(delimiter)
											if (entry.length == 1) {
												delimiter = " "
												entry = ((String)(it)).split(delimiter)
												
											}
										}
										int numEntries = entry.length
										String[] colHeaders = new String[numEntries]
										int i = 0
										entry.each() { 
											colHeaders[i] = "COLUMN${i + 1}" 
											i++
										}
										foreignKey.setListItems(colHeaders)
										includedFields.setListItems(colHeaders)
										firstLineRead = true
									}
		            			}
	            			}
	            		}
            		}
            	} 
            } as PropertyChangeListener
            headerCB.addPropertyChangeListener(lstrHeader)


			// resize the dialog to the standard size and display it
			sd.setSize(800, 400)
			sd.visible = true
		}
	}

	//@CompileStatic
	private void execute(String[] args) {
		try {
	  		int progress, oldProgress, line
	  		String delimiter = ","
	  		
			if (args.length != 5) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			
			// read the input parameters
			String[] inputData = args[0].split(";")
			if (inputData.length < 2) {
				pluginHost.showFeedback("Error reading one of the input parameters.")
				return
			}
			String databaseFile = inputData[0].trim().replace(".shp", ".dbf")
			String primaryKeyString = inputData[1].trim()
			String csvFile = args[1]
			boolean containsHeader = Boolean.parseBoolean(args[2])
			String foreignKeyString = args[3]
			String[] includeFieldString = args[4].split(";")
			
			int numIncludedFields = includeFieldString.length
			
			def file = new File(csvFile)
			if (!file.exists()) {
				pluginHost.showFeedback("The input CSV file could not be located.")
				return
			}

			// Count the number of lines in the file. This will be used to update progress.
			int numLinesInFile = 0
			boolean delimiterFound = false
			file.eachLine { 
				numLinesInFile++  
				if (!delimiterFound) {
					String[] entry = ((String)(it)).split(delimiter)
					if (entry.length == 1) {
						delimiter = "\t"
						entry = ((String)(it)).split(delimiter)
						if (entry.length == 1) {
							delimiter = " "
							entry = ((String)(it)).split(delimiter)
							if (entry.length == 1) {
								delimiter = ";"
								entry = ((String)(it)).split(delimiter)
							
							}
						}
					}
					delimiterFound = true
				}
			}
			
			// get the field numbers of the foreign key and include fields
			int foreignKey = -1
			int[] includedFields = new int[numIncludedFields]
			if (!containsHeader) {
				int k = 0
				includeFieldString.each() {
					String str = ((String)it).toLowerCase().replace("column", "").replace("col","")
					if (str.isNumber()) {
						includedFields[(int)k] = Integer.parseInt(str.trim()) - 1
					}
					k++
				}
				String str = ((String)foreignKeyString).toLowerCase().replace("column", "").replace("col","")
				if (str.isNumber()) {
					foreignKey = Integer.parseInt(str.trim()) - 1
				}
			} else {
				boolean firstLineRead = false
    			file.eachLine {
					String[] entry = ((String)(it)).split(delimiter)
					if (!firstLineRead) {
						for (int b in 0..<entry.length) {
							for (int a in 0..<numIncludedFields) {
								if (entry[b].trim().equals(includeFieldString[a])) {
									includedFields[a] = b
								}
							}
							if (entry[b].trim().equals(foreignKeyString)) {
								foreignKey = b
							}
						}
						firstLineRead = true
					}
    			}
			}
			
			if (foreignKey == -1) {
				pluginHost.showFeedback("Could not locate the foreign key (unique ID field).")
				return
			}
			
			
			String[] includedFieldNames = new String[numIncludedFields]
			for (int i in 0..<numIncludedFields) {
				includedFieldNames[i] = "NEWFIELD${i + 1}"
			}
			int[] fieldDataTypes = new int[numIncludedFields]
			int[] fieldMaxLength = new int[numIncludedFields]
			int[] fieldMaxDecimalCount = new int[numIncludedFields]
			HashMap<String, Object[]> hm = new HashMap<String, Object[]>()
			boolean checkedFirstLine = !containsHeader
			boolean checkedDataType = false
			line = 0
			oldProgress = -1
			file.eachLine {
				String[] entry = ((String)(it)).split(delimiter)
				if (checkedFirstLine) {
					String keyValue = entry[foreignKey].trim()
					Object[] values = new Object[numIncludedFields]
					for (int i in 0..<numIncludedFields) {
						String str = entry[includedFields[i]].replace("\"", "")
						if (!checkedDataType) {
							if (str.isNumber()) {
			            		fieldDataTypes[i] = 0 //"Double"
			            	} else if (str.toLowerCase().equals("true") || 
			            	  str.toLowerCase().equals("false") ||
			            	  str.toLowerCase().equals("yes") ||
			            	  str.toLowerCase().equals("no")) {
			            		fieldDataTypes[i] = 1 //"Boolean"
			            	} else {
			            		fieldDataTypes[i] = 2 //"String"
			            	}
						}
						switch (fieldDataTypes[i]) {
							case 0: // Double
								values[i] = Double.parseDouble(str)
								int j = str.lastIndexOf(".")
								if (j >= 0) {
									int numDec = str.length() - j - 1
									if (numDec > fieldMaxDecimalCount[i]) { fieldMaxDecimalCount[i] = numDec }
								}
								break
							case 1: // Boolean
								values[i] = Boolean.parseBoolean(str)
								break
							case 2: // String
								values[i] = str
								break
						}
						if (fieldMaxLength[i] < str.length()) { fieldMaxLength[i] = str.length() }
					}
					hm.put(keyValue, values)
					checkedDataType = true
				} else {
					// this will only occur if the file contains a header
					for (int i in 0..<numIncludedFields) {
						includedFieldNames[i] = entry[includedFields[i]].trim()
					}
				}
				checkedFirstLine = true

				line++
				progress = (int)(100f * line / numLinesInFile)
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
			
			AttributeTable table = new AttributeTable(databaseFile)
			
			int initialFieldCount = table.getFieldCount()

			// find the primary key
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

			DBFField[] fields = table.getAllFields()
        	int primaryKey = table.getFieldColumnNumberFromName(primaryKeyString)
        	
        	if (primaryKey == null || primaryKey < 0) {
        		pluginHost.showFeedback("Could not locate the primary key (unique ID field). Check your spelling.")
        		return
        	}
			if (fields[primaryKey].getDataType() == DBFDataType.NUMERIC || 
        	     fields[primaryKey].getDataType() == DBFDataType.FLOAT) {
        	    pluginHost.showFeedback("The primary and foreign keys must be text strings and not numeric.\nYou can convert the primary key using the Field Calculator in the Attribute Table Viewer.")
        		return
        	}
			
			// append the include fields to the table
			int[] outputFieldNums = new int[numIncludedFields]
			for (int i in 0..<numIncludedFields) {
				DBFField field = new DBFField()
				field.setName(includedFieldNames[i])
				switch (fieldDataTypes[i]) {
					case 0: // Numeric
						field.setDataType(DBFField.DBFDataType.FLOAT)
						field.setDecimalCount(fieldMaxDecimalCount[i])
						break
					case 1: // Boolean
						field.setDataType(DBFField.DBFDataType.BOOLEAN)
						break
					case 2: // String
						field.setDataType(DBFField.DBFDataType.STRING)
						break
				}
				field.setFieldLength(fieldMaxLength[i])
						
				table.addField(field)
				outputFieldNums[i] = initialFieldCount + i 
			}
			
			int numFeatures = table.getNumberOfRecords()
			Object[] rec;
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
	def tdf = new MergeTableWithCSV(pluginHost, args, name, descriptiveName)
}
