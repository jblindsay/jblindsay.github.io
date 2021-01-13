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
import java.awt.*
import java.text.DecimalFormat
import java.io.File
import java.util.Date
import java.util.ArrayList
import java.util.concurrent.atomic.AtomicInteger
import java.beans.PropertyChangeEvent
import java.beans.PropertyChangeListener
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.ui.plugin_dialog.*
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.geospatialfiles.shapefile.attributes.DBFField.DBFDataType
import whitebox.geospatialfiles.shapefile.ShapeFileRecord
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "ListUniqueValues"
def descriptiveName = "List Unique Values"
def description = "Lists the unique values contained in a field witin a vector's attribute table."
def toolboxes = ["DatabaseTools", "StatisticalTools"]

public class ListUniqueValues implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName

    public ListUniqueValues(WhiteboxPluginHost pluginHost, 
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
        	//DialogFile dfIn = sd.addDialogFile("Input file", "Input Vector File:", "open", "Vector Files (*.shp), SHP", true, false)
            DialogFieldSelector dfs = sd.addDialogFieldSelector("Input file and attribute field.", "Input Attribute Field:", false)
            
            // resize the dialog to the standard size and display it
            sd.setSize(800, 400)
            sd.visible = true
        }
    }

    @CompileStatic
    private void execute(String[] args) {
        try {

			if (args.length < 1) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			
			String[] inputData = args[0].split(";")
			if (inputData.length < 2) {
				pluginHost.showFeedback("Error reading one of the input parameters.")
				return
			}
			String inputFile = inputData[0].trim()
			
        	
        	ShapeFile shape = new ShapeFile(inputFile)
			
			String fieldName = inputData[1]

			AttributeTable table = shape.getAttributeTable()
			int numRecords = shape.getNumberOfRecords()

//			boolean isAttributeNumeric = true
			DBFField[] fields = table.getAllFields()
        	def fieldNum = table.getFieldColumnNumberFromName(fieldName)
        	if (fieldNum == null || fieldNum < 0) {
        		pluginHost.showFeedback("Could not locate the specified field in the attribute table. Check your spelling.")
        		return
        	}
//			if (fields[fieldNum].getDataType() != DBFDataType.NUMERIC && 
//        	     fields[fieldNum].getDataType() != DBFDataType.FLOAT) {
//        	    isAttributeNumeric = false
//        	}
			
			double[] data1 = new double[numRecords]
			Object[] rowData
			int progress
			int oldProgress = -1
			
			// there is a categorical attribute.
			Map<Object, AtomicInteger> hm = new TreeMap<Object, AtomicInteger>();
			oldProgress = -1
			for (int rec in 0..<numRecords) {
				AtomicInteger value = hm.get(table.getValue(rec, fieldName));
			    if (value == null) {
			       hm.put(table.getValue(rec, fieldName), new AtomicInteger(1))
				} else {
			        value.incrementAndGet()
				}
				progress = (int)(100f * rec / (numRecords - 1))
				if (progress != oldProgress) {
					pluginHost.updateProgress(progress)
					oldProgress = progress
				}
			}

			StringBuilder ret = new StringBuilder()
			ret.append("<!DOCTYPE html>")
			ret.append('<html lang="en">')
			
			ret.append("<head>")
			ret.append("<title>List Unique Values</title>").append("\n")
			
			ret.append("<style  type=\"text/css\">")
			ret.append("table {margin-left: 15px;} ")
			ret.append("h1 {font-size: 14pt; margin-left: 15px; margin-right: 15px; text-align: center; font-family: Helvetica, Verdana, Geneva, Arial, sans-serif;} ")
			ret.append("p {font-size: 12pt; font-family: Helvetica, Verdana, Geneva, Arial, sans-serif; margin-left: 15px; margin-right: 15px;} ")
			ret.append("table {font-size: 12pt; font-family: Helvetica, Verdana, Geneva, Arial, sans-serif;}")
			ret.append("table th {border-width: 1px; padding: 8px; border-style: solid; border-color: #666666; background-color: #dedede; }")
			ret.append("table td {border-width: 1px; padding: 8px; border-style: solid; border-color: #666666; background-color: #ffffff; }")
			ret.append("caption {font-family: Helvetica, Verdana, Geneva, Arial, sans-serif; margin-left: 15px; margin-right: 15px;} ")
			ret.append(".numberCell { text-align: right; }") 
			ret.append("</style></head>").append("\n")
			ret.append("<body><h1>List Unique Values</h1>").append("\n")

			ret.append("<p>Input file: ${inputFile}</p>")
			ret.append("<p>Field Name: ${fieldName}</p>")
			
			ret.append("<p><table border=\"1\" cellspacing=\"0\" cellpadding=\"3\">").append("\n")
			ret.append("<caption>Category Data</caption>")
			ret.append("<tr><th>Category</th><th>Frequency</th></tr>")
			
			for (Object k : hm.keySet()) {
				String val = String.valueOf((hm.get(k)).intValue())
				String key = k.toString()
				if (k == null || key.trim().isEmpty()) key = "null"
				ret.append("<tr><td>${key}</td><td class=\"numberCell\">${val}</td></tr>")
			}

			ret.append("</table></p>")
			ret.append("</body></html>")
			
			pluginHost.returnData(ret.toString())

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
    def f = new ListUniqueValues(pluginHost, args, name, descriptiveName)
}
