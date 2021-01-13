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
import whitebox.geospatialfiles.VectorLayerInfo
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.geospatialfiles.shapefile.attributes.DBFField.DBFDataType
import whitebox.geospatialfiles.shapefile.ShapeFileRecord
import groovy.transform.CompileStatic

import java.awt.Color
import java.awt.Dimension
import java.awt.Font
import java.awt.RenderingHints
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.attributes.DBFField
import whitebox.geospatialfiles.shapefile.attributes.*
import org.jfree.chart.title.LegendTitle
import org.jfree.chart.block.BlockBorder
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.IntervalMarker;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.Layer;
import org.jfree.ui.RectangleAnchor;
import org.jfree.ui.RefineryUtilities;
import org.jfree.ui.TextAnchor;
import org.jfree.chart.*;
import org.jfree.data.statistics.*;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CategoryPlot;

import org.jfree.chart.plot.PlotOrientation;

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "AttributeScattergram"
def descriptiveName = "Attribute Scattergram"
def description = "Creates a scattergram from two attributes in a vector's attribute table."
def toolboxes = ["DatabaseTools", "StatisticalTools"]

public class AttributeScattergram implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName

    public AttributeScattergram(WhiteboxPluginHost pluginHost, 
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
        	DialogFieldSelector dfs = sd.addDialogFieldSelector("Input file and two attribute fields.", "Select 2 Attribute Fields:", true)

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
			if (inputData.length != 3) {
				pluginHost.showFeedback("Error reading one of the input parameters.")
				return
			}
			String inputFile = inputData[0].trim()
			
        	ShapeFile shape = new ShapeFile(inputFile)
			
			String fieldName1 = inputData[1].trim()
			String fieldName2 = inputData[2].trim()

			AttributeTable table = shape.getAttributeTable()
			int numRecords = shape.getNumberOfRecords()

			// see if the specified field is numerical or categorical
			boolean isAttributeNumeric = true
			DBFField[] fields = table.getAllFields()
        	def fieldNum1 = table.getFieldColumnNumberFromName(fieldName1)
        	if (fieldNum1 == null || fieldNum1 < 0) {
        		pluginHost.showFeedback("Could not locate the specified field in the attribute table. Check your spelling.")
        		return
        	}
			if (fields[fieldNum1].getDataType() != DBFDataType.NUMERIC && 
        	     fields[fieldNum1].getDataType() != DBFDataType.FLOAT) {
        	    pluginHost.showFeedback("Both of the input attribute fields must be numerical.")
        		return
        	}

        	def fieldNum2 = table.getFieldColumnNumberFromName(fieldName2)
        	if (fieldNum2 == null || fieldNum2 < 0) {
        		pluginHost.showFeedback("Could not locate the specified field in the attribute table. Check your spelling.")
        		return
        	}
			if (fields[fieldNum2].getDataType() != DBFDataType.NUMERIC && 
        	     fields[fieldNum2].getDataType() != DBFDataType.FLOAT) {
        	    pluginHost.showFeedback("Both of the input attribute fields must be numerical.")
        		return
        	}
			
			XYSeriesCollection xyCollection = new XYSeriesCollection();
        	XYSeries series = new XYSeries("Data")
        	//Object[] rowData
			int progress
			int oldProgress = -1
			
			int numBins = (int)(Math.ceil(Math.log(numRecords) / Math.log(2) + 1))
			if (args.length >= 2 && args[1].isNumber()) {
				numBins = Integer.parseInt(args[1])
			}
		
//			double minVal = Double.POSITIVE_INFINITY
//			double maxVal = Double.NEGATIVE_INFINITY
		
			for (int rec in 0..<numRecords) {
				double x = (double)table.getValue(rec, fieldName1)
				double y = (double)table.getValue(rec, fieldName2)

				series.add(x, y)
				
//				if (x < minVal) minVal = x
//				if (y > maxVal) maxVal = y
			
				progress = (int)(100f * rec / (numRecords - 1))
				if (progress != oldProgress) {
					pluginHost.updateProgress(progress)
					oldProgress = progress
				}
			}

			xyCollection.addSeries(series)
			
			// create the chart...
			JFreeChart chart = ChartFactory.createScatterPlot(
				"",
	            fieldName1, fieldName2, 
	            xyCollection, 
	            PlotOrientation.VERTICAL,
	            false, 
	            true, 
	            false
	        )

        	XYPlot plot = (XYPlot)chart.getPlot();
			plot.setBackgroundPaint(Color.WHITE); 
			plot.setDomainGridlinePaint(Color.BLACK); 
			plot.setRangeGridlinePaint(Color.BLACK); 
			
			ChartPanel chartPanel = new ChartPanel(chart);
			chartPanel.setFillZoomRectangle(true);
			chartPanel.setMouseWheelEnabled(true);
			chartPanel.setPreferredSize(new Dimension(520, 470))
			
			pluginHost.returnData(chartPanel)

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
    def f = new AttributeScattergram(pluginHost, args, name, descriptiveName)
}
