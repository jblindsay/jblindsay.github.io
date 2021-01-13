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
import org.jfree.data.xy.IntervalXYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.Layer;
import org.jfree.ui.RectangleAnchor;
import org.jfree.ui.RefineryUtilities;
import org.jfree.ui.TextAnchor;
import org.jfree.chart.*;
import org.jfree.data.statistics.*;
import org.jfree.chart.renderer.xy.StandardXYBarPainter;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CategoryPlot;

import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.category.BarRenderer;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.chart.renderer.category.CategoryItemRenderer
import org.jfree.chart.renderer.category.BarRenderer
import org.jfree.chart.renderer.category.StandardBarPainter

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "AttributeHistogram"
def descriptiveName = "Attribute Histogram"
def description = "Creates a histogram for an attribute in a vector's attribute table."
def toolboxes = ["DatabaseTools", "StatisticalTools"]

public class AttributeHistogram implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName

    public AttributeHistogram(WhiteboxPluginHost pluginHost, 
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
            DialogDataInput din = sd.addDialogDataInput("Enter the number of bins", "Number of Bins:", "", true, true)
			din.setVisible(false)
            def lstr = { evt -> if (evt.getPropertyName().equals("value")) { 
            		String value = dfs.getValue()
            		if (value != null && !value.isEmpty()) { 
            			String[] inputData = value.split(";")
        				if (inputData.length > 1 && inputData[0] != null && !inputData[0].isEmpty()) {
        					if (inputData[1] != null && !inputData[1].isEmpty()) {
	        					String inputFile = inputData[0].trim()
	        					String fieldName = inputData[1].trim()
	        					if ((new File(inputFile)).exists()) {
	    							ShapeFile shape = new ShapeFile(inputFile)
	    							
	    							AttributeTable table = shape.getAttributeTable()
									
									// see if the specified field is numerical or categorical
									boolean isAttributeNumeric = true
									DBFField[] fields = table.getAllFields()
						        	def fieldNum = table.getFieldColumnNumberFromName(fieldName)
									if (fields[fieldNum].getDataType() != DBFDataType.NUMERIC && 
						        	     fields[fieldNum].getDataType() != DBFDataType.FLOAT) {
						        	    isAttributeNumeric = false
						        	}
						        	if (isAttributeNumeric) {
		    							int numRecords = shape.getNumberOfRecords()
		    							int numBins = (int)(Math.ceil(Math.log(numRecords) / Math.log(2) + 1))
										din.setValue(numBins.toString())
										din.setVisible(true)
									} else {
										din.setValue("")
										din.setVisible(false)
									}
	        					}
        					}
        				}
            		}
            	}
            } as PropertyChangeListener
            dfs.addPropertyChangeListener(lstr)

            DialogCheckBox cb = sd.addDialogCheckBox("Plot as a cumulative distribution function (CDF)", "Plot as a cumulative distribution function (CDF)?", false)
            def lstr2 = { evt -> if (evt.getPropertyName().equals("value")) { 
            	if (Boolean.parseBoolean(cb.getValue())) {
            		din.setVisible(false);
            	} else {
            		String[] inputData = dfs.getValue().split(";")
    				if (inputData.length > 1 && inputData[0] != null && !inputData[0].isEmpty()) {
    					if (inputData[1] != null && !inputData[1].isEmpty()) {
		            		String inputFile = inputData[0].trim()
							String fieldName = inputData[1].trim()
							if ((new File(inputFile)).exists()) {
								ShapeFile shape = new ShapeFile(inputFile)
								
								AttributeTable table = shape.getAttributeTable()
								
								// see if the specified field is numerical or categorical
								boolean isAttributeNumeric = true
								DBFField[] fields = table.getAllFields()
					        	def fieldNum = table.getFieldColumnNumberFromName(fieldName)
								if (fields[fieldNum].getDataType() != DBFDataType.NUMERIC && 
					        	     fields[fieldNum].getDataType() != DBFDataType.FLOAT) {
					        	    isAttributeNumeric = false
					        	}
					        	println isAttributeNumeric
					        	if (isAttributeNumeric) {
									int numRecords = shape.getNumberOfRecords()
									int numBins = (int)(Math.ceil(Math.log(numRecords) / Math.log(2) + 1))
									din.setValue(numBins.toString())
									din.setVisible(true)
								} else {
									din.setValue("")
									din.setVisible(false)
								}
							}
    					}
    				}
            	}
            } } as PropertyChangeListener
            cb.addPropertyChangeListener(lstr2)
			
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
			String inputFile = inputData[0].trim();
			boolean outputCDF = false;
			if (args.length >= 3) {
				if (!args[2].trim().toLowerCase().isEmpty()) {
					outputCDF = Boolean.parseBoolean(args[2])
				}
			}
			
        	ShapeFile shape = new ShapeFile(inputFile)
			
			String fieldName = inputData[1]

			AttributeTable table = shape.getAttributeTable()
			int numRecords = shape.getNumberOfRecords()

			// see if the specified field is numerical or categorical
			boolean isAttributeNumeric = true
			DBFField[] fields = table.getAllFields()
        	def fieldNum = table.getFieldColumnNumberFromName(fieldName)
        	if (fieldNum == null || fieldNum < 0) {
        		pluginHost.showFeedback("Could not locate the specified field in the attribute table. Check your spelling.")
        		return
        	}
			if (fields[fieldNum].getDataType() != DBFDataType.NUMERIC && 
        	     fields[fieldNum].getDataType() != DBFDataType.FLOAT) {
        	    isAttributeNumeric = false
        	}
			
			double[] data1 = new double[numRecords]
			Object[] rowData
			int progress
			int oldProgress = -1
			
			JFreeChart chart
			
			if (isAttributeNumeric) {
				int numBins = (int)(Math.ceil(Math.log(numRecords) / Math.log(2) + 1))
				if (args.length >= 2 && args[1].isNumber()) {
					numBins = Integer.parseInt(args[1])
				}
			
				double minVal = Double.POSITIVE_INFINITY
				double maxVal = Double.NEGATIVE_INFINITY
			
				for (int rec in 0..<numRecords) {
					double val = (double)table.getValue(rec, fieldName)
					data1[rec] = val
					if (val < minVal) minVal = val
					if (val > maxVal) maxVal = val
				
					progress = (int)(100f * rec / (numRecords - 1))
					if (progress != oldProgress) {
						pluginHost.updateProgress(progress)
						oldProgress = progress
					}
				}

				HistogramDataset dataset = new HistogramDataset();
				dataset.setType(HistogramType.RELATIVE_FREQUENCY);
				
				dataset.addSeries(fieldName, data1, numBins);
				
				// create the chart...
				chart = ChartFactory.createHistogram(
					"",
					fieldName,
					"Relative Frequency", 
					dataset,
					PlotOrientation.VERTICAL,
					false,
					true,
					false
				);

				// get a reference to the plot for further customisation...
				XYPlot plot = (XYPlot) chart.getPlot();
				plot.setBackgroundPaint(Color.WHITE);
				
				XYBarRenderer renderer = (XYBarRenderer)plot.getRenderer();
				renderer.setShadowVisible(false);
				renderer.setBarPainter(new StandardXYBarPainter()); 
				renderer.setDrawBarOutline(true);
				renderer.setSeriesOutlinePaint(0, Color.black);
				renderer.setSeriesPaint(0, Color.lightGray);
							
				NumberAxis domain = (NumberAxis) plot.getDomainAxis()
				domain.setRange(minVal, maxVal)
							
			} else {
				if (outputCDF) {
					pluginHost.showFeedback("The attribute is categorial and therefore \ncannot be plotted as a CDF.");
				}
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

				DefaultCategoryDataset dataset = new DefaultCategoryDataset();
				for (Object k : hm.keySet()) {
					int val = (hm.get(k)).intValue()
					String key = k.toString()
					if (k == null || key.trim().isEmpty()) key = "null"
					dataset.addValue(val, "data", key.trim())
				}

				chart = ChartFactory.createBarChart(
		            "",       // chart title
		            "Category",               // domain axis label
		            "Frequency",                  // range axis label
		            dataset,                  // data
		            PlotOrientation.VERTICAL, // orientation
		            false,
		            true,                     // tooltips?
		            false                     // URLs?
		        );

		        CategoryPlot plot = (CategoryPlot) chart.getPlot();

				plot.setBackgroundPaint(Color.WHITE);
				
				CategoryItemRenderer renderer = (CategoryItemRenderer)plot.getRenderer();
				renderer.setSeriesOutlinePaint(0, Color.black);
				renderer.setSeriesPaint(0, Color.lightGray);
				renderer.setSeriesOutlineStroke(0, new BasicStroke(0.5))

				BarRenderer renderer2 = (BarRenderer) plot.getRenderer();
				renderer2.setDrawBarOutline(true);
				renderer2.setBarPainter(new StandardBarPainter())
				
			}
			
			Font font = new Font("Sans-serif", Font.PLAIN, 12);
			StandardChartTheme chartTheme = new StandardChartTheme("");
			chartTheme.setExtraLargeFont(font.deriveFont(12f));
			chartTheme.setLargeFont(font.deriveFont(12f));
			chartTheme.setRegularFont(font.deriveFont(12f));
			chartTheme.setSmallFont(font.deriveFont(12f));
			ChartFactory.setChartTheme(chartTheme);
			
			
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
    def f = new AttributeHistogram(pluginHost, args, name, descriptiveName)
}
