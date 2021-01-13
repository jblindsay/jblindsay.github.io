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
import java.util.Arrays
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
import whitebox.geospatialfiles.LASReader
import whitebox.geospatialfiles.LASReader.PointRecord

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
def name = "LidarHistogram"
def descriptiveName = "LiDAR Histogram"
def description = "Creates a histogram for the attributes of a LiDAR (LAS) file."
def toolboxes = ["LidarTools"]

public class LidarHistogram implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName

    public LidarHistogram(WhiteboxPluginHost pluginHost, 
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
        	DialogFile df = sd.addDialogFile("Select the input LiDAR (.las) file", "Input LAS File:", "open", "LAS Files (*.las), LAS", true, false)
        	//DialogFieldSelector dfs = sd.addDialogFieldSelector("Input file and attribute field.", "Input Attribute Field:", false)
            String[] listItems = ["Elevation", "Intensity", "Scan Angle"]
            DialogComboBox zKey = sd.addDialogComboBox("Select the attribute to plot.", "Attribute", listItems, 0)
			DialogDataInput din = sd.addDialogDataInput("Enter the number of bins", "Number of Bins:", "", true, true)
            def lstr = { evt -> if (evt.getPropertyName().equals("value")) { 
            		String value = df.getValue()
            		if (value != null && !value.isEmpty()) { 
            			if ((new File(value)).exists()) {
            				LASReader las = new LASReader(value);
	        				int numPoints = (int) las.getNumPointRecords();
							int numBins = (int)(Math.ceil(Math.log(numPoints) / Math.log(2) + 1))
							din.setValue(numBins.toString())
            			}
            		}
            	}
            } as PropertyChangeListener
            df.addPropertyChangeListener(lstr)

            // resize the dialog to the standard size and display it
            sd.setSize(800, 400)
            sd.visible = true
        }
    }

    @CompileStatic
    private void execute(String[] args) {
        try {

			if (args.length < 2) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			PointRecord point;
			
			String inputFile = args[0].trim();
			LASReader las = new LASReader(inputFile);
	        int numPoints = (int) las.getNumPointRecords();
	        String fieldName = args[1].trim();
			int numBins = (int)(Math.ceil(Math.log(numPoints) / Math.log(2) + 1))
			if (args[2].isNumber()) {
				numBins = Integer.parseInt(args[2])
			}

			
			double[] data1 = new double[numPoints]
			int progress
			int oldProgress = -1
			
			JFreeChart chart
			
			double minVal = Double.POSITIVE_INFINITY
			double maxVal = Double.NEGATIVE_INFINITY

			if (fieldName.toLowerCase().contains("elev")) {
				for (int a = 0; a < numPoints; a++) {
					point = las.getPointRecord(a);
	                if (!point.isPointWithheld()) {
	                    double val = point.getZ();
	                    data1[a] = val;
						if (val < minVal) minVal = val;
						if (val > maxVal) maxVal = val;
	                }
					progress = (int)(100f * a / (numPoints - 1))
					if (progress != oldProgress) {
						pluginHost.updateProgress(progress)
						oldProgress = progress
						if (pluginHost.isRequestForOperationCancelSet()) {
	                        pluginHost.showFeedback("Operation cancelled")
							return
	                    }
					}
				}
			} else if (fieldName.toLowerCase().contains("int")) {
				for (int a = 0; a < numPoints; a++) {
					point = las.getPointRecord(a);
	                if (!point.isPointWithheld()) {
	                    double val = point.getIntensity();
	                    data1[a] = val;
						if (val < minVal) minVal = val;
						if (val > maxVal) maxVal = val;
	                }
					progress = (int)(100f * a / (numPoints - 1))
					if (progress != oldProgress) {
						pluginHost.updateProgress(progress)
						oldProgress = progress
						if (pluginHost.isRequestForOperationCancelSet()) {
	                        pluginHost.showFeedback("Operation cancelled")
							return
	                    }
					}
				}
			} else if (fieldName.toLowerCase().contains("scan")) {
				for (int a = 0; a < numPoints; a++) {
					point = las.getPointRecord(a);
	                if (!point.isPointWithheld()) {
	                    double val = point.getScanAngle();
	                    data1[a] = val;
						if (val < minVal) minVal = val;
						if (val > maxVal) maxVal = val;
	                }
					progress = (int)(100f * a / (numPoints - 1))
					if (progress != oldProgress) {
						pluginHost.updateProgress(progress)
						oldProgress = progress
						if (pluginHost.isRequestForOperationCancelSet()) {
	                        pluginHost.showFeedback("Operation cancelled")
							return
	                    }
					}
				}
			} else {
				pluginHost.showFeedback("The plotting attribute has not been specified correctly.\nTry either 'Elevation', 'Intensity', or 'Scan Angle'.");
				return 
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

			Arrays.parallelSort(data1)

			DecimalFormat df = new DecimalFormat("###,###,###,###.###")
			DecimalFormat df2 = new DecimalFormat("##0.0##%")
			
			StringBuilder ret = new StringBuilder()
			ret.append("<!DOCTYPE html>")
			ret.append('<html lang="en">')
			ret.append("<head>")
	        ret.append("<title>LiDAR Histogram Summary Statistics</title>").append("\n")
	        ret.append("<style>")
			ret.append("table {margin-left: 15px;} ")
			ret.append("h1 {font-size: 14pt; margin-left: 15px; margin-right: 15px; text-align: center; font-family: Helvetica, Verdana, Geneva, Arial, sans-serif;} ")
			ret.append("p {font-size: 12pt; font-family: Helvetica, Verdana, Geneva, Arial, sans-serif; margin-left: 15px; margin-right: 15px;} ")
			ret.append("table {font-size: 12pt; font-family: Helvetica, Verdana, Geneva, Arial, sans-serif;}")
			ret.append("table th {border-width: 1px; padding: 8px; border-style: solid; border-color: #666666; background-color: #dedede; }")
			ret.append("table td {border-width: 1px; padding: 8px; border-style: solid; border-color: #666666; background-color: #ffffff; }")
			ret.append("caption {font-family: Helvetica, Verdana, Geneva, Arial, sans-serif; margin-left: 15px; margin-right: 15px;} ")
			ret.append(".numberCell { text-align: right; }") 
            ret.append("</style></head>").append("\n")
            ret.append("<body><h1>LiDAR Histogram Summary Statistics</h1>").append("\n")

            
			ret.append("<p><b>File Name:</b> &nbsp ").append(new File(inputFile).getName()).append("</p>\n")

			ret.append("<br><table border=\"1\" cellspacing=\"0\" cellpadding=\"3\">").append("\n")
			ret.append("<caption>Percentile Statistics</caption>")
			
			ret.append("<tr><th>Percentile</th><th>Value</th></tr>")

			long[] pointsByRet = las.getNumPointsByReturn()
			double percentile = 0;
			for (int n = 0; n < 21; n++) {
				int entryNum = (int)(percentile / 100.0 * (numPoints - 1));
				if (entryNum < numPoints) {
					double value = data1[entryNum]
					ret.append("<tr><td class=\"numberCell\">").append(percentile)
					ret.append("</td><td class=\"numberCell\">")
					ret.append(df.format(value))
					ret.append("</td></tr>").append("\n")
					percentile += 5;
				}
			}
			
			ret.append("</table>")

			ret.append("</body></html>")
			pluginHost.returnData(ret.toString());
			
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
    def f = new LidarHistogram(pluginHost, args, name, descriptiveName)
}
