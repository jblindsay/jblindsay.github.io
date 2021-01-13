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
import java.io.File
import java.util.Date
import java.util.ArrayList
import java.util.Arrays
import java.text.DecimalFormat
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.attributes.AttributeTable
import whitebox.geospatialfiles.shapefile.attributes.DBFField
import whitebox.geospatialfiles.shapefile.attributes.DBFField.DBFDataType
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.FileUtilities;
import groovy.transform.CompileStatic

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.StandardChartTheme;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.time.Month;
import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;
import org.jfree.data.xy.XYDataset;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RectangleInsets;
import org.jfree.ui.RefineryUtilities;
import org.jfree.chart.renderer.xy.XYItemRenderer
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.chart.renderer.xy.*;
import org.jfree.chart.plot.*;
import org.jfree.chart.renderer.LookupPaintScale;
import org.jfree.data.xy.DefaultXYZDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.ui.RectangleAnchor;


// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "ClusterAttributes"
def descriptiveName = "Cluster Attributes"
def description = "Performs k-means clustering on selected attributes associated with a vector file."
def toolboxes = ["VectorTools", "StatisticalTools"]

public class ClusterAttributes implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
	
    public ClusterAttributes(WhiteboxPluginHost pluginHost, 
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
            sd.addDialogFieldSelector("Input file and attributes field.", "Select the Attributes to include:", true)															
			sd.addDialogComboBox("Should the data be rescaled prior to the analysis?", "Data Rescaling Method:", ["Standardize", "Normalize", "Do Nothing"], 0)	
			sd.addDialogDataInput("Enter the number of clusters", "Number of Clusters:", "", true, false)
            sd.addDialogDataInput("Enter the maximum number of iterations", "Maximum Iterations:", "100", true, false)
            sd.addDialogDataInput("Enter the minimum percent of records changing clusters between iterations to stop.", "Minimum % Change:", "0.5", true, false)
            
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
            int i, j, m, progress, oldProgress, fieldNum
            double val, dist, minDist
            
            if (args.length < 5) {
                pluginHost.showFeedback("Incorrect number of arguments given to tool.")
                return
            }
            
            // read the input parameters
            String[] inputData = args[0].split(";")
			if (inputData[0] == null || inputData[0].isEmpty()) {
				pluginHost.showFeedback("Input shapefile and attribute not specified.")
				return
			}
			if (inputData.length < 2 || inputData[1] == null || inputData[1].isEmpty()) {
				pluginHost.showFeedback("Input shapefile and attribute not specified.")
				return
			}
			String inputFile = inputData[0]

			int numAttributes = inputData.length - 1
			if (numAttributes < 2) {
				pluginHost.showFeedback("At least two attributes are needed to perform clusering.")
				return
			}
			
            String rescaleMode = "raw"
			if (args[1].toLowerCase().contains("norm")) {
				rescaleMode = "normalize"
			} else if (args[1].toLowerCase().contains("stand")) {
				rescaleMode = "standardize"
			}

			int k = Integer.parseInt(args[2])
			if (k < 2) {
				pluginHost.showFeedback("At least two clusters are needed to perform clusering.")
				k = 2
			}
			
			int numIterations = Integer.parseInt(args[3])
			double minChange = Double.parseDouble(args[4])
			
			
            // see if each of the input attributes are numerical
            def table = new AttributeTable(inputFile.replace(".shp", ".dbf"))
            int numRecords = table.getNumberOfRecords()
            DBFField[] fields = table.getAllFields()
            int[] attributeFieldNums = new int[numAttributes]
			for (i = 0; i < numAttributes; i++) {
				fieldNum = table.getFieldColumnNumberFromName(inputData[i + 1])
				attributeFieldNums[i] = fieldNum
				if (fields[fieldNum].getDataType() != DBFDataType.NUMERIC && 
            	     fields[fieldNum].getDataType() != DBFDataType.FLOAT) {
            	    pluginHost.showFeedback("The input attributes must be of numeric type. If you have \n" +
            	                            "categorical attributes, use the Field Calculator to create new \n" + 
            	                            "Dummy variables from these attributes (i.e. convert them to \n" + 
            	                            "numerical data).")
            	    return
            	}
			}

			// read the data in
			double[][] data = new double[numRecords][numAttributes]
			
			Object[] recData
			pluginHost.updateProgress("Reading attributes:", 0)
			oldProgress = -1
			for (i = 0; i < numRecords; i++) {
				recData = table.getRecord(i)
				for (j = 0; j < numAttributes; j++) {
					data[i][j] = (double)(recData[attributeFieldNums[j]])
				}
				progress = (int)(100f * i / (numRecords - 1))
				if (progress != oldProgress) {
					pluginHost.updateProgress("Reading attributes:", progress)
					oldProgress = progress
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}

			pluginHost.updateProgress("Rescaling the data...", 0)

			double[] minVals = new double[numAttributes]
			double[] maxVals = new double[numAttributes]
			
			double[] mean = new double[numAttributes]
			double[] stDev = new double[numAttributes]
				
			// rescale the input data if necessary
			if (rescaleMode.equals("normalize")) {
				// rescales each of the attributes to have a common 0-1 scale
				for (j = 0; j < numAttributes; j++) {
					minVals[j] = Double.POSITIVE_INFINITY
					maxVals[j] = Double.NEGATIVE_INFINITY
				}
				oldProgress = -1
				for (i = 0; i < numRecords; i++) {
					for (j = 0; j < numAttributes; j++) {
						val = data[i][j]
						if (val < minVals[j]) minVals[j] = val
						if (val > maxVals[j]) maxVals[j] = val
					}
					progress = (int)(100f * i / (numRecords - 1))
					if (progress != oldProgress) {
						pluginHost.updateProgress("Rescaling the data...", progress)
						oldProgress = progress
						// check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
					}
				}

				oldProgress = -1
				for (i = 0; i < numRecords; i++) {
					for (j = 0; j < numAttributes; j++) {
						val = data[i][j]
						data[i][j] = (val - minVals[j]) / (maxVals[j] - minVals[j])
					}
					progress = (int)(100f * i / (numRecords - 1))
					if (progress != oldProgress) {
						pluginHost.updateProgress("Rescaling the data...", progress)
						oldProgress = progress
						// check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
					}
				}
				
			} else if (rescaleMode.equals("standardize")) {
				// converts the attributes to z-scores

				// first compute the means and standard deviations of each field
				double[] M2 = new double[numAttributes]
				int[] n = new int[numAttributes]
				double delta
				
			 	for (i = 0; i < numRecords; i++) {
					for (j = 0; j < numAttributes; j++) {
						n[j]++
						val = data[i][j]
						delta = val - mean[j]
						mean[j] = mean[j] + delta / n[j]
						M2[j] = M2[j] + delta * (val - mean[j])
					}
				}
				for (j = 0; j < numAttributes; j++) {
					stDev[j] = Math.sqrt(M2[j] / (n[j] - 1))
				}

				// now convert the data to z-scores
			    for (i = 0; i < numRecords; i++) {
					for (j = 0; j < numAttributes; j++) {
						data[i][j] = (data[i][j] - mean[j]) / stDev[j]
					}
			    }
			}

			// ************************
			// perform the clustering
			// ************************
			
			double[][] clusterCentres = new double[k][numAttributes]
			int[] classNum = new int[numRecords]
			boolean[] deadCluster = new boolean[k]
			int[] clusterCentreCounts
			
			// initialize the clusterCentres with randomly selected points

			for (m = 0; m < k; m++) {
				i = (int)(Math.random() * numRecords)
				for (j = 0; j < numAttributes; j++) {
					clusterCentres[m][j] = data[i][j]
				}
			}

			double[][] clusterCentreTotals
			int actualIterationNumber = 0
			XYSeries series = new XYSeries("")
			pluginHost.updateProgress("Clustering...", 0)
			oldProgress = -1
			for (int iter = 0; iter < numIterations; iter++) {
				clusterCentreTotals = new double[k][numAttributes]
				clusterCentreCounts = new int[k]
				int numChanged = 0
			
				for (i = 0; i < numRecords; i++) {
					minDist = Double.POSITIVE_INFINITY
					int startingClass = classNum[i]
					for (m = 0; m < k; m++) {
						if (!deadCluster[m]) {
							dist = 0
							for (j = 0; j < numAttributes; j++) {
								dist += (clusterCentres[m][j] - data[i][j]) * (clusterCentres[m][j] - data[i][j])
							}
							if (dist < minDist) {
								minDist = dist
								classNum[i] = m
							}
						}
					}
					clusterCentreCounts[classNum[i]]++
					for (j = 0; j < numAttributes; j++) {
						clusterCentreTotals[classNum[i]][j] += data[i][j]
					}
					if (startingClass != classNum[i]) numChanged++
			    }

				// update the cluster centres
			    for (m = 0; m < k; m++) {
					for (j = 0; j < numAttributes; j++) {
						if (clusterCentreCounts[m] > 0) {
							clusterCentres[m][j] = clusterCentreTotals[m][j] / clusterCentreCounts[m]
						} else {
							deadCluster[m] = true
						}
					}
				}

				series.add(iter, 100f * numChanged / numRecords)
				if (100f * numChanged / numRecords <= minChange) {
					actualIterationNumber = iter
					break
				}
				
				progress = (int)(100f * iter / (numIterations - 1))
				if (progress != oldProgress) {
					pluginHost.updateProgress("Clustering...", progress)
					oldProgress = progress
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}			    
			}

			// output the chart of between iteration change
			XYSeriesCollection xyCollection = new XYSeriesCollection();
        	xyCollection.addSeries(series)

        	JFreeChart chart = ChartFactory.createXYLineChart(
            "Cluster Asignment Change Between Iterations",
            "Iteration Num.",
            "Percent Changed (% of records)",
            xyCollection,
            PlotOrientation.VERTICAL,  // Plot Orientation
            false,                      // Show Legend
            true,                      // Use tooltips
            false                      // Configure chart to generate URLs?
            );

			XYPlot plot = (XYPlot)chart.getPlot();
			plot.setBackgroundPaint(Color.WHITE); 
			plot.setDomainGridlinePaint(Color.BLACK); 
			plot.setRangeGridlinePaint(Color.BLACK); 

            ChartPanel chartPanel = new ChartPanel(chart)
			chartPanel.setPreferredSize(new Dimension(700, 300))
			
        	pluginHost.returnData(chartPanel)


			// figure out the mapping of input clusters to output cluster 
			// numbers allowing for dead clusters.
			int[] outClassNum = new int[numRecords]
			int ko = 1
			for (m = 0; m < k; m++) {
				if (!deadCluster[m]) {
					outClassNum[m] = ko
					ko++
				}
			}
			ko--


			// output the attribute table
			DBFField field = new DBFField();
            field = new DBFField();
            field.setName("CLUSTER");
            field.setDataType(DBFField.DBFDataType.NUMERIC);
            field.setFieldLength(10);
            field.setDecimalCount(0);
            table.addField(field);

			for (i = 0; i < numRecords; i++) {
		    	recData = table.getRecord(i);
                recData[recData.length - 1] = new Double(outClassNum[classNum[i]]);
                table.updateRecord(i, recData);
		    }

			pluginHost.returnData(table.getFileName())

			// produce the cluster report
			DecimalFormat df = new DecimalFormat("###,###,###,##0.00")
			
			StringBuilder ret = new StringBuilder()
			ret.append("<!DOCTYPE html>")
			ret.append('<html lang="en">')
			
			ret.append("<head>")
			ret.append("<title>Cluster Analysis Summary</title>").append("\n")
			
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
			ret.append("<body><h1><i>k</i>-Means Cluster Analysis Summary</h1>").append("\n")

			ret.append("<p>Input file: ${inputFile}</p>")
			ret.append("<p>Rescaling Method: ${rescaleMode}</p>")
			ret.append("<p>Number of Iterations: ${actualIterationNumber}</p>")
			ret.append("<p>Number of Starting Cluster: ${k}</p>")
			ret.append("<p>Number of Final Cluster: ${ko}</p>")
			
			ret.append("<br><table border=\"1\" cellspacing=\"0\" cellpadding=\"3\">").append("\n")
			ret.append("<caption>Cluster Centers</caption>")
			ret.append("<tr><th>Attribute</th>")
			for (m = 0; m < k; m++) {
				if (!deadCluster[m]) {
					ret.append("<th>C. ${outClassNum[m]}</th>")
				}
			}

			if (rescaleMode.equals("normalize")) {
				for (j = 0; j < numAttributes; j++) {
					ret.append("<tr><td>${inputData[j + 1]}</td>")
					for (m = 0; m < k; m++) {
						if (!deadCluster[m]) {
							ret.append("<td class=\"numberCell\">")
							val = clusterCentres[m][j] * (maxVals[j] - minVals[j]) + minVals[j]
							ret.append(df.format(val))
							ret.append("</td>")
						}
					}
					ret.append("</tr>")
				}
			} else if (rescaleMode.equals("standardize")) {
				for (j = 0; j < numAttributes; j++) {
					ret.append("<tr><td>${inputData[j + 1]}</td>")
					for (m = 0; m < k; m++) {
						if (!deadCluster[m]) {
							ret.append("<td class=\"numberCell\">")
							val = clusterCentres[m][j] * stDev[j] + mean[j]
							ret.append(df.format(val))
							ret.append("</td>")
						}
					}
					ret.append("</tr>")
				}
			} else { // raw data
				for (j = 0; j < numAttributes; j++) {
					ret.append("<tr><td>${inputData[j + 1]}</td>")
					for (m = 0; m < k; m++) {
						if (!deadCluster[m]) {
							ret.append("<td class=\"numberCell\">")
							ret.append(df.format(clusterCentres[m][j]))
							ret.append("</td>")
						}
					}
					ret.append("</tr>")
				}
			}

			ret.append("<tr><td>Percent of Recs.</td>")
			for (m = 0; m < k; m++) {
				if (!deadCluster[m]) {
					ret.append("<td class=\"numberCell\">")
					ret.append(df.format(100f * clusterCentreCounts[m] / numRecords))
					ret.append("</td>")
				}
			}
			ret.append("</tr>")
			ret.append("</table>")

			
			// Create a cluster centre distance table
			ret.append("<p><table border=\"1\" cellspacing=\"0\" cellpadding=\"3\">").append("\n")
			if (rescaleMode.equals("normalize")) {
				ret.append("<caption>Cluster Center Distances (Normalized)</caption>")
			} else if (rescaleMode.contains("standard")) {
				ret.append("<caption>Cluster Center Distances (Standardized)</caption>")
			} else {
				ret.append("<caption>Cluster Center Distances</caption>")
			}
			
			ret.append("<th></th>")
			for (m = 0; m < k; m++) {
				if (!deadCluster[m]) {
					ret.append("<th>C. ${outClassNum[m]}</th>")
				}
			}
			ret.append("\n")
			
			for (m = 0; m < k; m++) {
				if (!deadCluster[m]) {
					ret.append("<tr><th>C. ${outClassNum[m]}</th>")
					for (int p = 0; p < k; p++) {
						if (!deadCluster[p]) {
							if (p > m) {
								ret.append("<td class=\"numberCell\">")
								dist = 0
								for (j = 0; j < numAttributes; j++) {
									dist += (clusterCentres[m][j] - clusterCentres[p][j]) * (clusterCentres[m][j] - clusterCentres[p][j])
								}
								ret.append(df.format(Math.sqrt(dist)))
								ret.append("</td>")
							} else if (p == m) {
								ret.append("<td class=\"numberCell\">0.00</td>")
							} else {
								ret.append("<td></td>")
							}
						}
					}
					ret.append("</tr>")
				}
			}
			ret.append("</table></p>")


			// calculate the cluster variance
			double[] clusterVariance = new double[k]
			double meanVariance
			for (m = 0; m < k; m++) {
				if (!deadCluster[m]) {
					for (i = 0; i < numRecords; i++) {
						if (classNum[i] == m) {
							for (j = 0; j < numAttributes; j++) {
								clusterVariance[m] += (clusterCentres[m][j] - data[i][j]) * (clusterCentres[m][j] - data[i][j])
							}
						}
				    }
				    clusterVariance[m] = Math.sqrt(clusterVariance[m] / clusterCentreCounts[m])
				    meanVariance += clusterVariance[m]
				}
			}

			ret.append("<p><table border=\"1\" cellspacing=\"0\" cellpadding=\"3\">").append("\n")
			ret.append("<caption>Cluster Standard Deviation</caption>")
			ret.append("<th>Cluster</th><th>Std. Dev.</th>\n")
			for (m = 0; m < k; m++) {
				if (!deadCluster[m]) {
					ret.append("<tr><td class=\"numberCell\">${outClassNum[m]}</td>")
					ret.append("<td class=\"numberCell\">${df.format(clusterVariance[m])}</td></tr>")
				}
			}
			ret.append("</table></p>")
			ret.append("<p>Mean Cluster Standard Deviation: ${df.format(meanVariance / ko)}</p>")


			

			ret.append("</body></html>")
			pluginHost.returnData(ret.toString());
			
        } catch (OutOfMemoryError oe) {
            pluginHost.showFeedback("An out-of-memory error has occurred during operation.")
	    } catch (Exception e) {
	        pluginHost.showFeedback("An error has occurred during operation. See log file for details.")
	        pluginHost.logException("Error in " + descriptiveName, e)
        } finally {
        	// reset the progress bar
        	pluginHost.updateProgress("", 0)
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
    def f = new ClusterAttributes(pluginHost, args, name, descriptiveName)
}
