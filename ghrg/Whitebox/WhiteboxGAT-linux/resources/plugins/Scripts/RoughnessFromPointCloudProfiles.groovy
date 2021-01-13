/*
 * Copyright (C) 2016 Dr. John Lindsay <jlindsay@uoguelph.ca> and
 * Melanie Chabot.
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
import java.util.Date
import java.beans.PropertyChangeEvent
import java.beans.PropertyChangeListener
import java.awt.Dimension
import java.awt.Color;
import javax.swing.JPanel;
import java.text.DecimalFormat
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.ui.plugin_dialog.*
import groovy.transform.CompileStatic
import whitebox.structures.BoundingBox
import whitebox.geospatialfiles.LASReader
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.ShapeType
import whitebox.geospatialfiles.LASReader.PointRecord
import whitebox.structures.KdTree
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.attributes.*
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.StandardChartTheme;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.plot.PolarPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.renderer.DefaultPolarItemRenderer;
import org.jfree.data.time.Month;
import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;
import org.jfree.data.xy.XYDataset;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RectangleInsets;
import org.jfree.ui.RefineryUtilities;
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
import org.apache.commons.math3.linear.*;



// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "RoughnessFromPointCloudProfiles"
def descriptiveName = "Roughness From Point Cloud Profiles (RFPCP)"
def description = "Extracts profile roughness from LiDAR point cloud data."
def toolboxes = ["LidarTools"]

public class RoughnessFromPointCloudProfiles implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private String descriptiveName
	private ScriptDialog sd;
	
	public RoughnessFromPointCloudProfiles(WhiteboxPluginHost pluginHost, 
		String[] args, String name, String descriptiveName) {
		this.pluginHost = pluginHost;
		this.descriptiveName = descriptiveName;
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
			sd.addDialogFile("Input LAS file", "Input LAS File:", "open", "LAS Files (*.las), LAS", true, false)
			sd.addDialogFile("Input Vector Polygon for clipping point cloud (optional)", "Input Vector Clipping Polygon File (Optional):", "open", "Vector Files (*.shp), SHP", true, true)
			sd.addDialogDataInput("Enter the spacing between profiles", "Between-Profile Spacing (m):", "0.50", true, false)
            sd.addDialogDataInput("Enter spacing of points along the profile", "Profile Sampling Density (m; Along-Profile Spacing):", "0.25", true, false)
            //sd.addDialogDataInput("Moving average detrend filter length (metres)", "Moving Average Detrend Filter Length (m):", "", true, false)
            sd.addDialogDataInput("De-trending polynomial order (1-10)", "De-trending polynomial order (1 to 10):", "3", true, false)
            DialogComboBox comboBox = sd.addDialogComboBox("Enter Roughness Statistic ", "Roughness Statistic:", ["RMS Height", "Correlation Length", "Shadow Roughness Index"], 0)
            DialogDataInput incidenceBox = sd.addDialogDataInput("Enter Satellite Incidence Angle", "Incidence Angle (degrees):", "30", true, true)
            incidenceBox.visible = false
			DialogDataInput wavelengthBox = sd.addDialogDataInput("Enter Satellite Wavelength", "Wavelength (cm):", "5", true, true)
            wavelengthBox.visible = false
			DialogCheckBox modeCB = sd.addDialogCheckBox("Would you like to run in single azimuth mode or multi azimuth mode(0 to 360 degrees)?", "Run in Multi Azimuth Mode?", true)
			DialogDataInput azimuth = sd.addDialogDataInput("Enter Profile Azimuth, in degrees 0-360", "Profile Azimuth (degrees):", "90.0", true, true)
            azimuth.visible = false
            
            def lstr = { evt -> if (evt.getPropertyName().equals("value")) { 
            		String value = comboBox.getValue()
            		if (!value.isEmpty() && value != null) { 
            			if (value.toLowerCase().contains("shadow")) {
            				incidenceBox.visible = true
            				wavelengthBox.visible = true
		            	} else {
		            		incidenceBox.visible = false
		            		wavelengthBox.visible = false
		            	}
            		} 
            	} 
            } as PropertyChangeListener
            comboBox.addPropertyChangeListener(lstr)
           
            def lstr2 = { evt -> if (evt.getPropertyName().equals("value")) { 
            		String value = modeCB.getValue()
            		if (!value.isEmpty() && value != null) { 
            			if (value.toLowerCase().contains("t")) {
            				azimuth.visible = false
            				
		            	} else {
		            		azimuth.visible = true
		            		
		            	}
            		} 
            	} 
            } as PropertyChangeListener
            modeCB.addPropertyChangeListener(lstr2)

            
			// resize the dialog to the standard size and display it
			sd.setSize(800, 400)
			sd.visible = true
		}
	}



	// The execute function is the main part of the tool, where the actual
	// work is completed.
	@CompileStatic
	private void execute(String[] args) {
		try {
			int progress, oldProgress
			int a 
	  		double x, y, z, rotatedX, rotatedY
	  		int filterSize = 10
	  		
			if (args.length != 10) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String inputFile = args[0]
			boolean isClipFileSpecified = false
			String clipFile 
			if (!args[1].toLowerCase().contains("not specified") && args[1].trim().length() > 0) {
				isClipFileSpecified = true
				clipFile = args[1]
			}
			double profileSpacing = Double.parseDouble(args[2])
			double samplingDensity = Double.parseDouble(args[3])
//			double detrendDist = Double.parseDouble(args[4]);
//			filterSize = (int)((detrendDist / samplingDensity) / 2.0);
//			if (filterSize < 1) {
//				filterSize = 1;
//			}
			int polyOrder = (int)Integer.parseInt(args[4]);
			if (polyOrder > 10) { polyOrder = 10 }
			if (polyOrder < 1) { polyOrder = 1 }
			String statType = "rms"
			if (args[5].toLowerCase().contains("cor")) {
				statType = "cor"
			} else if (args[5].toLowerCase().contains("shadow")) {
				statType = "shadow"
			}
			double incidenceAngle = 0
			if (!args[6].toLowerCase().contains("not specified")) {
				incidenceAngle = Math.toRadians(Double.parseDouble(args[6]))
			}
			double wavelength = 0
			if (!args[7].toLowerCase().contains("not specified")) {
				wavelength = (Double.parseDouble(args[7]) / 100)
			}
			boolean multiAzMode = true
			double profileAzimuth = 0
			if (args[8].toLowerCase().contains("f")) {
				multiAzMode = false;
				if (!args[8].toLowerCase().contains("not specified")) {
			
					profileAzimuth = -Math.toRadians(Double.parseDouble(args[9]))
				}
			}
			
			if (samplingDensity >= profileSpacing) { 
				pluginHost.showFeedback("Warning: The Profile Spacing should be larger than the Sampling Density")
			}

			pluginHost.updateProgress("Reading LAS Dataset:", 0)
			LASReader las = new LASReader(inputFile)
			int numPoints = (int)las.getNumPointRecords()
			double unrotatedMinX = las.getMinX()
			double unrotatedMaxX = las.getMaxX()
			double unrotatedMinY = las.getMinY()
			double unrotatedMaxY = las.getMaxY()
			BoundingBox extent = new BoundingBox(unrotatedMinX, unrotatedMinY, unrotatedMaxX, unrotatedMaxY);
			double midX = unrotatedMinX + (unrotatedMaxX - unrotatedMinX) / 2.0
			double midY = unrotatedMinY + (unrotatedMaxY - unrotatedMinY) / 2.0
			KdTree<Double> pointsTree = new KdTree.SqrEuclid<Double>(2, new Integer(numPoints))
			List<KdTree.Entry<Double>> results;
		
			double[] entry;
			double minXRotated = Double.POSITIVE_INFINITY
			double maxXRotated = Double.NEGATIVE_INFINITY
			double minYRotated = Double.POSITIVE_INFINITY
			double maxYRotated = Double.NEGATIVE_INFINITY
			ArrayList<PointRecord> points = las.getPointRecordsInBoundingBox(extent)
	 		if (!isClipFileSpecified) {
	 			oldProgress = 0
	 			a = 0
		 		for (PointRecord point : points) {
		 			a++
					entry = [point.getY(), point.getX()]
					pointsTree.addPoint(entry, point.getZ());
		 			progress = (int)(100f * a / numPoints)
					if (progress != oldProgress) {
						pluginHost.updateProgress("Reading LAS Dataset:", progress)
						oldProgress = progress
						if (pluginHost.isRequestForOperationCancelSet()) {
	                        pluginHost.showFeedback("Operation cancelled")
							return
	                    }
					}
		 		}
	 		} else {
	 			
				def clipPoly = new ShapeFile(clipFile)
				ShapeType shapeType = clipPoly.getShapeType()

				unrotatedMinX = Double.POSITIVE_INFINITY
				unrotatedMaxX = Double.NEGATIVE_INFINITY
				unrotatedMinY = Double.POSITIVE_INFINITY
				unrotatedMaxY = Double.NEGATIVE_INFINITY
	
				if (shapeType.getBaseType() != ShapeType.POLYGON) {
					pluginHost.showFeedback("The input shapefile must be of a POLYGON base ShapeType.")
	                return
				}
		
				double[][] nodes = clipPoly.getRecord(0).getGeometry().getPoints()
				int numNodes = nodes.length;
				int m;
			    int n;
			    boolean result = false;
				oldProgress = 0
				a = 0
	 			for (PointRecord point : points) {
	 				a++
		 			x = point.getX()
		 			y = point.getY()
		 			z = point.getZ()
		 			n = numNodes - 1
		 			result = false;
		 			for (m = 0; m < numNodes; n = m++) {
				    	if ((nodes[m][1] > y) != (nodes[n][1] > y) &&
				            (x < (nodes[n][0] - nodes[m][0]) * (y - nodes[m][1]) / (nodes[n][1]-nodes[m][1]) + nodes[m][0])) {
				          result = !result;
				        }
				    }
				    if (result) {
						entry = [y, x]
						pointsTree.addPoint(entry, z);
						if (x < unrotatedMinX) { unrotatedMinX = x }
						if (x > unrotatedMaxX) { unrotatedMaxX = x }
						if (y < unrotatedMinY) { unrotatedMinY = y }
						if (y > unrotatedMaxY) { unrotatedMaxY = y }
				    }
				    progress = (int)(100f * a / numPoints)
					if (progress != oldProgress) {
						pluginHost.updateProgress("Reading LAS Dataset:", progress)
						oldProgress = progress
						if (pluginHost.isRequestForOperationCancelSet()) {
	                        pluginHost.showFeedback("Operation cancelled")
							return
	                    }
					}
		 		}

		 		midX = unrotatedMinX + (unrotatedMaxX - unrotatedMinX) / 2.0
				midY = unrotatedMinY + (unrotatedMaxY - unrotatedMinY) / 2.0
				
				extent = new BoundingBox(unrotatedMinX, unrotatedMinY, unrotatedMaxX, unrotatedMaxY);
	 		}

			if (multiAzMode) {

				double interval = 7.5 // the interval, in degrees, at which the profile azimuth is sampled
				int numIntervals = (int)(360 / interval)
				double deg;
				
				double[] meanData = new double[numIntervals]
				double[] sDData = new double[numIntervals]
				double[] medianData = new double[numIntervals]
				double[] lowerQuartileData = new double[numIntervals]
				double[] upperQuartileData = new double[numIntervals]
				double[] interQuartileData = new double[numIntervals]

				for (int degInt = 0; degInt < numIntervals; degInt++) {
					deg = degInt * interval
					profileAzimuth = -Math.toRadians(deg)
					//Rotate the corners of the unrotated bounding box to find new bounding box
					minXRotated = Double.POSITIVE_INFINITY
					maxXRotated = Double.NEGATIVE_INFINITY
					minYRotated = Double.POSITIVE_INFINITY
					maxYRotated = Double.NEGATIVE_INFINITY

					//Northwest corner
					x = unrotatedMinX - midX;
	                y = unrotatedMaxY - midY;
	                rotatedX = midX + (x * Math.cos(profileAzimuth)) - (y * Math.sin(profileAzimuth));
	                rotatedY = midY + (x * Math.sin(profileAzimuth)) + (y * Math.cos(profileAzimuth));
					if (rotatedX < minXRotated) { minXRotated = rotatedX }
					if (rotatedX > maxXRotated) { maxXRotated = rotatedX }
					if (rotatedY < minYRotated) { minYRotated = rotatedY }
					if (rotatedY > maxYRotated) { maxYRotated = rotatedY }
					
			    	//Northeast corner
					x = unrotatedMaxX - midX;
	                y = unrotatedMaxY - midY;
	                rotatedX = midX + (x * Math.cos(profileAzimuth)) - (y * Math.sin(profileAzimuth));
	                rotatedY = midY + (x * Math.sin(profileAzimuth)) + (y * Math.cos(profileAzimuth));
					if (rotatedX < minXRotated) { minXRotated = rotatedX }
					if (rotatedX > maxXRotated) { maxXRotated = rotatedX }
					if (rotatedY < minYRotated) { minYRotated = rotatedY }
					if (rotatedY > maxYRotated) { maxYRotated = rotatedY }

			    	//Southwest corner
					x = unrotatedMinX - midX;
	                y = unrotatedMinY - midY;
	                rotatedX = midX + (x * Math.cos(profileAzimuth)) - (y * Math.sin(profileAzimuth));
	                rotatedY = midY + (x * Math.sin(profileAzimuth)) + (y * Math.cos(profileAzimuth));
					if (rotatedX < minXRotated) { minXRotated = rotatedX }
					if (rotatedX > maxXRotated) { maxXRotated = rotatedX }
					if (rotatedY < minYRotated) { minYRotated = rotatedY }
					if (rotatedY > maxYRotated) { maxYRotated = rotatedY }
			    
					//Southeast corner
					x = unrotatedMaxX - midX;
	                y = unrotatedMinY - midY;
	                rotatedX = midX + (x * Math.cos(profileAzimuth)) - (y * Math.sin(profileAzimuth));
	                rotatedY = midY + (x * Math.sin(profileAzimuth)) + (y * Math.cos(profileAzimuth));
					if (rotatedX < minXRotated) { minXRotated = rotatedX }
					if (rotatedX > maxXRotated) { maxXRotated = rotatedX }
					if (rotatedY < minYRotated) { minYRotated = rotatedY }
					if (rotatedY > maxYRotated) { maxYRotated = rotatedY }
					
			    	double deltaXRotated = maxXRotated - minXRotated;
					double deltaYRotated = maxYRotated - minYRotated;
					int numProfiles = (int) (deltaXRotated / profileSpacing) + 1;
					int numSamples = (int) (deltaYRotated / samplingDensity) + 1;
					double dist
					double noData = -32768.0
					double[] statValues = new double[numProfiles]

					int[] profLength = new int[numProfiles];
					oldProgress = -1
					for (int column = 0; column < numProfiles; column++) {
						int i = 0
						for (int row = 0; row < numSamples; row++) {
							x = minXRotated + column * profileSpacing - midX;
							y = minYRotated + row * samplingDensity - midY;
							rotatedX = midX + (x * Math.cos(profileAzimuth)) - (y * Math.sin(profileAzimuth));
			                rotatedY = midY + (x * Math.sin(profileAzimuth)) + (y * Math.cos(profileAzimuth));
							if (extent.isPointInBox(rotatedX, rotatedY)) {
								i++
							}
						}
						profLength[column] = i
						progress = (int)(100f * column / (numProfiles - 1))
						if (progress != oldProgress) {
							pluginHost.updateProgress("Calculating Profile Lengths (loop ${degInt + 1} of $numIntervals):", progress)
							oldProgress = progress
							if (pluginHost.isRequestForOperationCancelSet()) {
		                        pluginHost.showFeedback("Operation cancelled")
								return
		                    }
						}
					}

					oldProgress = -1;
					for (int column = 0; column < numProfiles; column++) {
						if (profLength[column] > 2) {
							int i = 0;
							double[] zValues = new double[profLength[column]] //numSamples]
							for (int row = 0; row < numSamples; row++) {
								x = minXRotated + column * profileSpacing - midX;
								y = minYRotated + row * samplingDensity - midY;
								rotatedX = midX + (x * Math.cos(profileAzimuth)) - (y * Math.sin(profileAzimuth));
				                rotatedY = midY + (x * Math.sin(profileAzimuth)) + (y * Math.cos(profileAzimuth));
				                if (extent.isPointInBox(rotatedX, rotatedY)) {
									entry = [rotatedY, rotatedX];
				                    results = pointsTree.nearestNeighbor(entry, 1, false)
				                    dist = results.get(0).distance
				                    if (samplingDensity > dist) {
					                    zValues[i] = (double) results.get(0).value
				                    } else {
				                    	zValues[i] = noData
				                    }
				                    i++;
				                }
							}
			
							switch (statType) {
								case "rms":
									//calculate RMS for profile
									//double[] detrendedProfile = detrendProfile(zValues, filterSize)
									double[] detrendedProfile = detrendProfilePolynomial(zValues, samplingDensity, polyOrder);
									statValues[column] = getRmsh(detrendedProfile, samplingDensity)
									break;
			
								case "cor":
									//calculate CL for profile
									//double[] detrendedProfile = detrendProfile(zValues, filterSize)
									double[] detrendedProfile = detrendProfilePolynomial(zValues, samplingDensity, polyOrder);
									statValues[column] = getCL(detrendedProfile, samplingDensity)
									break;
			
								case "shadow":
									// calculate the proportional shadow based rougness for profile
									statValues[column] = getSri(zValues, samplingDensity, incidenceAngle, wavelength)
							}
						}
						progress = (int)(100f * column / (numProfiles - 1))
						if (progress != oldProgress) {
							pluginHost.updateProgress("Calculating Statistic (loop ${degInt + 1} of $numIntervals):", progress)
							oldProgress = progress
							if (pluginHost.isRequestForOperationCancelSet()) {
		                        pluginHost.showFeedback("Operation cancelled")
								return
		                    }
						}
					}
				
			    	// calculate the overall stats
			    	double overallMean, overallSD, median, lowerQuartile, upperQuartile, interquartileRange;
			    	double sumX = 0;
			    	double sumXX = 0;
			    	int N = 0;
			    	ArrayList<Double> statValues2 = new ArrayList<Double>(numProfiles)
			    	for (int row = 0; row < numProfiles; row++) {
			            x = statValues[row];
			            if (x != noData) {
//			            	if (Double.isNaN(x)) {
//			            		pluginHost.showFeedback("$overallMean, $sumX, $N")
//			            	}
			            	statValues2.add(x)
			            	sumX += x;
			            	sumXX += x * x;
			            	N++;
			            }
			    	}
			    	
			    	overallMean = sumX / N;
			    	overallSD = Math.sqrt(sumXX / N - overallMean * overallMean);

					Collections.sort(statValues2);
			    	int len = statValues2.size();
					if (len % 2 == 0) {
						median = (statValues2.get((int)(len / 2)) + statValues2.get((int)(len/2 - 1)))/2;
					} else {
						median = statValues2.get((int)(len / 2));
					}
					lowerQuartile = statValues2.get((int)(Math.round((double)(len * 0.25))))
					upperQuartile = statValues2.get((int)Math.round((double)(len * 0.75)))
					interquartileRange = upperQuartile - lowerQuartile

					meanData[degInt] = overallMean
					sDData[degInt] = overallSD
					medianData[degInt] = median
					lowerQuartileData[degInt] = lowerQuartile
					upperQuartileData[degInt] = upperQuartile
					interQuartileData[degInt] = interquartileRange
				}

				//Output the table and polar chart
				String chartTitle
		        switch (statType) {
					case "rms":
						chartTitle = "Profile Root-Mean-Square Height"
						break
		
					case "cor":
						chartTitle = "Profile Roughness Correlation Length"
						break

					default:
						chartTitle = "Profile Shadow Roughness Index"
				}
				XYSeriesCollection xyCollection = new XYSeriesCollection();
				XYSeries series1 = new XYSeries("Lower Quartile")
				XYSeries series2 = new XYSeries("Median")
				XYSeries series3 = new XYSeries("Mean")
				XYSeries series4 = new XYSeries("Upper Quartile")
			
				DecimalFormat df = new DecimalFormat("###,###,###,##0.0000")
				DecimalFormat df2 = new DecimalFormat("###,###,###,##0.0")
				StringBuilder ret = new StringBuilder()
				ret.append("<!DOCTYPE html>")
				ret.append('<html lang="en">')
	
				ret.append("<head>")
	            ret.append("<title>Roughness Statistics</title>").append("\n")
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
	            ret.append("<body>").append("\n")
				ret.append("<br><table border=\"1\" cellspacing=\"0\" cellpadding=\"3\">").append("\n")
				ret.append("<caption>$chartTitle</caption>")
				
				ret.append("<tr><th>Az</th><th>Mean</th></th><th>SD</th><th>25%</th><th>50%</th><th>75%</th><th>IQR</th></tr>")
				for (int degInt = 0; degInt < numIntervals; degInt++) {
					deg = degInt * interval
					ret.append("<tr>")
					ret.append("<td>${df2.format(deg)}</td>")
					ret.append("<td>${df.format(meanData[degInt])}</td>")
					ret.append("<td>${df.format(sDData[degInt])}</td>")
					ret.append("<td>${df.format(lowerQuartileData[degInt])}</td>")
					ret.append("<td>${df.format(medianData[degInt])}</td>")
					ret.append("<td>${df.format(upperQuartileData[degInt])}</td>")
					ret.append("<td>${df.format(interQuartileData[degInt])}</td>")
					
					ret.append("</tr>")

					series1.add(deg, lowerQuartileData[degInt])
					series2.add(deg, medianData[degInt])
					series3.add(deg, meanData[degInt])
					series4.add(deg, upperQuartileData[degInt])
					
				}
				ret.append("</table>")
				ret.append("</body></html>")
				pluginHost.returnData(ret.toString())

				xyCollection.addSeries(series1);
				xyCollection.addSeries(series2);
				xyCollection.addSeries(series3);
				xyCollection.addSeries(series4);

				JFreeChart chart = ChartFactory.createPolarChart(
					chartTitle, xyCollection, true, true, false
		        );
		
				PolarPlot plot = (PolarPlot)chart.getPlot();
				plot.setBackgroundPaint(Color.white);
			    plot.setRadiusGridlinePaint(Color.gray);
			    DefaultPolarItemRenderer renderer = (DefaultPolarItemRenderer) plot.getRenderer()
			    renderer.setShapesVisible(false);
			    plot.setRadiusGridlinesVisible(true);
			    plot.setAngleGridlinesVisible(true);
			    plot.setAngleLabelsVisible(true);
			    plot.setOutlineVisible(false);
			    plot.setAngleGridlinePaint(Color.BLACK);
			    plot.setRenderer(renderer);
		
				ChartPanel chartPanel = new ChartPanel(chart)
				chartPanel.setPreferredSize(new Dimension(450, 450))
				
		    	pluginHost.returnData(chartPanel)
	
				
			} else {

				minXRotated = Double.POSITIVE_INFINITY
				maxXRotated = Double.NEGATIVE_INFINITY
				minYRotated = Double.POSITIVE_INFINITY
				maxYRotated = Double.NEGATIVE_INFINITY

				//Northwest corner
				x = unrotatedMinX  - midX;
                y = unrotatedMaxY - midY;
                rotatedX = midX + (x * Math.cos(profileAzimuth)) - (y * Math.sin(profileAzimuth));
                rotatedY = midY + (x * Math.sin(profileAzimuth)) + (y * Math.cos(profileAzimuth));
				if (rotatedX < minXRotated) { minXRotated = rotatedX }
				if (rotatedX > maxXRotated) { maxXRotated = rotatedX }
				if (rotatedY < minYRotated) { minYRotated = rotatedY }
				if (rotatedY > maxYRotated) { maxYRotated = rotatedY }
				
		    	//Northeast corner
				x = unrotatedMaxX - midX;
                y = unrotatedMaxY - midY;
                rotatedX = midX + (x * Math.cos(profileAzimuth)) - (y * Math.sin(profileAzimuth));
                rotatedY = midY + (x * Math.sin(profileAzimuth)) + (y * Math.cos(profileAzimuth));
				if (rotatedX < minXRotated) { minXRotated = rotatedX }
				if (rotatedX > maxXRotated) { maxXRotated = rotatedX }
				if (rotatedY < minYRotated) { minYRotated = rotatedY }
				if (rotatedY > maxYRotated) { maxYRotated = rotatedY }

		    	//Southwest corner
				x = unrotatedMinX - midX;
                y = unrotatedMinY - midY;
                rotatedX = midX + (x * Math.cos(profileAzimuth)) - (y * Math.sin(profileAzimuth));
                rotatedY = midY + (x * Math.sin(profileAzimuth)) + (y * Math.cos(profileAzimuth));
				if (rotatedX < minXRotated) { minXRotated = rotatedX }
				if (rotatedX > maxXRotated) { maxXRotated = rotatedX }
				if (rotatedY < minYRotated) { minYRotated = rotatedY }
				if (rotatedY > maxYRotated) { maxYRotated = rotatedY }
		    
				//Southeast corner
				x = unrotatedMaxX - midX;
                y = unrotatedMinY - midY;
                rotatedX = midX + (x * Math.cos(profileAzimuth)) - (y * Math.sin(profileAzimuth));
                rotatedY = midY + (x * Math.sin(profileAzimuth)) + (y * Math.cos(profileAzimuth));
				if (rotatedX < minXRotated) { minXRotated = rotatedX }
				if (rotatedX > maxXRotated) { maxXRotated = rotatedX }
				if (rotatedY < minYRotated) { minYRotated = rotatedY }
				if (rotatedY > maxYRotated) { maxYRotated = rotatedY }

				double deltaXRotated = maxXRotated - minXRotated
				double deltaYRotated = maxYRotated - minYRotated
				int numProfiles = (int) (deltaXRotated / profileSpacing) + 1
				int numSamples = (int) (deltaYRotated / samplingDensity) + 1
				double dist
				double noData = -32768.0
				double[] statValues = new double[numProfiles]
				oldProgress = 0

				// Create the output file
				DBFField[] fields = new DBFField[3];
	
	            fields[0] = new DBFField();
	            fields[0].setName("ID");
	            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
	            fields[0].setFieldLength(10);
	            fields[0].setDecimalCount(0);

	            fields[1] = new DBFField();
	            fields[1].setName("PROFILE");
	            fields[1].setDataType(DBFField.DBFDataType.NUMERIC);
	            fields[1].setFieldLength(10);
	            fields[1].setDecimalCount(0);

	            fields[2] = new DBFField();
	            fields[2].setName("PROF_PT");
	            fields[2].setDataType(DBFField.DBFDataType.NUMERIC);
	            fields[2].setFieldLength(10);
	            fields[2].setDecimalCount(0);
	            
	            String outputFile = pluginHost.getWorkingDirectory() + "RoughProf.shp"
				ShapeFile output = new ShapeFile(outputFile, ShapeType.POINT, fields);
				int j = 1
				whitebox.geospatialfiles.shapefile.Point wbGeometry
				int[] profLength = new int[numProfiles];
				oldProgress = -1;
				for (int column = 0; column < numProfiles; column++) {
					int i = 1
					for (int row = 0; row < numSamples; row++) {
						x = minXRotated + column * profileSpacing - midX;
						y = minYRotated + row * samplingDensity - midY;
						rotatedX = midX + (x * Math.cos(profileAzimuth)) - (y * Math.sin(profileAzimuth));
		                rotatedY = midY + (x * Math.sin(profileAzimuth)) + (y * Math.cos(profileAzimuth));
						if (extent.isPointInBox(rotatedX, rotatedY)) {
							wbGeometry = new whitebox.geospatialfiles.shapefile.Point(rotatedX, rotatedY);                  
	                		Object[] rowData = new Object[3]
	                		rowData[0] = new Double(j)
	                		j++
	                		rowData[1] = new Double(column)
	                		rowData[2] = new Double(i)
	                		i++
	                		output.addRecord(wbGeometry, rowData);
						}
					}
					profLength[column] = i - 1
					progress = (int)(100f * column / (numProfiles - 1))
					if (progress != oldProgress) {
						pluginHost.updateProgress("Calculating Profile Lengths:", progress)
						oldProgress = progress
						if (pluginHost.isRequestForOperationCancelSet()) {
	                        pluginHost.showFeedback("Operation cancelled")
							return
	                    }
					}
				}

				output.write()
				pluginHost.returnData(outputFile);
				
				oldProgress = -1;
				for (int column = 0; column < numProfiles; column++) {
					if (profLength[column] > 2) {
						int i = 0;
						double[] zValues = new double[profLength[column]] //numSamples]
						for (int row = 0; row < numSamples; row++) {
							x = minXRotated + column * profileSpacing - midX;
							y = minYRotated + row * samplingDensity - midY;
							rotatedX = midX + (x * Math.cos(profileAzimuth)) - (y * Math.sin(profileAzimuth));
			                rotatedY = midY + (x * Math.sin(profileAzimuth)) + (y * Math.cos(profileAzimuth));
			                if (extent.isPointInBox(rotatedX, rotatedY)) {
								entry = [rotatedY, rotatedX];
			                    results = pointsTree.nearestNeighbor(entry, 1, false)
			                    dist = results.get(0).distance
			                    if (samplingDensity > dist) {
				                    zValues[i] = (double) results.get(0).value
			                    } else {
			                    	zValues[i] = noData
			                    }
			                    i++;
			                }
						}
		
						switch (statType) {
							case "rms":
								//calculate RMS for profile
								//double[] detrendedProfile = detrendProfile(zValues, filterSize);
								double[] detrendedProfile = detrendProfilePolynomial(zValues, samplingDensity, polyOrder);
								statValues[column] = getRmsh(detrendedProfile, samplingDensity);
								break;
		
							case "cor":
								//calculate CL for profile
								//double[] detrendedProfile = detrendProfile(zValues, filterSize);
								double[] detrendedProfile = detrendProfilePolynomial(zValues, samplingDensity, polyOrder);
								statValues[column] = getCL(detrendedProfile, samplingDensity);
								break;
		
							case "shadow":
								// calculate the proportional shadow based rougness for profile
								statValues[column] = getSri(zValues, samplingDensity, incidenceAngle, wavelength)
								break;
						}
					}
					progress = (int)(100f * column / (numProfiles - 1))
					if (progress != oldProgress) {
						pluginHost.updateProgress("Calculating Statistic:", progress)
						oldProgress = progress
						if (pluginHost.isRequestForOperationCancelSet()) {
	                        pluginHost.showFeedback("Operation cancelled")
							return
	                    }
					}
				}
				
				
	
			XYSeriesCollection xyCollection = new XYSeriesCollection();
			XYSeries series = new XYSeries("")
			//String outStr = ""
			for (int row = 0; row < numProfiles; row++) {
	            x = row;
	            y = statValues[row];
	            //outStr += statValues[a]+ "\t" + x + "\t" + y + "\n"
	            if (y != noData) {
	            	series.add(x, y)	
	            } else {
//	            	if (series.getItemCount() > 0) {
//		            	xyCollection.addSeries(series);	
//		            	series = new XYSeries("")
//	            	}	
	            }
	        }
	        xyCollection.addSeries(series);
	        //pluginHost.returnData(outStr)
	
	       String chartTitle
	       String rangeLabel
	       String domainLabel = "Profile Number"
			switch (statType) {
				case "rms":
					chartTitle = "Root-Mean-Square Height"
					rangeLabel = "RMS Height (Z units)"
					break
	
				case "cor":
					chartTitle = "Roughness Correlation Length"
					rangeLabel = "Correlaion Length (Map units)"
					break
					
				case "shadow":
					chartTitle = "Shadow Roughness Index"
					rangeLabel = "Proportion of Profile Length"
					break
					
			}
	
	       
			JFreeChart chart = ChartFactory.createXYLineChart(
	        chartTitle,
	        domainLabel,
	        rangeLabel,
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
	
			// TODO: change the line colour of each series to red.
	
			ChartPanel chartPanel = new ChartPanel(chart)
			chartPanel.setPreferredSize(new Dimension(700, 300))
			
	    	pluginHost.returnData(chartPanel)
	
	    	// calculate the overall stats
	    	double overallMean, overallSD, median, interquartileRange;
	    	double sumX = 0;
	    	double sumXX = 0;
	    	int N = 0;
	    	ArrayList<Double> statValues2 = new ArrayList<Double>(numProfiles)
	    	for (int row = 0; row < numProfiles; row++) {
	            x = statValues[row];
	            if (x != noData) {
	            	statValues2.add(x)
	            	sumX += x;
	            	sumXX += x * x;
	            	N++;
	            }
	    	}
	    	
	    	overallMean = sumX / N;
	    	overallSD = Math.sqrt(sumXX / N - overallMean * overallMean);
	
			Collections.sort(statValues2);
	    	int len = statValues2.size();
			if (len % 2 == 0) {
				median = (statValues2.get((int)(len/2)) + statValues2.get((int)(len/2 - 1)))/2;
			} else {
				median = statValues2.get((int)(len/2));
			}
			double lowerQuartile = statValues2.get((int)(Math.round((double)(len * 0.25))))
			double upperQuartile = statValues2.get((int)Math.round((double)(len * 0.75)))
			 interquartileRange = upperQuartile - lowerQuartile
	        	
			}

			
					
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

	@CompileStatic
    public double[] detrendProfilePolynomial(double[] zValues, double samplingDensity, int polyOrder) {
		int N = zValues.length;
		int validN = 0;
		double x, y;
		double noData = -32768.0;
		double[] zValuesDetrended = new double[N];
		int i, j, n;
		double sumX, sumY, sumXY, sumXX, sumYY;
	
		//regressionData = new double[5];
	    for (i = 0; i < N; i++) {
			if (zValues[i] != noData) {
				validN++
			}
	    }
		double[] yVals = new double[validN];
		int numCoefficients = polyOrder + 1;
	
		// Solve the forward transformation equations
	    double[][] forwardCoefficientMatrix = new double[validN][numCoefficients];
	    n = 0;
	    for (i = 0; i < N; i++) {
	    	if (zValues[i] != noData) {
	    		x = i * samplingDensity;
	    		y = zValues[i];
	    		yVals[n] = y;
		        for (j = 0; j <= polyOrder; j++) {
	                forwardCoefficientMatrix[n][j] = Math.pow(x, j);
		        }
		        n++
	    	}
	    }
	
	    RealMatrix coefficients = new Array2DRowRealMatrix(forwardCoefficientMatrix, false);
	    DecompositionSolver solver = new SingularValueDecomposition(coefficients).getSolver();
	    
	    RealVector constants = new ArrayRealVector(yVals, false);
	    RealVector solution = solver.solve(constants);
	    double[] forwardRegressCoeffX = new double[numCoefficients];
	    //println("Polynomial co-efficients")
	    for (int a = 0; a < numCoefficients; a++) {
	        forwardRegressCoeffX[a] = solution.getEntry(a);
	        //println("c$a: ${forwardRegressCoeffX[a]}");
	    }
	    //println("")
	
	    //double[] residualsX = new double[n];
	    //double SSresidX = 0;
	    n = 0
	    for (i = 0; i < N; i++) {
	    	if (zValues[i] != noData) {
		        double yHat = 0.0;
		        for (j = 0; j < numCoefficients; j++) {
		            yHat += forwardCoefficientMatrix[n][j] * forwardRegressCoeffX[j];
		        }
		        zValuesDetrended[i] = zValues[i] - yHat;
		        n++
	    	} else {
	    		 zValuesDetrended[i] = noData
	    	}
	    }
		
		return zValuesDetrended;
	}

//    @CompileStatic
//    public double[] detrendProfile(double[] zValues, int filterSize) {
//    	int N = zValues.length;
//    	double noData = -32768.0;
//    	double[] zValuesDetrended = new double[N];
//    	int i, j, n;
//    	double sum;
//    	for (i = 0; i < N; i++) {
//    		if (zValues[i] != noData) {
//
//				sum = zValues[i];
//				n = 1
//				// get the filterSize number of points to the left
//				j = i - 1;
//				while (j >= 0 && n < filterSize + 1) {
//					if (zValues[j] != noData) {
//						sum += zValues[j];
//	    				n++;
//					}
//					j--;
//				}
//
//				// get the filterSize number of points to the right
//				j = i + 1;
//				while (j < N && n < filterSize * 2 + 1) {
//					if (zValues[j] != noData) {
//						sum += zValues[j];
//	    				n++;
//					}
//					j++;
//				}
//				
//	    		zValuesDetrended[i] = zValues[i] - sum / n;
//    		} else {
//    			zValuesDetrended[i] = noData;
//    		}
//    	}
//
//    	return zValuesDetrended;
//    }
    
	@CompileStatic
    public double getRmsh(double[] zValues, double samplingDensity) {
    	double result;
    	double sum, z, mean, sumZZ;
    	double noData = -32768.0;
    	int row;
    	long N = 0;
    	for (row = 0; row < zValues.length; row++) {
            z = zValues[row];
            if (z != noData) {
                sum += z;
                N++;
            }
        }
        if (N > 1) {
        	mean = sum / N;
        	sumZZ = 0;
        	for (row = 0; row < zValues.length; row++) {
    			if (zValues[row] != noData) {
    				sumZZ += (zValues[row] - mean) * (zValues[row] - mean)
    			}
    		}
    		result = Math.sqrt(sumZZ / (N - 1))
        } else {
        	result = noData
        }
//        double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0, sumYY = 0;
//        long N = 0;
//        double noData = -32768.0
//        double x, y
//        int row
//        double slope, intercept, yHat
//           
//		for (row = 0; row < zValues.length; row++) {
//            x = row * samplingDensity;
//            y = zValues[row];
//            if (y != noData) {
//                sumX += x;
//                sumY += y;
//                sumXY += x * y;
//                sumXX += x * x;
//                sumYY += y * y;
//                N++;
//            }
//        }
//        if (N > 1) {
//        	slope = (N * sumXY - (sumX * sumY)) / (N * sumXX - (sumX * sumX));
//        	intercept = (sumY - slope * sumX) / N;
//    		sumXX = 0
//    		for (row = 0; row < zValues.length; row++) {
//    			if (zValues[row] != noData) {
//    				yHat = (slope * row * samplingDensity + intercept)
//    				sumXX += (zValues[row] - yHat) * (zValues[row] - yHat)
//    			}
//    		}
//    		result = Math.sqrt(sumXX / (N - 1))
//        } else {
//        	result = noData
//        }
       return result
    }
    
    @CompileStatic
    public double getCL(double[] zValues, double samplingDensity) {
		double result
	    double mean, sumZ, variance;
	    long N = 0;
	    double noData = -32768.0
	    double x, y, z
	    int row
	    double slope, intercept, yHat
	
	    // Calculate the sum, N, and the mean.
		for (row = 0; row < zValues.length; row++) {
	        z = zValues[row];
	        if (z != noData) {
	            sumZ += z;
	            N++;
	        }
	    }
	    mean = sumZ / N
	
		// Calculate the variance.
	//    for (row = 0; row < zValues.length; row++) {
	//        z = zValues[row];
	//        if (z != noData) {
	//            variance += (z - mean) * (z - mean)
	//        }
	//        	
	//    }
	    //variance = variance / N
	
		// Calculate the autocorrelation function (acf).
	    double[] acf = new double[zValues.length]
	    double[] tmp = new double[zValues.length]
	    int[] numPairs = new int[zValues.length]
	    double z1, z2
	    int lag
	    
	    for (int a = 0; a < zValues.length - 1; a++) {
	    	z1 = zValues[a]
	    	 if (z1 != noData){
	    	 	for (int b = a + 1; b < zValues.length; b++) {
	    	 		z2 = zValues[b]
	    	 		 if (z2 != noData) {
	    	 		 	lag = b - a;
	    	 		 	acf[lag - 1] += (z1 - mean) * (z2 - mean)
	    	 		 	numPairs[lag - 1]++
	    	 		 	tmp[lag - 1] += (z1 - mean) * (z1 - mean)
	    	 		 }
	    	 	}
	    	 }
	    }
	
	    // Normalize the acf by the variance
	     for (row = 0; row < zValues.length; row++) {
	     	if (numPairs[row] > 0) {
	     		acf[row] = acf[row] / tmp[row] //variance
	     		println("${acf[row]} \t${numPairs[row]}")
	     	}
	     }
	
		// Find the lag distance at which acf drops to 1 / e
	     double threshold = 1.0 / Math.E
	     int lagCrit
	     for (row = 0; row < zValues.length; row++) {
	     	if (numPairs[row] > 0) {
	     		if (acf[row] <= threshold) {
	     			lagCrit = row;
	     			if (row > 0) {
	     				result = (lagCrit - (threshold - acf[lagCrit]) / (acf[lagCrit - 1] - acf[lagCrit])) * samplingDensity
	     				//println("$row of ${zValues.length}")
	     			} else {
	     				result = 0;
	     			}
	     			break;
	     			
	     		}
	     	}
	      }
	     
	   return result
	}
    
	@CompileStatic
    public double getSri(double[] zValues, double samplingDensity, double incidenceAngle, double wavelength) {
    	double noData = -32768.0
        double x, y, z, zN, xN
        int row, rowN
        double dist, deltaZ, angleN
		double[] horizonAngle = new double[zValues.length]
		double[] hALocation = new double[zValues.length]
		double[] elevation = new double[zValues.length]
		
        // Calculate horizon angle for each location
        for (row = 0; row < zValues.length; row++) {
        	horizonAngle[row] = -99;
            x = row * samplingDensity;
            z = zValues[row];
            if (z != noData) {
               for (rowN = row - 1; rowN >= 0; rowN--) {
           			if (rowN > 0) {
               	 		zN = zValues[rowN];
               	 		if (zN != noData) {
		               	 	xN = rowN * samplingDensity;
		               	 	dist = (row - rowN) * samplingDensity; //x - xN;
		               	 	deltaZ = zN - z;
		               	 	angleN = Math.tan(deltaZ/dist)
		               	 	if (angleN > horizonAngle[row]) {
		               	 		horizonAngle[row] = angleN;
		               	 		hALocation[row] = xN
		               	 		elevation[row] = zN
               	 			}
               	 		}

               		}
               }
            }
        }

		// Calculate roughness element height
		double numValidElement, numShadowElement
		
		double[] roughnessElementHeight = new double[zValues.length]
		int searchWindowSize = (int)((10 * (wavelength / samplingDensity)) / 2) 
			if (searchWindowSize < 1) {
				searchWindowSize = 1
			}
		for (row = 0; row < zValues.length; row++) {
			z = zValues[row];
			if (z != noData) {
				numValidElement++
				for (rowN = row - searchWindowSize; rowN <= row + searchWindowSize; rowN++) {
					if (rowN >= 0 && rowN < zValues.length) {
						if (zValues[rowN] != noData) {
							if (zValues[rowN] < z) {
								z = zValues[rowN]
							}
						}
					}
				}
				roughnessElementHeight[row] = zValues[row] - z
			}

		}
		// Find Profile locations within shadows
		 boolean[] inShadow = new boolean[zValues.length]
		 for (row = 0; row < zValues.length; row++) {
		 	if (horizonAngle[row] > incidenceAngle && zValues[row] != noData) {
		 		if (roughnessElementHeight[row] >= wavelength) {
		 			inShadow[row] = true
		 			numShadowElement++
		 		}
		 	}
		 }

		if (numValidElement > 2) {
			return numShadowElement / numValidElement
		} else {
			return noData
		}
    }
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def myTool = new RoughnessFromPointCloudProfiles(pluginHost, args, name, descriptiveName)
}
