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
import java.text.DecimalFormat
import java.io.File
import java.util.Date
import java.util.ArrayList
import java.util.PriorityQueue
import java.util.Arrays
import java.util.Collections
import java.beans.PropertyChangeEvent
import java.beans.PropertyChangeListener
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.*
import whitebox.utilities.FileUtilities;
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase
import whitebox.geospatialfiles.VectorLayerInfo
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.geospatialfiles.shapefile.ShapeFileRecord
import whitebox.utilities.Topology
import whitebox.structures.BoundingBox
import whitebox.structures.KdTree
import java.util.Arrays
import org.apache.commons.math3.distribution.NormalDistribution
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
def name = "NearestNeighbourAnalysis"
def descriptiveName = "Nearest-Neighbour Analysis"
def description = "Performs a nearest-neighbour analysis on a vector points coverage"
def toolboxes = ["VectorTools", "StatisticalTools"]

public class NearestNeighbourAnalysis implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
    private String basisFunctionType = ""
	
    public NearestNeighbourAnalysis(WhiteboxPluginHost pluginHost, 
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
        	sd.addDialogFile("Input file", "Input Vector Points File:", "open", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogComboBox("Strategy for handle edge-effects", "Edge-Effect Strategy:", ["Tile data", "Do nothing"], 0)
	
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
        	int progress, oldProgress, i
        	double x, y, z
        	double north, south, east, west
        	List<KdTree.Entry<Integer>> results
            ArrayList<Double> xList = new ArrayList<>()
            ArrayList<Double> yList = new ArrayList<>()
            
			if (args.length != 2) {
                pluginHost.showFeedback("Incorrect number of arguments given to tool.")
                return
            }
            
        	// read the input parameters
            String inputFile = args[0]
			String edgeEffectStrategy = args[1].toLowerCase()
            
			ShapeFile input = new ShapeFile(inputFile)
			ShapeType shapeType = input.getShapeType()
            if (shapeType.getBaseType() != ShapeType.POINT) {
            	pluginHost.showFeedback("The input file must be of a POINT base ShapeType.")
            	return
            }

            north = input.getyMax()
	        south = input.getyMin()
	        east = input.getxMax()
	        west = input.getxMin()
	        double xRange = Math.abs(east - west)
	        double yRange = Math.abs(north - south)

			int numFeatures = input.getNumberOfRecords()
			double[][] point
			for (ShapeFileRecord record : input.records) {
				point = record.getGeometry().getPoints()
				for (int p = 0; p < point.length; p++) {
					xList.add(point[p][0])
					yList.add(point[p][1])
				}
				i++
                progress = (int)(100f * i / numFeatures)
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Reading Points:", progress)
            		oldProgress = progress
            	}
            	// check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
			}
			
			int numPoints = xList.size()
			int numPointsBuffered = numPoints
			if (!edgeEffectStrategy.contains("nothing")) {
				numPointsBuffered = numPoints * 9
			}
			
			KdTree<Integer> pointsTree = new KdTree.SqrEuclid<Integer>(2, new Integer(numPointsBuffered));
			
			for (i = 0; i < numPoints; i++) {
				double[] entry = new double[2]
				entry[0] = xList.get(i)
				entry[1] = yList.get(i)
				pointsTree.addPoint(entry, i)
			}

			// Now buffer the area with points all around 
			// to compensate for edge effects.
			if (edgeEffectStrategy.contains("tile")) {
				for (i = 0; i < numPoints; i++) {
					double[] entry = new double[2]
					entry[0] = xList.get(i)
					entry[1] = yList.get(i) + yRange
					pointsTree.addPoint(entry, i)
				}
	
				for (i = 0; i < numPoints; i++) {
					double[] entry = new double[2]
					entry[0] = xList.get(i)
					entry[1] = yList.get(i) - yRange
					pointsTree.addPoint(entry, i)
				}
	
				for (i = 0; i < numPoints; i++) {
					double[] entry = new double[2]
					entry[0] = xList.get(i) + xRange
					entry[1] = yList.get(i)
					pointsTree.addPoint(entry, i)
				}
	
				for (i = 0; i < numPoints; i++) {
					double[] entry = new double[2]
					entry[0] = xList.get(i) - xRange
					entry[1] = yList.get(i)
					pointsTree.addPoint(entry, i)
				}
				
				for (i = 0; i < numPoints; i++) {
					double[] entry = new double[2]
					entry[0] = xList.get(i) - xRange
					entry[1] = yList.get(i) + yRange
					pointsTree.addPoint(entry, i)
				}

				for (i = 0; i < numPoints; i++) {
					double[] entry = new double[2]
					entry[0] = xList.get(i) + xRange
					entry[1] = yList.get(i) + yRange
					pointsTree.addPoint(entry, i)
				}

				for (i = 0; i < numPoints; i++) {
					double[] entry = new double[2]
					entry[0] = xList.get(i) - xRange
					entry[1] = yList.get(i) - yRange
					pointsTree.addPoint(entry, i)
				}

				for (i = 0; i < numPoints; i++) {
					double[] entry = new double[2]
					entry[0] = xList.get(i) + xRange
					entry[1] = yList.get(i) - yRange
					pointsTree.addPoint(entry, i)
				}
			}

			double totalDist = 0
            double dist1, dist2
			oldProgress = -1
			for (i = 0; i < numPoints; i++) {
				double[] entry = new double[2]
                entry[0] = xList.get(i)
                entry[1] = yList.get(i)
                
                results = pointsTree.nearestNeighbor(entry, 2, true)
				dist1 = Math.sqrt(results.get(0).distance)
				totalDist += dist1
				
				progress = (int)(100f * i / (numPoints - 1))
				if (progress > oldProgress) {
					pluginHost.updateProgress("Calculating Mean Nearest-Neighbour Distance:", progress)
					oldProgress = progress
				}
				// check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
			}

			// calculate the area of the minimum bounding box
			double psi = 0
			double boxCentreX, boxCentreY
			double midX = west + (east - west) / 2.0
            double midY = south + (north - south) / 2.0
			double DegreeToRad = Math.PI / 180.0
			double RadToDegree = 180.0 / Math.PI
			double[][] verticesRotated = new double[numPoints][2]
            double[] axes = new double[2]
            axes[0] = 9999999
            axes[1] = 9999999
            double[] newBoundingBox = new double[4]
            double newXAxis = 0
        	double newYAxis = 0
        	int incrementVal = 2
	  		double slope
	  		final double rightAngle = Math.toRadians(90)
	  		
			for (int m = 0; m <= 180; m++) {
                psi = -m * 0.5 * DegreeToRad; // rotation in clockwise direction
                // Rotate each edge cell in the array by m degrees.
                for (int n = 0; n < numPoints; n++) {
                    x = xList.get(n) - midX;
                    y = yList.get(n) - midY;
                    verticesRotated[n][0] = (x * Math.cos(psi)) - (y * Math.sin(psi));
                    verticesRotated[n][1] = (x * Math.sin(psi)) + (y * Math.cos(psi));
                }
                // calculate the minimum bounding box in this coordinate 
                // system and see if it is less
                newBoundingBox[0] = Double.MAX_VALUE; // west
                newBoundingBox[1] = Double.MIN_VALUE; // east
                newBoundingBox[2] = Double.MAX_VALUE; // north
                newBoundingBox[3] = Double.MIN_VALUE; // south
                for (int n = 0; n < numPoints; n++) {
                    x = verticesRotated[n][0];
                    y = verticesRotated[n][1];
                    if (x < newBoundingBox[0]) {
                        newBoundingBox[0] = x;
                    }
                    if (x > newBoundingBox[1]) {
                        newBoundingBox[1] = x;
                    }
                    if (y < newBoundingBox[2]) {
                        newBoundingBox[2] = y;
                    }
                    if (y > newBoundingBox[3]) {
                        newBoundingBox[3] = y;
                    }
                }
                newXAxis = newBoundingBox[1] - newBoundingBox[0];
                newYAxis = newBoundingBox[3] - newBoundingBox[2];

                if ((newXAxis * newYAxis) < (axes[0] * axes[1])) { // minimize the area of the bounding box.
                    axes[0] = newXAxis;
                    axes[1] = newYAxis;

                    if (axes[0] > axes[1]) {
                        slope = -psi;
                    } else {
                        slope = -(rightAngle + psi);
                    }
                    x = newBoundingBox[0] + newXAxis / 2;
                    y = newBoundingBox[2] + newYAxis / 2;
                    boxCentreX = midX + (x * Math.cos(-psi)) - (y * Math.sin(-psi));
                    boxCentreY = midY + (x * Math.sin(-psi)) + (y * Math.cos(-psi));
                }
            }
            double A = axes[0] * axes[1]
            double P = axes[0] * 2 + axes[1] * 2
            
			
			DecimalFormat df = new DecimalFormat("###,###,##0.000")
			StringBuilder ret = new StringBuilder()
			ret.append("<!DOCTYPE html>").append("\n")
			ret.append('<html lang="en">').append("\n")

			ret.append("<head>").append("\n")
//            ret.append("<meta content=\"text/html; charset=iso-8859-1\" http-equiv=\"content-type\">")
            ret.append("<title>Nearest-Neighbour Analysis Output</title>").append("\n")
            
            ret.append("<style>")
			ret.append("table {margin-left: 15px;} ")
			ret.append("h1 {font-size: 14pt; margin-left: 15px; margin-right: 15px; text-align: center; font-family: Helvetica, Verdana, Geneva, Arial, sans-serif;} ")
			ret.append("p {font-size: 12pt; font-family: Helvetica, Verdana, Geneva, Arial, sans-serif; margin-left: 15px; margin-right: 15px;} ")
			ret.append("table {font-size: 12pt; font-family: Helvetica, Verdana, Geneva, Arial, sans-serif;}")
			ret.append("table th {border-width: 1px; padding: 8px; border-style: solid; border-color: #666666; background-color: #dedede; }")
			ret.append("table td {border-width: 1px; padding: 8px; border-style: solid; border-color: #666666; background-color: #ffffff; }")
			ret.append(".numberCell { text-align: right; }") 
            ret.append("</style>")
            
            ret.append("</head>").append("\n")
            ret.append("<body><h1>Nearest-Neighbour Analysis</h1>").append("\n")
            
            ret.append("<p><b>Input File:</b> &nbsp " + new File(inputFile).getName()).append("</p>\n")
            ret.append("<br><table border=\"1\" cellspacing=\"0\" cellpadding=\"3\">").append("\n")

            ret.append("<tr><th>Statistic</th><th>Value</th></tr>")

            double meanNNDist = totalDist / numPoints
			ret.append("<tr><td>Observed Mean NN Distance:</td><td class=\"numberCell\">")
			ret.append(df.format(meanNNDist)).append("</td></tr>\n")

			/* expected mean distance; Note, this analysis is based on Davis
			 *  2002, Statistics and Data Analysis in Geology. There is an 
			 *  error in his eq 5.27 in that it should be se = 0.26136 * (A / n^2)^0.5.
			 */
			//double expectedDist = 0.5 * Math.sqrt(A / numPoints) + (0.514 + 0.412 / Math.sqrt(numPoints)) * (P / numPoints)
			double expectedDist = 0.5 * Math.sqrt(A / numPoints) // Davis' as well as Walford's equation
			//double expectedDist = 1.0d / (2.0d * Math.sqrt(numPoints / A)) // jensen and jensen formula
			ret.append("<tr><td>Expected Mean NN Distance:</td><td class=\"numberCell\">")
			ret.append(df.format(expectedDist)).append("</td></tr>").append("\n")
			
			// standard error of mean distance
			//double se = 0.070 * A / (numPoints * numPoints) + 0.035 * P * (Math.sqrt(A) / Math.pow(numPoints, (5.0 / 2.0)))
			//double se = 0.26136 * Math.sqrt(A / (numPoints * numPoints)) // Davis' equation modified
			//double se = Math.sqrt((0.06831 * A) / (numPoints * numPoints))
			double se = 0.26136 / Math.sqrt(numPoints * numPoints / A) // Jensen and Jensen as well as Walford's equation.
			ret.append("<tr><td>Standard Error:</td><td class=\"numberCell\">")
			ret.append(df.format(se)).append("</td></tr>").append("\n")

			// NN statistic
			double nnStat = meanNNDist / expectedDist
			ret.append("<tr><td>NN Ratio:</td><td class=\"numberCell\">")
			ret.append(df.format(nnStat)).append("</td></tr>")
			
			// z value
			double zVal = (meanNNDist - expectedDist) / se
			ret.append("<tr><td><i>z</i>-score:</td><td class=\"numberCell\">")
			ret.append(df.format(zVal)).append("</td></tr>").append("\n")

			NormalDistribution nd = new NormalDistribution()
			double pVal = nd.cumulativeProbability(-Math.abs(zVal)) * 2.0
			ret.append("<tr><td><i>p</i>-value:</td><td class=\"numberCell\">")
			ret.append(df.format(pVal)).append("</td></tr>")

			ret.append("</table>")
			
			// resulting statement
			if (pVal <= 0.1 && zVal < 0) {
				ret.append("<p>The distribution of points appears to be <b>significantly clustered</b>.</p>")
			} else if (pVal <= 0.1 && zVal > 0) {
				ret.append("<p>The distribution of points appears to be <b>significantly dispersed</b> in a regular grid.</p>")
			} else {
				ret.append("<p>The distribution of points appears to be <b>randomly distributed</b> and not significantly clustered nor dispersed.</p>")
			}

			ret.append("<p><i>Notes:</i><br><br>1) The nearest-neighbour ratio ranges from 0.0 for a ")
			ret.append("distribution of coincident points, to 1.0 for a random point distribution, ")
			ret.append("to a maximum of 2.15 for a hexagonal grid with maximal dispersion. Points ")
			ret.append("that are situated along a regular square grid should have a NN ratio near 2.0.<br><br>")

			if (edgeEffectStrategy.contains("nothing")) {
				ret.append("2) Nearest-neighbour analysis is strongly affected by edge effects ")
				ret.append("because points located at the edges of the extent may have their NN ")
				ret.append("situated beyond the edge. If the study ")
				ret.append("area is anomolous within the larger region, the results of this analysis ")
				ret.append("will be significantly influenced by edge effects.</p>")
			
			} else if (edgeEffectStrategy.contains("tile")) {
				ret.append("2) Nearest-neighbour analysis is strongly affected by edge effects ")
				ret.append("because points located at the edges of the extent may have their NN ")
				ret.append("situated beyond the edge. To compensate for this bias, the data were ")
				ret.append("tiled in a 3 x 3 grid and the statistic was calculated for ")
				ret.append("the central area, a procedure recommended by Davis (2002) ")
				ret.append("<i>Statistics and Data Analysis In Geology, 3rd Ed.</i> If the study ")
				ret.append("area is anomolous within the larger region, the results of this analysis ")
				ret.append("will be significantly influenced by edge effects.</p>")
			}
			
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
    def f = new NearestNeighbourAnalysis(pluginHost, args, name, descriptiveName)
}
