/*
 * Copyright (C) 2013 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
import java.util.concurrent.Future
import java.util.concurrent.*
import java.util.Date
import java.util.ArrayList
import java.util.Map
import java.util.Random
import jopensurf.*;
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import whitebox.structures.XYPoint
import whitebox.stats.PolynomialLeastSquares2DFitting
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "FindTiePoints"
def descriptiveName = "Find Tie Points"
def description = "Locate tie points between an image pair."
def toolboxes = ["Photogrammetry"]

public class FindTiePoints implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public FindTiePoints(WhiteboxPluginHost pluginHost, 
		String[] args, def descriptiveName) {
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
			def helpFile = "FindTiePoints"
			sd.setHelpFile(helpFile)
		
			// Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.
			def pathSep = File.separator
			def scriptFile = pluginHost.getResourcesDirectory() + 
			 "plugins" + pathSep + "Scripts" + pathSep + 
			   "FindTiePoints.groovy"
			sd.setSourceFile(scriptFile)
			
			// add some components to the dialog
			sd.addDialogFile("Left image", "Left Image:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Right image", "Right Image:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Output tiepoints file left", "Left Output Points File:", "save", "Shapefiles (*.shp), SHP", true, false)
			sd.addDialogFile("Output tiepoints file right", "Right Output Points File:", "save", "Shapefiles (*.shp), SHP", true, false)
			
			sd.addDialogLabel(" ")
			sd.addDialogLabel("<html><b><i>SURF Feature Detection Parameters:</i></b></html>")
			sd.addDialogDataInput("Threshold value", "Threshold Value:", "4.0", true, false)
            sd.addDialogDataInput("Number of octaves", "Number of Octaves:", "4", true, false)
            //sd.addDialogDataInput("Balance value", "Balance Value:", "0.81", true, false)
			sd.addDialogDataInput("SURF point matching threshold value", "Matching Threshold Value:", "0.6", true, false)

            sd.addDialogLabel(" ")
			sd.addDialogLabel("<html><b><i>Model-Fitting Parameters:</i></b></html>")
			sd.addDialogDataInput("Threshold value for determining when a datum fits a model", "Threshold Value For Removal:", "2.0", true, false)
            sd.addDialogDataInput("Polynomial order", "Polynomial Order:", "1", true, false)
            
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
	  	
			if (args.length != 9) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}

			whitebox.geospatialfiles.shapefile.Point outputPoint;
			Object[] rowData
			int progress, oldProgress
			int iteration, a, b, i, k
			int n = 10
			double x, y, x2, y2
			
			// read the input parameters
			String leftFile = args[0]
			String rightFile = args[1]
			String leftOutputFile = args[2]
			String rightOutputFile = args[3]
			float surfThreshold = (float)(Double.parseDouble(args[4]) / 1000.0)
			int surfNumOctaves = Integer.parseInt(args[5])
			float surfBalanceValue = 0.81; //Float.parseFloat(args[6])
			double matchingThreshold = Double.parseDouble(args[6])
			double ransacThreshold = Double.parseDouble(args[7])
			int polyOrder = Integer.parseInt(args[8])
			if (polyOrder > 5) { polyOrder = 5; }

			// read in the left and right images
			WhiteboxRaster leftImage = new WhiteboxRaster(leftFile, "r");
			double leftNodata = leftImage.getNoDataValue();
			WhiteboxRaster rightImage = new WhiteboxRaster(rightFile, "r");
			double rightNodata = rightImage.getNoDataValue();

			// perform the SURF analysis
			pluginHost.updateProgress("Performing SURF analysis on left image...", 0);
            Surf leftSurf = new Surf(leftImage, surfBalanceValue, surfThreshold, surfNumOctaves);
            if (leftSurf.getNumberOfPoints() == 0) {
                pluginHost.showFeedback("Number of points equals zero, reset threshold: " + leftSurf.getNumberOfPoints());
                return;
            }

            pluginHost.updateProgress("Performing SURF analysis on right image...", 0);
            Surf rightSurf = new Surf(rightImage, surfBalanceValue, surfThreshold, surfNumOctaves);
            if (rightSurf.getNumberOfPoints() == 0) {
                pluginHost.showFeedback("Number of points equals zero, reset threshold: " + leftSurf.getNumberOfPoints());
                return;
            }

			// find the matching points between each SURF point set
            pluginHost.updateProgress("Matching points of interest...", 0);
            Map<SURFInterestPoint, SURFInterestPoint> matchingPoints =
                    leftSurf.getMatchingPoints(rightSurf, matchingThreshold, false);

            int numTiePoints = matchingPoints.size();
            if (numTiePoints < 3) {
                pluginHost.showFeedback("The number of potential tie points is less than 3. Adjust your threshold parameters and retry.");
                return;
            }

            // use RANSAC to find cooresponding points and eliminate outliers
			pluginHost.updateProgress("Trimming outlier tie points...", 0);

			XYPoint[][] tiePointsList = new XYPoint[numTiePoints][2]
			ArrayList<XYPoint> leftTiePointsList = new ArrayList<>();
            ArrayList<XYPoint> rightTiePointsList = new ArrayList<>();
            ArrayList<XYPoint> leftConsensusSet = new ArrayList<>();
            ArrayList<XYPoint> rightConsensusSet = new ArrayList<>();
            boolean[] includedList;

			// place the tie points into an array 
            i = 0;
            for (SURFInterestPoint point : matchingPoints.keySet()) {
                x = point.getX();
                y = point.getY();

                SURFInterestPoint target = matchingPoints.get(point);
                x2 = target.getX();
                y2 = target.getY();

                tiePointsList[i][0] = new XYPoint(x, y); // left point
                tiePointsList[i][1] = new XYPoint(x2, y2); // right point
                i++;
            }

			PolynomialLeastSquares2DFitting bestModel
			double overallError = Double.POSITIVE_INFINITY
			// add the points to the consensus set
			for (a = 0; a < numTiePoints; a++) {
				leftConsensusSet.add(tiePointsList[a][0]);
				rightConsensusSet.add(tiePointsList[a][1]);
			}

			double maxError = Double.POSITIVE_INFINITY
			double[] residuals
			while (maxError > ransacThreshold) {
				if (leftConsensusSet.size() <= n) {
					pluginHost.showFeedback("The error threshold appears to be too high. Please reset it and try again.")
					break
				}
				// create the model
				bestModel = new PolynomialLeastSquares2DFitting(leftConsensusSet, rightConsensusSet, polyOrder);
				residuals = bestModel.getResidualsXY()

				// find the maximum error
				int maxErrorIndex = -1
				maxError = 0
				for (a = 0; a < residuals.length; a++) {
					if (residuals[a] > maxError) {
						maxError = residuals[a]
						maxErrorIndex = a
					}
				}

				if (maxError > ransacThreshold) {
					// remove the point with the greatest error
					leftConsensusSet.remove(maxErrorIndex)
					rightConsensusSet.remove(maxErrorIndex)
				}
			}

			overallError = bestModel.getOverallRMSE()
			
//			PolynomialLeastSquares2DFitting bestModel
//			int largestConsensusSet = 0;
//			double overallError = Double.POSITIVE_INFINITY;
//            
//			oldProgress = -1
//            for (iteration = 0; iteration < numIterations; iteration++) {
//            	leftTiePointsList.clear();
//            	rightTiePointsList.clear();
//				includedList = new boolean[numTiePoints]
//				
//            	// choose some points at random to create a model from
//            	a = 0
//            	k = 0
//            	while (a < n) {
//            		b = (int)(Math.random() * numTiePoints)
//            		
//					if (!includedList[b]) {
//						leftTiePointsList.add(tiePointsList[b][0]);
//						rightTiePointsList.add(tiePointsList[b][1]);
//						includedList[b] = true;
//						a++
//					}
//					k++
//					if (k == numTiePoints) {
//						break
//					}
//            	}
//
//            	PolynomialLeastSquares2DFitting maybeModel = new PolynomialLeastSquares2DFitting(leftTiePointsList, rightTiePointsList, polyOrder);
//
//				for (a = 0; a < numTiePoints; a++) {
//					if (!includedList[a]) {
//						includedList[a] = true;
//						x = tiePointsList[a][0].x;
//	                    y = tiePointsList[a][0].y;
//	                    XYPoint pt1 = maybeModel.getForwardCoordinates(x, y);
//	                    XYPoint pt2 = new XYPoint(tiePointsList[a][1].x, tiePointsList[a][1].y)
//	                    double error = pt1.getSquareDistance(pt2);
//		                if (error < ransacThreshold) {
//		                    // add the point to the list
//		                    leftTiePointsList.add(tiePointsList[a][0]);
//							rightTiePointsList.add(tiePointsList[a][1]);
//		                }
//					}
//				}
//			 	maybeModel = new PolynomialLeastSquares2DFitting(leftTiePointsList, rightTiePointsList, polyOrder);
//            	double modelError = maybeModel.getOverallRMSE();
//            	
//				if (leftTiePointsList.size() > largestConsensusSet) { //modelError < overallError) {
//					largestConsensusSet = leftTiePointsList.size()
//					overallError = modelError
//					bestModel = new PolynomialLeastSquares2DFitting(leftTiePointsList, rightTiePointsList, polyOrder);
//        	  		leftConsensusSet = new ArrayList<XYPoint>(leftTiePointsList.size());
//					for (XYPoint item: leftTiePointsList) leftConsensusSet.add(new XYPoint(item.x, item.y));
//					rightConsensusSet = new ArrayList<XYPoint>(rightTiePointsList.size());
//					for (XYPoint item: rightTiePointsList) rightConsensusSet.add(new XYPoint(item.x, item.y));
//				}
//				
//				progress = (int)(100f * iteration / (numIterations - 1))
//				if (progress > oldProgress) {
//					oldProgress = progress
//					pluginHost.updateProgress("Trimming outlier tie points...", progress)
//				}
//            }
            
            // output the final points to the output shapefiles
			pluginHost.updateProgress("Outputing tie point files...", 0);

			DBFField[] fields = new DBFField[2];
            fields[0] = new DBFField();
            fields[0].setName("FID");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setDecimalCount(4);
            fields[0].setFieldLength(10);

            fields[1] = new DBFField();
            fields[1].setName("RESID");
            fields[1].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[1].setDecimalCount(4);
            fields[1].setFieldLength(10);


            ShapeFile leftOutput = new ShapeFile(leftOutputFile, ShapeType.POINT, fields);
            ShapeFile rightOutput = new ShapeFile(rightOutputFile, ShapeType.POINT, fields);

			for (a = 0; a < leftConsensusSet.size(); a++) {
				XYPoint p1 = leftConsensusSet.get(a)
                x = leftImage.getXCoordinateFromColumn((int) p1.x);
                y = leftImage.getYCoordinateFromRow((int) p1.y);

                XYPoint p2 = rightConsensusSet.get(a)
                x2 = rightImage.getXCoordinateFromColumn((int) p2.x);
                y2 = rightImage.getYCoordinateFromRow((int) p2.y);

                outputPoint = new whitebox.geospatialfiles.shapefile.Point(x, y);
                rowData = new Object[2];
                rowData[0] = new Double(a);
                rowData[1] = new Double(residuals[a]);
                leftOutput.addRecord(outputPoint, rowData);

                outputPoint = new whitebox.geospatialfiles.shapefile.Point(x2, y2);
                rowData = new Object[2];
                rowData[0] = new Double(a);
                rowData[1] = new Double(residuals[a]);
                rightOutput.addRecord(outputPoint, rowData);
            }

            leftOutput.write();
            rightOutput.write();

			n = rightConsensusSet.size()
			
            pluginHost.showFeedback(n + " points located, RMSE = " + overallError + ", Max Error = " + maxError);
		
	  	} catch (OutOfMemoryError oe) {
            pluginHost.showFeedback("An out-of-memory error has occurred during operation.")
	    } catch (Exception e) {
	        pluginHost.showFeedback("An error has occurred during operation. See log file for details.")
	        pluginHost.logException("Error in " + descriptiveName, e)
        } finally {
        	// reset the progress bar
        	pluginHost.updateProgress("Progress:", 0)
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
	def f = new FindTiePoints(pluginHost, args, descriptiveName)
}
