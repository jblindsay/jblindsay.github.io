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
import java.io.File
import java.util.Date
import java.util.ArrayList
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.FileUtilities
import whitebox.structures.BoundingBox
import whitebox.geospatialfiles.LASReader
import whitebox.geospatialfiles.LASReader.PointRecord
import whitebox.geospatialfiles.VectorLayerInfo
import whitebox.structures.KdTree
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "Hexbinning"
def descriptiveName = "Hex-binning"
def description = "Creates a hex-grid density plot (heat map)"
def toolboxes = ["VectorTools", "LidarTools"]

public class Hexbinning implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
	
    public Hexbinning(WhiteboxPluginHost pluginHost, 
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
            sd.addDialogFile("Input points file", "Input Points File:", "open", "Whitebox Files (*.las; *.shp), LAS, SHP", true, false)
            sd.addDialogFile("Output file", "Output Vector File:", "close", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogDataInput("Enter a grid resolution", "Hexagon Width:", "", true, false)
            sd.addDialogComboBox("Grid direction", "Grid Direction:", ["horizontal", "vertical"], 0)
            
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

            if (args.length != 4) {
                pluginHost.showFeedback("Incorrect number of arguments given to tool.")
                return
            }
            // read the input parameters
            String baseFile = args[0]
            String outputFile = args[1]
            double width = Double.parseDouble(args[2])
            String orientation = "h"
            if (args[3].toLowerCase().contains("v")) {
            	orientation = "v"
            }

            final double sixtyDegrees = Math.PI / 6

            double halfWidth = 0.5 * width
			double size = halfWidth / Math.cos(sixtyDegrees)
			double height = size * 2
			double threeQuarterHeight = 0.75 * height
			int row, col
			int FID = 0
			int progress
			int oldProgress = -1
				

            int[] parts = [0];

			// figure out the extent of the base data set
			BoundingBox extent
			if (baseFile.toLowerCase().contains(".las")) {

				// set up the output files of the shapefile and the dbf
	            DBFField[] fields = new DBFField[6];
	
	            fields[0] = new DBFField();
	            fields[0].setName("FID");
	            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
	            fields[0].setFieldLength(10);
	            fields[0].setDecimalCount(0);
	
	            fields[1] = new DBFField();
	            fields[1].setName("COUNT");
	            fields[1].setDataType(DBFField.DBFDataType.NUMERIC);
	            fields[1].setFieldLength(10);
	            fields[1].setDecimalCount(0);
	
	            fields[2] = new DBFField();
	            fields[2].setName("ROW");
	            fields[2].setDataType(DBFField.DBFDataType.NUMERIC);
	            fields[2].setFieldLength(10);
	            fields[2].setDecimalCount(0);
	
	            fields[3] = new DBFField();
	            fields[3].setName("COLUMN");
	            fields[3].setDataType(DBFField.DBFDataType.NUMERIC);
	            fields[3].setFieldLength(10);
	            fields[3].setDecimalCount(0);

	            fields[4] = new DBFField();
	            fields[4].setName("MIN_Z");
	            fields[4].setDataType(DBFField.DBFDataType.NUMERIC);
	            fields[4].setFieldLength(12);
	            fields[4].setDecimalCount(5);

	            fields[5] = new DBFField();
	            fields[5].setName("MAX_Z");
	            fields[5].setDataType(DBFField.DBFDataType.NUMERIC);
	            fields[5].setFieldLength(12);
	            fields[5].setDecimalCount(5);
	            
	            ShapeFile output = new ShapeFile(outputFile, ShapeType.POLYGON, fields);


				LASReader las = new LASReader(baseFile)
				extent = new BoundingBox(las.getMinX(), las.getMinY(), las.getMaxX(), las.getMaxY());
			
	            if (orientation.equals("h")) { // horizontal orientation
		            
		            double center_x_0 = extent.getMinX() + halfWidth
					double center_y_0 = extent.getMaxY() - 0.25 * height
					double center_x, center_y
					int rows = (int)Math.ceil(extent.getHeight() / threeQuarterHeight)
					int cols = (int)Math.ceil(extent.getWidth() / width)
					if (rows * cols > 100000) {
						int ret = pluginHost.showFeedback("This operation will produce a vector file with a very large number of \npolygons which may be difficult to handle. Do you want to proceed?", 1, 1)
						if (ret == 1) {
							return
						}
					}
					
					int numHexagons = 0;
					for (row = 0; row < rows; row++) {
						cols = (int)Math.ceil((extent.getWidth() + halfWidth * (row % 2)) / width)
						for (col = 0; col < cols; col++) {
							numHexagons++
						}
					}
					
					KdTree<Integer> pointsTree = new KdTree.SqrEuclid<Integer>(2, new Integer(numHexagons))
					int index = 0
					for (row = 0; row < rows; row++) {
						center_y = center_y_0 - row * threeQuarterHeight
						cols = (int)Math.ceil((extent.getWidth() + halfWidth * (row % 2)) / width)
						for (col = 0; col < cols; col++) {
							center_x = (center_x_0 - halfWidth * (row % 2)) + col * width
							double[] entry = [center_y, center_x]
							pointsTree.addPoint(entry, new Integer(index));
							index++;
						}
					}
					int[] counts = new int[numHexagons]
					double[] minZ = new double[numHexagons]
					double[] maxZ = new double[numHexagons]
					for (int a = 0; a < numHexagons; a++) { 
						minZ[a] = Double.POSITIVE_INFINITY;
						maxZ[a] = Double.NEGATIVE_INFINITY; 
					}
					
					int numPoints = (int) las.getNumPointRecords();
					PointRecord point;
	    			//List<KdTree.Entry<Integer>> result
	    			KdTree.Entry<Integer> result
	    			double x, y, z
					for (int a = 0; a < numPoints; a++) {
		                point = las.getPointRecord(a);
		                if (!point.isPointWithheld()) {
		                    x = point.getX();
		                    y = point.getY();
		                    z = point.getZ();
		                    double[] entry = [y, x];
		                    //result = pointsTree.nearestNeighbor(entry), 1, true, false);
		                    result = pointsTree.nearestNeighbor(entry);
		                    //index = (int)result.get(0).value;
		                    index = (int)result.value;
		                    counts[index]++
		                    if (z < minZ[index]) { minZ[index] = z; }
		                    if (z > maxZ[index]) { maxZ[index] = z; }
		                }
		                progress = (int) (100f * (a + 1) / numPoints);
		                if (progress != oldProgress) {
		                    oldProgress = progress;
		                    pluginHost.updateProgress("Reading point data:", progress);
		                    if (pluginHost.isRequestForOperationCancelSet()) {
		                        pluginHost.showFeedback("Operation cancelled")
								return
		                    }
		                }
		            }
		            
					for (row = 0; row < rows; row++) {
						center_y = center_y_0 - row * threeQuarterHeight
						cols = (int)Math.ceil((extent.getWidth() + halfWidth * (row % 2)) / width)
						for (col = 0; col < cols; col++) {
							if (counts[FID] > 0) {
								PointsList points = new PointsList();
								center_x = (center_x_0 - halfWidth * (row % 2)) + col * width
								for (int i = 6; i >= 0; i--) {
					    			double angle = 2 * sixtyDegrees * (i + 0.5)
					    			double x_i = center_x + size * Math.cos(angle)
					    			double y_i = center_y + size * Math.sin(angle)
					
					    			points.addPoint(x_i, y_i);
					            }
					            
					            Polygon poly = new Polygon(parts, points.getPointsArray());
					            Object[] rowData = new Object[6];
					            rowData[0] = new Double(FID);
					            rowData[1] = new Double(counts[FID]);
					            rowData[2] = new Double(row);
					            rowData[3] = new Double(col);
					            rowData[4] = new Double(minZ[FID]);
					            rowData[5] = new Double(maxZ[FID]);
					            output.addRecord(poly, rowData);
							}
				            FID++
						}
						progress = (int)(100f * row / (rows - 1))
	        			if (progress != oldProgress) {
							pluginHost.updateProgress(progress)
	            			oldProgress = progress
	            		}
	            		// check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
					}
	            } else { // vertical orientation
	            	double center_x_0 = extent.getMinX() + 0.5 * size
					double center_y_0 = extent.getMaxY() - halfWidth
					double center_x, center_y
					int rows = (int)Math.ceil(extent.getHeight() / width)
					int cols = (int)Math.ceil(extent.getWidth() / threeQuarterHeight)
					if (rows * cols > 100000) {
						int ret = pluginHost.showFeedback("This operation will produce a vector file with a very large number of \npolygons which may be difficult to handle. Do you want to proceed?", 1, 1)
						if (ret == 1) {
							return
						}
					}

					int numHexagons = 0;
					for (col = 0; col < cols; col++) {
						rows = (int)Math.ceil((extent.getHeight() + ((col % 2) * halfWidth)) / width)
						for (row = 0; row < rows; row++) {					
							numHexagons++
						}
					}
					
					KdTree<Integer> pointsTree = new KdTree.SqrEuclid<Integer>(2, new Integer(numHexagons))
					int index = 0
					for (col = 0; col < cols; col++) {
						rows = (int)Math.ceil((extent.getHeight() + ((col % 2) * halfWidth)) / width)
						for (row = 0; row < rows; row++) {					
							center_y = center_y_0 - row * width + ((col % 2) * halfWidth)
							center_x = center_x_0 + col * threeQuarterHeight
							double[] entry = [center_y, center_x]
							pointsTree.addPoint(entry, new Integer(index));
							index++;
						}
					}
					int[] counts = new int[numHexagons]
					double[] minZ = new double[numHexagons]
					double[] maxZ = new double[numHexagons]
					for (int a = 0; a < numHexagons; a++) { 
						minZ[a] = Double.POSITIVE_INFINITY;
						maxZ[a] = Double.NEGATIVE_INFINITY; 
					}

					int numPoints = (int) las.getNumPointRecords();
					PointRecord point;
	    			//List<KdTree.Entry<Integer>> result
	    			KdTree.Entry<Integer> result
	    			double x, y, z
					for (int a = 0; a < numPoints; a++) {
		                point = las.getPointRecord(a);
		                if (!point.isPointWithheld()) {
		                    x = point.getX();
		                    y = point.getY();
		                    
		                    double[] entry = [y, x];
		                    //result = pointsTree.nearestNeighbor(entry, 1, true, false);
		                    //index = (int)result.get(0).value;
		                    result = pointsTree.nearestNeighbor(entry);
		                    index = (int)result.value;
		                    counts[index]++
		                    if (z < minZ[index]) { minZ[index] = z; }
		                    if (z > maxZ[index]) { maxZ[index] = z; }
		                }
		                progress = (int) (100f * (a + 1) / numPoints);
		                if (progress != oldProgress) {
		                    oldProgress = progress;
		                    pluginHost.updateProgress("Reading point data:", progress);
		                    if (pluginHost.isRequestForOperationCancelSet()) {
		                        pluginHost.showFeedback("Operation cancelled")
								return
		                    }
		                }
		            }
					
					for (col = 0; col < cols; col++) {
						rows = (int)Math.ceil((extent.getHeight() + ((col % 2) * halfWidth)) / width)
						for (row = 0; row < rows; row++) {	
							if (counts[FID] > 0) {				
								center_y = center_y_0 - row * width + ((col % 2) * halfWidth)
								
								PointsList points = new PointsList();
								center_x = center_x_0 + col * threeQuarterHeight
								for (int i = 6; i >= 0; i--) {
					    			double angle = 2 * sixtyDegrees * (i + 0.5) - sixtyDegrees
					    			double x_i = center_x + size * Math.cos(angle)
					    			double y_i = center_y + size * Math.sin(angle)
					
					    			points.addPoint(x_i, y_i);
					            }
					            
					            Polygon poly = new Polygon(parts, points.getPointsArray());
					            Object[] rowData = new Object[6];
					            rowData[0] = new Double(FID);
					            rowData[1] = new Double(counts[FID]);
					            rowData[2] = new Double(row);
					            rowData[3] = new Double(col);
					            rowData[4] = new Double(minZ[FID]);
					            rowData[5] = new Double(maxZ[FID]);
					            output.addRecord(poly, rowData);
							}
				            FID++
						}
						progress = (int)(100f * col / (cols - 1))
	        			if (progress != oldProgress) {
							pluginHost.updateProgress(progress)
	            			oldProgress = progress
	            		}
	            		// check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
					}
	            }
				
	            output.write();
		           
//	            // display the output image
//	            pluginHost.returnData(outputFile)			
			
			} else if (baseFile.toLowerCase().contains(".shp")) {
				// set up the output files of the shapefile and the dbf
	            DBFField[] fields = new DBFField[4];
	
	            fields[0] = new DBFField();
	            fields[0].setName("FID");
	            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
	            fields[0].setFieldLength(10);
	            fields[0].setDecimalCount(0);
	
	            fields[1] = new DBFField();
	            fields[1].setName("COUNT");
	            fields[1].setDataType(DBFField.DBFDataType.NUMERIC);
	            fields[1].setFieldLength(10);
	            fields[1].setDecimalCount(0);
	
	            fields[2] = new DBFField();
	            fields[2].setName("ROW");
	            fields[2].setDataType(DBFField.DBFDataType.NUMERIC);
	            fields[2].setFieldLength(10);
	            fields[2].setDecimalCount(0);
	
	            fields[3] = new DBFField();
	            fields[3].setName("COLUMN");
	            fields[3].setDataType(DBFField.DBFDataType.NUMERIC);
	            fields[3].setFieldLength(10);
	            fields[3].setDecimalCount(0);
	            
	            ShapeFile output = new ShapeFile(outputFile, ShapeType.POLYGON, fields);

				ShapeFile sf = new ShapeFile(baseFile)
				ShapeType shapeType = sf.getShapeType()
				if (shapeType.getBaseType() != ShapeType.POINT &&
				 shapeType.getBaseType() != ShapeType.MULTIPOINT) {
				 	pluginHost.showFeedback("The input shapefile must be of a point or multipoint shapetype.")
				 	return
				}
				
				extent = new BoundingBox(sf.getxMin(), sf.getyMin(), sf.getxMax(), sf.getyMax())

				if (orientation.equals("h")) { // horizontal orientation
		            
		            double center_x_0 = extent.getMinX() + halfWidth
					double center_y_0 = extent.getMaxY() - 0.25 * height
					double center_x, center_y
					int rows = (int)Math.ceil(extent.getHeight() / threeQuarterHeight)
					int cols = (int)Math.ceil(extent.getWidth() / width)
					if (rows * cols > 100000) {
						int ret = pluginHost.showFeedback("This operation will produce a vector file with a very large number of \npolygons which may be difficult to handle. Do you want to proceed?", 1, 1)
						if (ret == 1) {
							return
						}
					}
					
					int numHexagons = 0;
					for (row = 0; row < rows; row++) {
						cols = (int)Math.ceil((extent.getWidth() + halfWidth * (row % 2)) / width)
						for (col = 0; col < cols; col++) {
							numHexagons++
						}
					}
					
					KdTree<Integer> pointsTree = new KdTree.SqrEuclid<Integer>(2, new Integer(numHexagons))
					int index = 0
					for (row = 0; row < rows; row++) {
						center_y = center_y_0 - row * threeQuarterHeight
						cols = (int)Math.ceil((extent.getWidth() + halfWidth * (row % 2)) / width)
						for (col = 0; col < cols; col++) {
							center_x = (center_x_0 - halfWidth * (row % 2)) + col * width
							double[] entry = [center_y, center_x]
							pointsTree.addPoint(entry, new Integer(index));
							index++;
						}
					}
					int[] counts = new int[numHexagons]
					
					int numPoints;
					PointRecord point;
	    			List<KdTree.Entry<Integer>> result
	    			double x, y
	    			if (shapeType.getBaseType() == ShapeType.POINT) {
	    				numPoints = sf.getNumberOfRecords();
		    			for (ShapeFileRecord record : sf.records) {
		    				int a = record.getRecordNumber()
			                double[][] points = record.getGeometry().getPoints()
			                x = points[0][0]
			                y = points[0][1]
			                double[] entry = [y, x];
		                    result = pointsTree.nearestNeighbor(entry, 1, true, false);
		                    index = (int)result.get(0).value;
		                    counts[index]++
			                progress = (int) (100f * (a + 1) / numPoints);
			                if (progress != oldProgress) {
			                    oldProgress = progress;
			                    pluginHost.updateProgress("Reading point data:", progress);
			                    if (pluginHost.isRequestForOperationCancelSet()) {
			                        pluginHost.showFeedback("Operation cancelled")
									return
			                    }
			                }
			            }
	    			} else if (shapeType.getBaseType() == ShapeType.MULTIPOINT) {
		    			for (ShapeFileRecord record : sf.records) {
			                double[][] points = record.getGeometry().getPoints()
			                numPoints = points.length
			                for (int i = 0; i < numPoints; i++) {
			                	x = points[i][0]
			                	y = points[i][1]
			               		double[] entry = [y, x];
		                        result = pointsTree.nearestNeighbor(entry, 1, true, false);
			                    index = (int)result.get(0).value;
			                    counts[index]++

			                    progress = (int) (100f * (i + 1) / numPoints);
				                if (progress != oldProgress) {
				                    oldProgress = progress;
				                    pluginHost.updateProgress("Reading point data:", progress);
				                    if (pluginHost.isRequestForOperationCancelSet()) {
				                        pluginHost.showFeedback("Operation cancelled")
										return
				                    }
				                }
			                }
			            }
	    			}
		            
					for (row = 0; row < rows; row++) {
						center_y = center_y_0 - row * threeQuarterHeight
						cols = (int)Math.ceil((extent.getWidth() + halfWidth * (row % 2)) / width)
						for (col = 0; col < cols; col++) {
							if (counts[FID] > 0) {
								PointsList points = new PointsList();
								center_x = (center_x_0 - halfWidth * (row % 2)) + col * width
								for (int i = 6; i >= 0; i--) {
					    			double angle = 2 * sixtyDegrees * (i + 0.5)
					    			double x_i = center_x + size * Math.cos(angle)
					    			double y_i = center_y + size * Math.sin(angle)
					
					    			points.addPoint(x_i, y_i);
					            }
					            
					            Polygon poly = new Polygon(parts, points.getPointsArray());
					            Object[] rowData = new Object[4];
					            rowData[0] = new Double(FID);
					            rowData[1] = new Double(counts[FID]);
					            rowData[2] = new Double(row);
					            rowData[3] = new Double(col);
					            output.addRecord(poly, rowData);
							}
				            FID++
						}
						progress = (int)(100f * row / (rows - 1))
	        			if (progress != oldProgress) {
							pluginHost.updateProgress(progress)
	            			oldProgress = progress
	            		}
	            		// check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
					}
	            } else { // vertical orientation
	            	double center_x_0 = extent.getMinX() + 0.5 * size
					double center_y_0 = extent.getMaxY() - halfWidth
					double center_x, center_y
					int rows = (int)Math.ceil(extent.getHeight() / width)
					int cols = (int)Math.ceil(extent.getWidth() / threeQuarterHeight)
					if (rows * cols > 100000) {
						int ret = pluginHost.showFeedback("This operation will produce a vector file with a very large number of \npolygons which may be difficult to handle. Do you want to proceed?", 1, 1)
						if (ret == 1) {
							return
						}
					}

					int numHexagons = 0;
					for (col = 0; col < cols; col++) {
						rows = (int)Math.ceil((extent.getHeight() + ((col % 2) * halfWidth)) / width)
						for (row = 0; row < rows; row++) {					
							numHexagons++
						}
					}
					
					KdTree<Integer> pointsTree = new KdTree.SqrEuclid<Integer>(2, new Integer(numHexagons))
					int index = 0
					for (col = 0; col < cols; col++) {
						rows = (int)Math.ceil((extent.getHeight() + ((col % 2) * halfWidth)) / width)
						for (row = 0; row < rows; row++) {					
							center_y = center_y_0 - row * width + ((col % 2) * halfWidth)
							center_x = center_x_0 + col * threeQuarterHeight
							double[] entry = [center_y, center_x]
							pointsTree.addPoint(entry, new Integer(index));
							index++;
						}
					}
					int[] counts = new int[numHexagons]

					int numPoints;
					PointRecord point;
	    			List<KdTree.Entry<Integer>> result
	    			double x, y
	    			if (shapeType.getBaseType() == ShapeType.POINT) {
	    				numPoints = sf.getNumberOfRecords();
		    			for (ShapeFileRecord record : sf.records) {
		    				int a = record.getRecordNumber()
			                double[][] points = record.getGeometry().getPoints()
			                x = points[0][0]
			                y = points[0][1]
			                double[] entry = [y, x];
		                    result = pointsTree.nearestNeighbor(entry, 1, true, false);
		                    index = (int)result.get(0).value;
		                    counts[index]++
			                progress = (int) (100f * (a + 1) / numPoints);
			                if (progress != oldProgress) {
			                    oldProgress = progress;
			                    pluginHost.updateProgress("Reading point data:", progress);
			                    if (pluginHost.isRequestForOperationCancelSet()) {
			                        pluginHost.showFeedback("Operation cancelled")
									return
			                    }
			                }
			            }
	    			} else if (shapeType.getBaseType() == ShapeType.MULTIPOINT) {
		    			for (ShapeFileRecord record : sf.records) {
			                double[][] points = record.getGeometry().getPoints()
			                numPoints = points.length
			                for (int i = 0; i < numPoints; i++) {
			                	x = points[i][0]
			                	y = points[i][1]
			               		double[] entry = [y, x];
		                        result = pointsTree.nearestNeighbor(entry, 1, true, false);
			                    index = (int)result.get(0).value;
			                    counts[index]++

			                    progress = (int) (100f * (i + 1) / numPoints);
				                if (progress != oldProgress) {
				                    oldProgress = progress;
				                    pluginHost.updateProgress("Reading point data:", progress);
				                    if (pluginHost.isRequestForOperationCancelSet()) {
				                        pluginHost.showFeedback("Operation cancelled")
										return
				                    }
				                }
			                }
			            }
	    			}
		            
					for (col = 0; col < cols; col++) {
						rows = (int)Math.ceil((extent.getHeight() + ((col % 2) * halfWidth)) / width)
						for (row = 0; row < rows; row++) {	
							if (counts[FID] > 0) {				
								center_y = center_y_0 - row * width + ((col % 2) * halfWidth)
								
								PointsList points = new PointsList();
								center_x = center_x_0 + col * threeQuarterHeight
								for (int i = 6; i >= 0; i--) {
					    			double angle = 2 * sixtyDegrees * (i + 0.5) - sixtyDegrees
					    			double x_i = center_x + size * Math.cos(angle)
					    			double y_i = center_y + size * Math.sin(angle)
					
					    			points.addPoint(x_i, y_i);
					            }
					            
					            Polygon poly = new Polygon(parts, points.getPointsArray());
					            Object[] rowData = new Object[4];
					            rowData[0] = new Double(FID);
					            rowData[1] = new Double(counts[FID]);
					            rowData[2] = new Double(row);
					            rowData[3] = new Double(col);
					            output.addRecord(poly, rowData);
							}
				            FID++
						}
						progress = (int)(100f * col / (cols - 1))
	        			if (progress != oldProgress) {
							pluginHost.updateProgress(progress)
	            			oldProgress = progress
	            		}
	            		// check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
					}
	            }
				
	            output.write();
		           
//	            
//	            pluginHost.returnData(outputFile)

			} else {
				pluginHost.showFeedback("The input base file is of an unrecognized type.")
				return
			}

			// display the output image
			String paletteDirectory = pluginHost.getResourcesDirectory() + "palettes" + File.separator;
			VectorLayerInfo vli = new VectorLayerInfo(outputFile, paletteDirectory, 255i, -1)
			vli.setPaletteFile(paletteDirectory + "light_quant.pal") // "imhof2.pal")
			vli.setFilledWithOneColour(false)
			vli.setFillAttribute("COUNT")
			vli.setPaletteScaled(true)
			vli.setNonlinearity(0.5)
			vli.setRecordsColourData()
			pluginHost.returnData(vli)
            
        } catch (Exception e) {
            pluginHost.showFeedback(e.getMessage())
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
    def f = new Hexbinning(pluginHost, args, name, descriptiveName)
}
