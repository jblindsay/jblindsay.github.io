/*
 * Copyright (C) 2015 Dr. John Lindsay <jlindsay@uoguelph.ca>
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

 /* The following is an experiment for a paper that I am drafting. It
  *  is not yet ready for general use.
  */
import java.awt.event.ActionListener
import java.awt.event.ActionEvent
import java.beans.PropertyChangeEvent
import java.beans.PropertyChangeListener
import java.util.Date
import java.util.ArrayList
import java.util.PriorityQueue
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.WhiteboxRasterBase.DataScale
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.ShapeFileRecord
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.ui.plugin_dialog.*
import whitebox.utilities.StringUtilities;
import whitebox.structures.KdTree;
import whitebox.structures.XYPoint;
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "VectorStreamNetworkAnalysis"
def descriptiveName = "Vector Stream Network Analysis"
//def description = "Vector stream network analysis"
//def toolboxes = ["StreamAnalysis"]

public class VectorStreamNetworkAnalysis implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public VectorStreamNetworkAnalysis(WhiteboxPluginHost pluginHost, 
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
			sd.addDialogFile("Input streams file", "Input Streams File:", "open", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogFile("Input DEM raster", "Input DEM Raster:", "open", "Raster Files (*.dep), DEP", true, false)
            sd.addDialogFile("Output file", "Output File:", "save", "Vector Files (*.shp), SHP", true, false)
			sd.addDialogDataInput("Maximum ridge-cutting height (z units)", "Maximum ridge-cutting height (z-units):", "10.0", true, false)
            sd.addDialogDataInput("Snap distance for dangling streams (xy units)", "Snap distance for dangling streams (xy units):", "0.001", true, false)
            
			// resize the dialog to the standard size and display it
			sd.setSize(800, 400)
			sd.visible = true
		}
	}

	@CompileStatic
	private void execute(String[] args) {
		try {
			
	  		int progress, oldProgress, col, row;
	  		int n, j;
	  		double x, y, z, z1, z2, x1, x2, y1, y2;
	  		double searchDist = 0.000000 // This has to be a small non-zero value and is used in the nearest-neighbour search.
	  		double length;
  		  	double distMultiplier = 1.0;
  		  	Object[] rowData;
			double[] linkMag;
			boolean[] isBeyondEdge;
			double[] linkLengths;
			int[] outletNums;
			int[] numDownstreamNodes;
			KdTree<Integer> pointsTree;

  		  	
	  		if (args.length < 5) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String streamsFile = args[0];
			String demFile = args[1];
			String outputFile = args[2];
			double maxRidgeCuttingHeight = Double.parseDouble(args[3]);
			double snapDistance = Double.parseDouble(args[4]);
			snapDistance = snapDistance * snapDistance;
			
			
			// read the input image
			WhiteboxRaster dem = new WhiteboxRaster(demFile, "r");
			dem.setForceAllDataInMemory(true)
			double nodata = dem.getNoDataValue();
			int rows = dem.getNumberRows()
			int cols = dem.getNumberColumns()
			int rowsLessOne = rows - 1
			int colsLessOne = cols - 1

			if (dem.getXYUnits().toLowerCase().contains("deg")) {
				double midLat = (dem.getNorth() - dem.getSouth()) / 2.0;
                if (midLat <= 90 && midLat >= -90) {
                    midLat = Math.toRadians(midLat);
                	double a = 6378137.0; 
					double b = 6356752.314;
					double e2 = (a*a - b*b) / (a*a)
					double num = (Math.PI * a * Math.cos(midLat));
					double denum = (180 * Math.sqrt((1 - e2 * Math.sin(midLat) * Math.sin(midLat))))
					double longDegDist = (num / denum);
					double latDegDist = 111132.954 - 559.822 * Math.cos(2.0 * midLat) + 1.175 * Math.cos(4.0 * midLat)
					distMultiplier = (longDegDist + latDegDist) / 2.0;
                }
			}

			ShapeFile input = new ShapeFile(streamsFile)
			ShapeType shapeType = input.getShapeType()
            if (shapeType.getBaseType() != ShapeType.POLYLINE) {
            	pluginHost.showFeedback("The input shapefile should be of a POLYLINE ShapeType.")
            	return
            }

            // create the output file
            DBFField[] fields = new DBFField[12];

            fields[0] = new DBFField();
            fields[0].setName("TUCL");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setFieldLength(10);
            fields[0].setDecimalCount(4);

			fields[1] = new DBFField();
            fields[1].setName("MIN_ELEV");
            fields[1].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[1].setFieldLength(10);
            fields[1].setDecimalCount(4);

			fields[2] = new DBFField();
            fields[2].setName("MAX_ELEV");
            fields[2].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[2].setFieldLength(10);
            fields[2].setDecimalCount(4);
            
            fields[3] = new DBFField();
            fields[3].setName("OUTLET");
            fields[3].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[3].setFieldLength(6);
            fields[3].setDecimalCount(0);

            fields[4] = new DBFField();
            fields[4].setName("STRAHLER");
            fields[4].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[4].setFieldLength(6);
            fields[4].setDecimalCount(0);

            fields[5] = new DBFField();
            fields[5].setName("SHREVE");
            fields[5].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[5].setFieldLength(6);
            fields[5].setDecimalCount(0);

            fields[6] = new DBFField();
            fields[6].setName("DIST2MOUTH");
            fields[6].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[6].setFieldLength(10);
            fields[6].setDecimalCount(4);

            fields[7] = new DBFField();
            fields[7].setName("DS_NODES");
            fields[7].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[7].setFieldLength(6);
            fields[7].setDecimalCount(0);

            fields[8] = new DBFField();
            fields[8].setName("IS_OUTLET");
            fields[8].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[8].setFieldLength(1);
            fields[8].setDecimalCount(0);

            fields[9] = new DBFField();
            fields[9].setName("DS_LINK_ID");
            fields[9].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[9].setFieldLength(6);
            fields[9].setDecimalCount(0);

            fields[10] = new DBFField();
            fields[10].setName("MAINSTEM");
            fields[10].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[10].setFieldLength(1);
            fields[10].setDecimalCount(0);

            fields[11] = new DBFField();
            fields[11].setName("TRIB_ID");
            fields[11].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[11].setFieldLength(6);
            fields[11].setDecimalCount(0);
            
      		ShapeFile output = new ShapeFile(outputFile, ShapeType.POLYLINE, fields);

      		fields = new DBFField[1];

            fields[0] = new DBFField();
            fields[0].setName("FID");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setFieldLength(10);
            fields[0].setDecimalCount(0);

            ShapeFile outputConfluences = new ShapeFile(outputFile.replace(".shp", "_confluences.shp"), ShapeType.POINT, fields);
            whitebox.geospatialfiles.shapefile.Point pointOfInterest

			ShapeFile outputChannelHeads = new ShapeFile(outputFile.replace(".shp", "_channelHeads.shp"), ShapeType.POINT, fields);

            ShapeFile outputOutlets = new ShapeFile(outputFile.replace(".shp", "_outlets.shp"), ShapeType.POINT, fields);


            // first enter the line end-nodes into a kd-tree
			int numFeatures = input.getNumberOfRecords()
        	int count = 0;
			double[][] points;
			int[] partData;
			int startingPointInPart, endingPointInPart
			int i, numParts, numPoints, recNum, part, p
			int outletNum = 1;
			int featureNum = 0;
			int totalNumParts = 0;
			boolean isBeyondEdgeLine
			boolean isInterior
			boolean flag = true
			List<KdTree.Entry<Integer>> results;
            double[] entry;
			
            pluginHost.updateProgress("Pre-processing...", 0)
           	// count the number of parts
           	for (ShapeFileRecord record : input.records) {
				totalNumParts += record.getGeometry().getParts().length
            }

			boolean[] crossesNoData = new boolean[totalNumParts];
            linkMag = new double[totalNumParts];
			isBeyondEdge = new boolean[totalNumParts];

			StreamLinkKeyPoints[] linkKeyPoints = new StreamLinkKeyPoints[totalNumParts];
			
			linkLengths = new double[totalNumParts];
			outletNums = new int[totalNumParts];
			numDownstreamNodes = new int[totalNumParts];
			pointsTree = new KdTree.SqrEuclid<Integer>(2, null)

			double[] linkMinElev = new double[totalNumParts];
			double[] linkMaxElev = new double[totalNumParts];
			int[] downstreamLink = new int[totalNumParts];
			for (i = 0; i < totalNumParts; i++) {
				linkMinElev[i] = Double.POSITIVE_INFINITY
				linkMaxElev[i] = Double.NEGATIVE_INFINITY
				downstreamLink[i] = -99;
			}
			
            // Read the end-nodes into the KD-tree. 
            featureNum = -1;
			for (ShapeFileRecord record : input.records) {
				recNum = record.getRecordNumber()
                points = record.getGeometry().getPoints()
				numPoints = points.length;
				partData = record.getGeometry().getParts()
				numParts = partData.length
				for (part = 0; part < numParts; part++) {
					featureNum++
					startingPointInPart = partData[part];
                    if (part < numParts - 1) {
                        endingPointInPart = partData[part + 1];
                    } else {
                        endingPointInPart = numPoints;
                    }

					// Is this line off the edge of the DEM or within an 
					// area of nodata?
					isBeyondEdgeLine = true

					for (i = startingPointInPart; i < endingPointInPart; i++) {
                    
						row = dem.getRowFromYCoordinate(points[i][1]);
						col = dem.getColumnFromXCoordinate(points[i][0]);
						z = dem.getValue(row, col)
						if (i == startingPointInPart) { z1 = z; }
						if (i == endingPointInPart - 1) { z2 = z; }
						
						if (z != nodata) {
						    isBeyondEdgeLine = false;
						    if (z < linkMinElev[featureNum]) { linkMinElev[featureNum] = z}
						    if (z > linkMaxElev[featureNum]) { linkMaxElev[featureNum] = z}
						} else {
							crossesNoData[featureNum] = true;
						}
					}
                    
                    if (isBeyondEdgeLine) {
						isBeyondEdge[featureNum] = true;
                    } else {
                    	// calculate the length of this line
	                    length = 0;
	                    for (i = startingPointInPart + 1; i < endingPointInPart; i++) {
	                    	length += distMultiplier * Math.sqrt((points[i][0] - points[i - 1][0]) * (points[i][0] - points[i - 1][0]) + (points[i][1] - points[i - 1][1]) * (points[i][1] - points[i - 1][1]))
	                    }
	                    linkLengths[featureNum] = length;
                    }

                    // add both the end points to the kd-tree
                    x1 = points[startingPointInPart][0]
                	y1 = points[startingPointInPart][1]
                	entry = [y1, x1]
					pointsTree.addPoint(entry, new Integer(featureNum));

					x2 = points[endingPointInPart-1][0]
                	y2 = points[endingPointInPart-1][1]
                	entry = [y2, x2]
					pointsTree.addPoint(entry, new Integer(featureNum));

					linkKeyPoints[featureNum] = new StreamLinkKeyPoints(x1, y1, z1, x2, y2, z2);
					
				}

				progress = (int)(100f * recNum / numFeatures)
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Building search tree", progress)
            		oldProgress = progress
            		// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}
            }

            /*
             * Now we must find y-junctions. This occurs where
             * a stream link's end node intersects with another 
             * stream link but not at one of its end-nodes. Instead,
             * it touches one of its intermediate nodes. We will
             * perform a NN search at the location of all 
             * intermediate nodes and wherever one is within the 
             * search distance of an end-node (already in the kd-tree)
             * then it will be added to the kd-tree as well.
             */
			 int numYJunctions = 0;
			 featureNum = -1;
			for (ShapeFileRecord record : input.records) {
				recNum = record.getRecordNumber()
                points = record.getGeometry().getPoints()
				numPoints = points.length;
				partData = record.getGeometry().getParts()
				numParts = partData.length
				for (part = 0; part < numParts; part++) {
					featureNum++
					startingPointInPart = partData[part];
                    if (part < numParts - 1) {
                        endingPointInPart = partData[part + 1];
                    } else {
                        endingPointInPart = numPoints;
                    }
					for (i = startingPointInPart + 1; i <= (endingPointInPart - 1) - 1; i++) {
                    	x = points[i][0];
		            	y = points[i][1];
		            	entry = [y, x];
		            	results = pointsTree.neighborsWithinRange(entry, searchDist);
						
						if (results.size() > 0) {
							
							// add it to the tree
							entry = [y, x]
							pointsTree.addPoint(entry, new Integer(featureNum));
							numYJunctions++;

							linkKeyPoints[featureNum].addIntermediatePoint(x, y);

							pointOfInterest = new whitebox.geospatialfiles.shapefile.Point(x, y);                  
                			rowData = new Object[1];
                			rowData[0] = new Double(2);
                			outputConfluences.addRecord(pointOfInterest, rowData);
						}
					}
				}

				progress = (int)(100f * recNum / numFeatures)
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Building search tree", progress)
            		oldProgress = progress
            		// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}
            }

            /*
             *  Exterior links can be identified 
             *  as lines that either do not connect to another
             *  or that have at least one end-node with a NoData
             *  elevation value. Exterior links include both 
             *  channel heads (first-order stream) and outlet links.
             *  Add each of these to a priority queue.
             */

            PriorityQueue<StreamLink> queue = new PriorityQueue<StreamLink>(totalNumParts);
			
            boolean[] isChannelHead = new boolean[totalNumParts];
            boolean[] isExteriorLink = new boolean[totalNumParts];
            boolean isExterior
            boolean[] isOutletLink = new boolean[totalNumParts];
            int id
			oldProgress = -1
            for (i = 0; i < totalNumParts; i++) {
            	if (!isBeyondEdge[i]) {
            		z = Double.POSITIVE_INFINITY;
					/*
					 * To be an exterior link, it must have 
					 * at least one end that either isn't connected
					 * to any other link, has one link end that 
					 * is nodata in the DEM, or
					 */
					isExterior = false
					x = linkKeyPoints[i].endPoint1.x;
		            y = linkKeyPoints[i].endPoint1.y;
		            entry = [y, x];
		            results = pointsTree.neighborsWithinRange(entry, searchDist);
					j = 0;
					for (n = 0; n < results.size(); n++) {
						id = (int)results.get(n).value;
		            	if (id != i && isBeyondEdge[id] == false) {
		            		j++;
		            		if (linkMinElev[id] < z) { z = linkMinElev[id]; }
		            	}
		            }
	
		            if (j == 0) {
		            	isExterior = true
		            	x1 = x;
		            	y1 = y;
		            }
	
//		            if (!isExterior) {
		            	x = linkKeyPoints[i].endPoint2.x;
			            y = linkKeyPoints[i].endPoint2.y;
			            entry = [y, x];
			            results = pointsTree.neighborsWithinRange(entry, searchDist);
						j = 0;
						for (n = 0; n < results.size(); n++) {
							id = (int)results.get(n).value;
			            	if (id != i && !isBeyondEdge[id]) {
			            		j++;
			            		if (linkMinElev[id] < z) { z = linkMinElev[id]; }
			            	}
			            }
		
			            if (j == 0) {
			            	isExterior = true
			            	x1 = x;
		            		y1 = y;
			            }
//		            }
	
		            if (isExterior || crossesNoData[i]) {
		            	isExteriorLink[i] = true;
		            	if (linkMinElev[i] <= z || crossesNoData[i]) { 
//		            		isOutletLink[i] = true;

		            		z = linkMinElev[i];
							queue.add(new StreamLink(i, z + maxRidgeCuttingHeight));
		            	
		            	}
//						z = linkMinElev[i];
//						queue.add(new StreamLink(i, z + maxRidgeCuttingHeight));

//						if (isExterior) {
//				        	pointOfInterest = new whitebox.geospatialfiles.shapefile.Point(x1, y1);                  
//		                	rowData = new Object[1];
//		                	rowData[0] = new Double(1);
//		                	outputChannelHeads.addRecord(pointOfInterest, rowData);
//				        }
		            }
            	}
				progress = (int)(100f * i / totalNumParts)
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Finding starting points", progress)
            		oldProgress = progress
            		// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}
			}

			// perform the priority-flood operation
			int numSnappedOutlets = 0;
			StreamLink sl;
			boolean[] haveVisited = new boolean[totalNumParts];
			boolean[] haveEnteredQueue = new boolean[totalNumParts];
			int[] numInfowingLinks = new int[totalNumParts];
			double[] distToOutlet = new double[totalNumParts];
			int[] tribNum = new int[totalNumParts];
			int link;
			int currentMaxOutletNum = 0;
			int dsn;
			boolean isConfluence;
			int numLinks;
			int totalNumLinks;
			int numLinksVisited = 0;
			XYPoint endPoint;
			oldProgress = -1;
			while (!queue.isEmpty()) {
				sl = queue.poll();
				link = sl.index;
				if (!haveVisited[link]) {
					haveVisited[link] = true;
					haveEnteredQueue[link] = true;

					distToOutlet[link] += linkLengths[link];

					// What is the downstream link?
					dsn = downstreamLink[link];
	
					// What outlet number does the DSN belong to?
					if (dsn >= 0) {
						outletNum = outletNums[dsn];				
					} else {
//						// There isn't a DSN and we need a new outlet number
//						currentMaxOutletNum++;
//						outletNum = currentMaxOutletNum;
//						outletNums[link] = outletNum;
//						isOutletLink[link] = true;
						
						// which end point is the downstream outlet node?
						endPoint = linkKeyPoints[link].endPoint1;

						x = endPoint.x;
				        y = endPoint.y;
				        entry = [y, x];
				        results = pointsTree.neighborsWithinRange(entry, searchDist);
				        numLinks = 0;
				        for (i = 0; i < results.size(); i++) {
				        	n = (int)results.get(i).value;
				        	if (!isBeyondEdge[n] && !haveVisited[n] && !isOutletLink[n]) {
				        		numLinks++;
				        	}
				        }

				        if (numLinks > 0) {
				        	// end point 2 is the downstream node
				        	x = linkKeyPoints[link].endPoint2.x;
				        	y = linkKeyPoints[link].endPoint2.y;
				        } else {
				        	// how many linking nodes are at end point 2?
				        	endPoint = linkKeyPoints[link].endPoint2;

							x = endPoint.x;
					        y = endPoint.y;
					        entry = [y, x];
					        results = pointsTree.neighborsWithinRange(entry, searchDist);
					        numLinks = 0;
					        for (i = 0; i < results.size(); i++) {
					        	n = (int)results.get(i).value;
					        	if (!isBeyondEdge[n] && !haveVisited[n] && !isOutletLink[n]) {
					        		numLinks++;
					        	}
					        }

							if (numLinks > 0) {
								// end point 1 is the downstream node
								x = linkKeyPoints[link].endPoint1.x;
					        	y = linkKeyPoints[link].endPoint1.y;
							} else { // it's a single channel stream, which end is lower?
								if (linkKeyPoints[link].z1 < linkKeyPoints[link].z2 || 
								(linkKeyPoints[link].z1 == nodata && linkKeyPoints[link].z2 != nodata)) {
									x = linkKeyPoints[link].endPoint1.x;
						        	y = linkKeyPoints[link].endPoint1.y;
								} else {
									x = linkKeyPoints[link].endPoint2.x;
						        	y = linkKeyPoints[link].endPoint2.y;
								}
							}
				        }

				        if (!crossesNoData[link]) {
				        	/* This is a dangling stream. First let's make 
				        	 *  sure that there isn't a link end node from 
				        	 *  a previously discovered outlet nearby that 
				        	 *  we could connect to this outlet point
				        	*/ 
				        	entry = [y, x];
					        results = pointsTree.nearestNeighbor(entry, 3, false);
					        int snappedNeighbour = -1;
					        for (i = 0; i < results.size(); i++) {
				        		n = (int)results.get(i).value;
					        	if (!isBeyondEdge[n] && haveVisited[n] && isExteriorLink[n] && n != link) {
				    	    		// Check to see if the distance is less than the specified
				    	    		// snap distance.
				    	    		if (results.get(i).distance < snapDistance) {
				    	    			snappedNeighbour = n;
				    	    			break;
				    	    		}
				        		}
				        	}
				        	if (snappedNeighbour >= 0) {
				        		// we found a neighbour to snap to
				        		dsn = snappedNeighbour;
				        		outletNum = outletNums[dsn];
				        		outletNums[link] = outletNum;
				        		downstreamLink[link] = dsn;
				        		numInfowingLinks[dsn]++;
				        		numDownstreamNodes[link] = numDownstreamNodes[dsn] + 1;
				        		distToOutlet[link] += distToOutlet[dsn];
				        		numSnappedOutlets++;
				        	} else {
				        		// it is a true outlet

				        		// There isn't a DSN and we need a new outlet number
								currentMaxOutletNum++;
								outletNum = currentMaxOutletNum;
								outletNums[link] = outletNum;
								isOutletLink[link] = true;
								
				        		pointOfInterest = new whitebox.geospatialfiles.shapefile.Point(x, y);                  
		                		rowData = new Object[1];
		                		rowData[0] = new Double(outletNum);
		                		outputOutlets.addRecord(pointOfInterest, rowData);
				        	}
				        } else {
				        	// There isn't a DSN and we need a new outlet number
							currentMaxOutletNum++;
							outletNum = currentMaxOutletNum;
							outletNums[link] = outletNum;
							isOutletLink[link] = true;
						
				        	pointOfInterest = new whitebox.geospatialfiles.shapefile.Point(x, y);                  
	                		rowData = new Object[1];
	                		rowData[0] = new Double(outletNum);
	                		outputOutlets.addRecord(pointOfInterest, rowData);
				        }
					}

					for (XYPoint pt : linkKeyPoints[link].getAllPoints()) {
						x = pt.x;
				        y = pt.y;
				        entry = [y, x];
				        results = pointsTree.neighborsWithinRange(entry, searchDist);
				        numLinks = 0;
				        for (i = 0; i < results.size(); i++) {
				        	n = (int)results.get(i).value;
				        	if (!isBeyondEdge[n] && !haveEnteredQueue[n]) {
				        		numLinks++;
				        	}
				        }
				        isConfluence = (numLinks > 1) ? true : false; 
				        if (isConfluence) {
				        	pointOfInterest = new whitebox.geospatialfiles.shapefile.Point(x, y);                  
		                	rowData = new Object[1];
		                	rowData[0] = new Double(1);
		                	outputConfluences.addRecord(pointOfInterest, rowData);
				        }
				        for (i = 0; i < results.size(); i++) {
				        	n = (int)results.get(i).value;
				        	if (!isBeyondEdge[n] && !haveEnteredQueue[n]) {
				        		// add the link to the queue
				        		z = linkMinElev[n];
				            	queue.add(new StreamLink(n, z));

				            	haveEnteredQueue[n] = true;

				            	// update the DSN for this link
				            	downstreamLink[n] = link;
								if (isConfluence) {
				            		numDownstreamNodes[n] = numDownstreamNodes[link] + 1;
								} else {
									numDownstreamNodes[n] = numDownstreamNodes[link];
								}

								distToOutlet[n] += distToOutlet[link];
								
				            	outletNums[n] = outletNum;

				            	numInfowingLinks[link]++;
				        	}
				        }
					}

					numLinksVisited++;
	        		progress = (int)(100.0 * numLinksVisited / totalNumParts)
	        		if (progress != oldProgress) {
						pluginHost.updateProgress("Priority-Flood Operation:", progress)
	            		oldProgress = progress
	            		// check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
	            	}
				}
			}

			// calculate the link mag variables
			int[] strahlerOrder = new int[totalNumParts];
			int[] shreveOrder = new int[totalNumParts];
			LinkedList<Integer> stack = new LinkedList<Integer>();
			boolean foundDownstreamEnd;
			for (i = 0; i < totalNumParts; i++) {
				if (numInfowingLinks[i] == 0 && !isBeyondEdge[i]) {
					stack.push(i);
					strahlerOrder[i] = 1;
					shreveOrder[i] = 1;

					// this is a headwater, find which end is the channel head
					foundDownstreamEnd = false;
					//outletNum = outletNums[i];
					dsn = downstreamLink[i];
					endPoint = linkKeyPoints[i].endPoint1;
					x = endPoint.x;
			        y = endPoint.y;
			        entry = [y, x];
			        results = pointsTree.neighborsWithinRange(entry, searchDist);
			        for (j = 0; j < results.size(); j++) {
			        	n = (int)results.get(j).value;
			        	if (n == dsn) {
			        		foundDownstreamEnd = true;
			        	}
			        }
			        if (!foundDownstreamEnd) {
			        	pointOfInterest = new whitebox.geospatialfiles.shapefile.Point(x, y);                  
            	    	rowData = new Object[1];
                		rowData[0] = new Double(1);
                		outputChannelHeads.addRecord(pointOfInterest, rowData);
			        } else {
				        endPoint = linkKeyPoints[i].endPoint2;
						x = endPoint.x;
				        y = endPoint.y;
				        pointOfInterest = new whitebox.geospatialfiles.shapefile.Point(x, y);                  
            	    	rowData = new Object[1];
                		rowData[0] = new Double(1);
                		outputChannelHeads.addRecord(pointOfInterest, rowData);
			        }
				}
			}
			
			count = 0;
			oldProgress = -1;
			while (!stack.isEmpty()) {
				i = stack.pop();
				linkMag[i] += linkLengths[i];
				dsn = downstreamLink[i];
				if (dsn >= 0) {
					// pass this downstream
					linkMag[dsn] += linkMag[i];
					numInfowingLinks[dsn]--;
					if (numInfowingLinks[dsn] == 0) {
						stack.push(dsn);
					}

					if (strahlerOrder[dsn] == strahlerOrder[i]) {
						strahlerOrder[dsn]++;
					} else if (strahlerOrder[i] > strahlerOrder[dsn]) {
						strahlerOrder[dsn] = strahlerOrder[i];
					}

					shreveOrder[dsn] += shreveOrder[i];
				}
				count++;
				progress = (int)(100f * count / totalNumParts)
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Accumulation operations:", progress)
            		oldProgress = progress
            		// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}
			}

			// perform the outlet-to-head ops like finding the main stem
			// and assign tributary numbers
            boolean[] isMainStem = new boolean[totalNumParts];
            stack = new LinkedList<Integer>();
            int currentTribNum = 0;

            for (i = 0; i < totalNumParts; i++) {
				if (isOutletLink[i]) {
					isMainStem[i] = true;
					stack.push(i);
					currentTribNum++;
					tribNum[i] = currentTribNum;
				}
				progress = (int)(100f * i / totalNumParts)
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Assigning tributary ID:", progress)
            		oldProgress = progress
            		// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}
            }

			List<Integer> neighbourList = new ArrayList<Integer>();
			count = 0;
			oldProgress = -1;
			while (!stack.isEmpty()) {
				i = stack.pop();
				neighbourList.clear();
				double maxTUCL = 0;
				int maxTUCLlink = -1;
				for (XYPoint pt : linkKeyPoints[i].getAllPoints()) {
					x = pt.x;
			        y = pt.y;
			        entry = [y, x];
			        results = pointsTree.neighborsWithinRange(entry, searchDist);
			        numLinks = 0;
			        for (j = 0; j < results.size(); j++) {
			        	n = (int)results.get(j).value;
			        	if (downstreamLink[n] == i) { 
			        		neighbourList.add(n);
			        		if (linkMag[n] > maxTUCL) {
			        			maxTUCL = linkMag[n];
			        			maxTUCLlink = n;
			        		}
			        	}
			        }
				}
				if (maxTUCLlink >= 0) {
					//isMainStem[maxTUCLlink] = true;
					for (Integer q : neighbourList) {
						n = (int)q;
						// add it to the stack
						stack.push(n);
						if (n != maxTUCLlink) {
							currentTribNum++;
							tribNum[n] = currentTribNum;
						} else {
							tribNum[n] = tribNum[i];
							if (isMainStem[downstreamLink[n]]) {
								isMainStem[n] = true;
							}
						}
					}
				}
						
            	count++;
				progress = (int)(100f * count / totalNumParts)
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Assigning tributary ID:", progress)
            		oldProgress = progress
            		// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}
			}
            
				
			// Output the data into the attribute table.
			featureNum = -1;
			oldProgress = -1;
			count = 0;
			for (ShapeFileRecord record : input.records) {
				recNum = record.getRecordNumber()
                points = record.getGeometry().getPoints()
				numPoints = points.length;
				partData = record.getGeometry().getParts()
				numParts = partData.length
				for (part = 0; part < numParts; part++) {
					featureNum++
					if (!isBeyondEdge[featureNum]) {
		                rowData = new Object[12];
	            		rowData[0] = new Double(linkMag[featureNum]);
	            		rowData[1] = new Double(linkMinElev[featureNum]);
	                	rowData[2] = new Double(linkMaxElev[featureNum]);
	                	rowData[3] = new Double(outletNums[featureNum]);
	                	rowData[4] = new Double(strahlerOrder[featureNum]);
	                	rowData[5] = new Double(shreveOrder[featureNum]);
	                	rowData[6] = new Double(distToOutlet[featureNum]);
	                	rowData[7] = new Double(numDownstreamNodes[featureNum]);
	                	if (isOutletLink[featureNum]) {
	                		rowData[8] = new Double(1.0);
	                	} else {
	                		rowData[8] = new Double(0.0);
	                	}
	                	rowData[9] = new Double(downstreamLink[featureNum]);
	                	if (isMainStem[featureNum]) {
	                		rowData[10] = new Double(1.0);
	                	} else {
	                		rowData[10] = new Double(0.0);
	                	}
	                	rowData[11] = new Double(tribNum[featureNum]);
	                	Geometry poly;
	                	poly = new PolyLine(partData, points);
	                	output.addRecord(poly, rowData);
					}
				}

				count++
                progress = (int)(100f * count / numFeatures)
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Writing data:", progress)
            		oldProgress = progress
            		// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}
            }

			output.write();
            
			// display the output streams
			pluginHost.returnData(outputFile)

			outputConfluences.write();
            pluginHost.returnData(outputFile.replace(".shp", "_confluences.shp"));

            outputChannelHeads.write();
			pluginHost.returnData(outputFile.replace(".shp", "_channelHeads.shp"));

			outputOutlets.write();
			pluginHost.returnData(outputFile.replace(".shp", "_outlets.shp"));

			
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

    @CompileStatic
    class StreamLinkKeyPoints {
    	public XYPoint endPoint1;
    	public XYPoint endPoint2;
    	public double z1, z2;
    	public ArrayList<XYPoint> intermediatePoints = new ArrayList<>();

    	public StreamLinkKeyPoints(double x1, double y1, double z1, double x2, double y2, double z2) {
    		endPoint1 = new XYPoint(x1, y1);
    		this.z1 = z1;
    		endPoint2 = new XYPoint(x2, y2);
    		this.z2 = z2;
    	}

    	public void addIntermediatePoint(double x, double y) {
    		intermediatePoints.add(new XYPoint(x, y));
    	}

    	public ArrayList<XYPoint> getAllPoints() {
    		ArrayList<XYPoint> points = new ArrayList<>();
    		points.add(endPoint1);
    		points.add(endPoint2);
    		for (XYPoint p : intermediatePoints) {
    			points.add(p);
    		}

    		return points;
    	}
    }

    @CompileStatic
    class StreamLink implements Comparable<StreamLink> {
    	public int index;
    	public double minZ;

    	public StreamLink(int index, double minZ) {
    		this.index = index;
    		this.minZ = minZ;
    	}

    	@Override
        public int compareTo(StreamLink other) {
        	if (this.minZ < other.minZ) {
        		return -1;
        	} else if (this.minZ > other.minZ) {
        		return 1;
        	} else {
        		if (this.index < other.index) {
        			return -1;
        		} else if (this.index > other.index) {
        			return 1;
        		} else {
        			return 0;
        		}
        	}
        }
    }
	
	@CompileStatic
    class Node {

        public int featureNum;
        public int recNum;
        public boolean isActive = true;

        public Node(int featureNum, int partNum) {
            this.featureNum = featureNum;
            this.recNum = recNum;
        }
    }
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def tdf = new VectorStreamNetworkAnalysis(pluginHost, args, name, descriptiveName)
}
