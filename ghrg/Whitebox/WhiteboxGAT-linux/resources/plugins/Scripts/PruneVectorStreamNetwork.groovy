/*
 * Copyright (C) 2016 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "PruneVectorStreamNetwork"
def descriptiveName = "Prune Vector Stream Network"
def description = "Prunes a vector stream network"
def toolboxes = ["StreamAnalysis"]

public class PruneVectorStreamNetwork implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public PruneVectorStreamNetwork(WhiteboxPluginHost pluginHost, 
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
			sd.addDialogDataInput("Minimum link upstream channel length (xy-units)", "Min. Upstream Channel Length:", "", true, true)
            sd.addDialogCheckBox("Extend stream main stems to sources?", "Extend stream main stems to sources?", true)
			
			// resize the dialog to the standard size and display it
			sd.setSize(800, 400)
			sd.visible = true
		}
	}

	boolean[] isInactive;
    double[] linkMag;
	boolean[] isBeyondEdge;
	double[][] endNodeCoordinates;
	double[] linkLengths;
	int[] outletNums;
	int[] numDownstreamNodes;
	KdTree<Integer> pointsTree;

	@CompileStatic
	private void accumulate(int featureNum, int callingNum, int outletNum, int nodeNum) {
		if (pluginHost.isRequestForOperationCancelSet()) {
			pluginHost.showFeedback("Operation cancelled")
			return
		}
		int validLinks
		double x, y;
		int n, i;
        double sum = linkLengths[featureNum];
		double searchDist = 0.0000001;
		double[] entry;
        List<KdTree.Entry<Integer>> results;
        boolean searchLink;
        isInactive[featureNum] = true;
        numDownstreamNodes[featureNum] = nodeNum
		
		// find the adjoining active links
		x = endNodeCoordinates[featureNum][0];
        y = endNodeCoordinates[featureNum][1];
        entry = [y, x];
        results = pointsTree.neighborsWithinRange(entry, searchDist);
		searchLink = true;
		validLinks = 0
		for (i = 0; i < results.size(); i++) {
        	n = (int)results.get(i).value;
        	if (n == callingNum) {
        		searchLink = false; // don't search the end connected to the origin node
        	}
        	if (!isBeyondEdge[n] && !isInactive[n]) {
        		validLinks++
        	}
        }
        if (searchLink) {
        	if (validLinks > 1) { nodeNum++ }
			for (i = 0; i < results.size(); i++) {
            	n = (int)results.get(i).value;
            	if (!isBeyondEdge[n] && !isInactive[n]) {
            		isInactive[n] = true;
            		outletNums[n] = outletNum;
            		accumulate(n, featureNum, outletNum, nodeNum);
            		sum += linkMag[n];
            	}
            }
        }

        x = endNodeCoordinates[featureNum][2];
        y = endNodeCoordinates[featureNum][3];
        entry = [y, x];
        results = pointsTree.neighborsWithinRange(entry, searchDist);
		searchLink = true;
		validLinks = 0;
		for (i = 0; i < results.size(); i++) {
        	n = (int)results.get(i).value;
        	if (n == callingNum) {
        		searchLink = false; // don't search the end connected to the origin node
        	}
        	if (!isBeyondEdge[n] && !isInactive[n]) {
        		validLinks++
        	}
        }
        if (searchLink) {
        	if (validLinks > 1) { nodeNum++ }
			for (i = 0; i < results.size(); i++) {
            	n = (int)results.get(i).value;
            	if (!isBeyondEdge[n] && !isInactive[n]) {
            		isInactive[n] = true;
            		outletNums[n] = outletNum;
            		accumulate(n, featureNum, outletNum, nodeNum);
            		sum += linkMag[n];
            	}
            }
        }
        
        linkMag[featureNum] = sum;
	}

	@CompileStatic
	private void extendMainStem(int featureNum, int callingNum) {
		if (pluginHost.isRequestForOperationCancelSet()) {
			pluginHost.showFeedback("Operation cancelled")
			return
		}
		double x, y;
		int n, i;
        double largestInflowingLinkVal;
        int largestInflowingLinkID;
		double searchDist = 0.0000001;
		double[] entry;
        List<KdTree.Entry<Integer>> results;
        boolean searchLink;
        
		isInactive[featureNum] = true;
		
		// find the adjoining active links
		x = endNodeCoordinates[featureNum][0];
        y = endNodeCoordinates[featureNum][1];
        entry = [y, x];
        results = pointsTree.neighborsWithinRange(entry, searchDist);
		searchLink = true;
		for (i = 0; i < results.size(); i++) {
        	if ((int)results.get(i).value == callingNum) {
        		searchLink = false; // don't search the end connected to the origin node
        	}
        }
        if (searchLink) {
        	// find the inflowing neighbour with the largest link mag
        	largestInflowingLinkVal = 0;
        	largestInflowingLinkID = -1;
        	for (i = 0; i < results.size(); i++) {
            	n = (int)results.get(i).value;
            	if (!isBeyondEdge[n] && !isInactive[n]) {
            		if (linkMag[n] > largestInflowingLinkVal) {
            			largestInflowingLinkID = n;
            			largestInflowingLinkVal = linkMag[n];
            		}
            	}
            }

            linkMag[largestInflowingLinkID] = linkMag[featureNum];
            
			for (i = 0; i < results.size(); i++) {
            	n = (int)results.get(i).value;
            	if (!isBeyondEdge[n] && !isInactive[n]) {
            		isInactive[n] = true;
            		extendMainStem(n, featureNum);
            	}
            }
        }

        x = endNodeCoordinates[featureNum][2];
        y = endNodeCoordinates[featureNum][3];
        entry = [y, x];
        results = pointsTree.neighborsWithinRange(entry, searchDist);
		searchLink = true;
		for (i = 0; i < results.size(); i++) {
        	if ((int)results.get(i).value == callingNum) {
        		searchLink = false; // don't search the end connected to the origin node
        	}
        }
        if (searchLink) {
			// find the inflowing neighbour with the larges link mag
        	largestInflowingLinkVal = 0;
        	largestInflowingLinkID = -1;
        	for (i = 0; i < results.size(); i++) {
            	n = (int)results.get(i).value;
            	if (!isBeyondEdge[n] && !isInactive[n]) {
            		if (linkMag[n] > largestInflowingLinkVal) {
            			largestInflowingLinkID = n;
            			largestInflowingLinkVal = linkMag[n];
            		}
            	}
            }

            linkMag[largestInflowingLinkID] = linkMag[featureNum];
            
			for (i = 0; i < results.size(); i++) {
            	n = (int)results.get(i).value;
            	if (!isBeyondEdge[n] && !isInactive[n]) {
            		isInactive[n] = true;
            		extendMainStem(n, featureNum);
            	}
            }
        }
	}
	
	@CompileStatic
	private void execute(String[] args) {
		try {
			/* 
			 *  The following code calculates the upstream
			 *  channel length for each link in the  vector
			 *  stream network. This task involves identifying
			 *  exterior stream links. These include both 
			 *  first-order and outlet links, on opposite
			 *  sides of the network. Each exterior link is
			 *  visited and a recursive accumulation algorithm
			 *  is applied to find all links connected to 
			 *  each exterior link. Notice exterior links 
			 *  are visited in order from lowest (with elevation
			 *  derived from an underlying DEM) to highest. 
			 *  Since outlet links will necessarily be lower
			 *  in elevation than all of the upstream 
			 *  first-order links, the algorithm implicitly 
			 *  distinguishes between outlet and first-order
			 *  streams. Link connections are identified by 
			 *  placing each link's end nodes coordinates 
			 *  into a kd-tree and performing nearest neighbour
			 *  searches. Also, all links that are either 
			 *  beyond the edges of the DEM or within areas of
			 *  nodata within the DEM are excluded from the
			 *  analysis. Also note that the DEM that is used
			 *  does not have to be particularly accurate, i.e.
			 *  it does not have to represent the entire stream
			 *  network well, since it is only being used 1) 
			 *  to distinguish between outlet and first-order 
			 *  exterior links, and 2) to exclude links beyond
			 *  the area of interest/within nodata areas.
			 *  
			 *  Once the upstream channel length (UCL) is measured,
			 *  the vector network can be pruned by parsing
			 *  links that have a UCL lower than a user-specified
			 *  threshold.
			 */
			 
	  		int progress, oldProgress, col, row;
	  		int n, j;
	  		double x, y, z, z1, z2;
	  		double searchDist = 0.0000001 // This has to be a small non-zero value and is used in the nearest-neighbour search.
	  		double length;
  		  	double distMultiplier = 1.0;
  		  	
	  		if (args.length < 5) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String streamsFile = args[0];
			String baseFile = args[1];
			String outputFile = args[2];
			double minLength = 0.0;
			if (!(args[3].toLowerCase().contains("not specified"))) {
				minLength = Double.parseDouble(args[3]);
			}
			boolean findMainStem = true;
			if (args[4].toLowerCase().contains("f")) {
				findMainStem = false;
			}
			
			// read the input image
			WhiteboxRaster baseRaster = new WhiteboxRaster(baseFile, "r");
			baseRaster.setForceAllDataInMemory(true)
			double nodata = baseRaster.getNoDataValue();
			int rows = baseRaster.getNumberRows()
			int cols = baseRaster.getNumberColumns()
			int rowsLessOne = rows - 1
			int colsLessOne = cols - 1

			if (baseRaster.getXYUnits().toLowerCase().contains("deg")) {
				double midLat = (baseRaster.getNorth() - baseRaster.getSouth()) / 2.0;
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
            DBFField[] fields = new DBFField[5];

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
            fields[4].setName("DS_NODES");
            fields[4].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[4].setFieldLength(6);
            fields[4].setDecimalCount(0);
            
      		ShapeFile output = new ShapeFile(outputFile, ShapeType.POLYLINE, fields);
			
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
				points = record.getGeometry().getPoints()
				totalNumParts += record.getGeometry().getParts().length
            }

			boolean[] crossesNoData = new boolean[totalNumParts];
            isInactive = new boolean[totalNumParts];
            linkMag = new double[totalNumParts];
			isBeyondEdge = new boolean[totalNumParts];
			endNodeCoordinates = new double[totalNumParts][4];
			linkLengths = new double[totalNumParts];
			outletNums = new int[totalNumParts];
			numDownstreamNodes = new int[totalNumParts];
			pointsTree = new KdTree.SqrEuclid<Integer>(2, new Integer(totalNumParts * 2))

			double[] linkMinElev = new double[totalNumParts];
			double[] linkMaxElev = new double[totalNumParts];
			for (i = 0; i < totalNumParts; i++) {
				linkMinElev[i] = Double.POSITIVE_INFINITY
				linkMaxElev[i] = Double.NEGATIVE_INFINITY
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
                    
						row = baseRaster.getRowFromYCoordinate(points[i][1]);
						col = baseRaster.getColumnFromXCoordinate(points[i][0]);
						z = baseRaster.getValue(row, col)
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
                    x = points[startingPointInPart][0]
                	y = points[startingPointInPart][1]
                	entry = [y, x]
					pointsTree.addPoint(entry, new Integer(featureNum));
					endNodeCoordinates[featureNum][0] = x;
					endNodeCoordinates[featureNum][1] = y;
					
					x = points[endingPointInPart-1][0]
                	y = points[endingPointInPart-1][1]
                	entry = [y, x]
					pointsTree.addPoint(entry, new Integer(featureNum));
					endNodeCoordinates[featureNum][2] = x;
					endNodeCoordinates[featureNum][3] = y;
					
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
             *  as lines that either do no connect to another
             *  or that have at least one end-node with a NoData
             *  elevation value. Exterior links include both 
             *  channel heads (first-order stream) and outlet links.
             *  Add each of these to a list, which will be sorted
             *  by elevation at the end.
             */
            List<IntegerDoublePair> exteriorLinks = new ArrayList<>();

      		//boolean[] isExteriorLink = new boolean[totalNumParts];
            boolean isExterior
            int id
			oldProgress = -1
            for (i = 0; i < totalNumParts; i++) {
            	if (!isBeyondEdge[i]) {
					/*
					 * To be an exterior link, it must have 
					 * at least one end that either isn't connected
					 * to any other link, has one link end that 
					 * is nodata in the DEM, or
					 */
					isExterior = false
					x = endNodeCoordinates[i][0];
		            y = endNodeCoordinates[i][1];
		            entry = [y, x];
		            results = pointsTree.neighborsWithinRange(entry, searchDist);
					j = 0;
					for (n = 0; n < results.size(); n++) {
						id = (int)results.get(n).value;
		            	if (id != i && isBeyondEdge[id] == false) {
		            		j++;
		            	}
		            }
	
		            if (j == 0) {
		            	isExterior = true
		            }
	
		            if (!isExterior) {
		            	x = endNodeCoordinates[i][2];
			            y = endNodeCoordinates[i][3];
			            entry = [y, x];
			            results = pointsTree.neighborsWithinRange(entry, searchDist);
						j = 0;
						for (n = 0; n < results.size(); n++) {
							id = (int)results.get(n).value;
			            	if (id != i && isBeyondEdge[id] == false) {
			            		j++;
			            	}
			            }
		
			            if (j == 0) {
			            	isExterior = true
			            }
		            }
	
		            if (isExterior) {
						z = linkMinElev[i];
		            	exteriorLinks.add(new IntegerDoublePair(i, z));
		            } else if (crossesNoData[i]) {
		            	z = linkMinElev[i];
		            	exteriorLinks.add(new IntegerDoublePair(i, z));
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

			// sort the list
			Collections.sort(exteriorLinks);

			// now visit each outlet link and accumulate
			j = 0;
			n = exteriorLinks.size() - 1;
			for (IntegerDoublePair link : exteriorLinks) {
				i = link.intValue;
				if (!isInactive[i]) { // && crossesNoData[i]) {
					outletNums[i] = outletNum;
					accumulate(i, -1, outletNum, 0);
					outletNum++;
				}
				progress = (int)(100f * j / n)
				j++;
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Accumulating...", progress)
            		oldProgress = progress
            		// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}
			}

			if (findMainStem) {
				isInactive = new boolean[totalNumParts];

				// now visit each outlet link and propagate it's value to its source
				j = 0;
				n = exteriorLinks.size() - 1;
				for (IntegerDoublePair link : exteriorLinks) {
					i = link.intValue;
					if (!isInactive[i]) { // && crossesNoData[i]) {
						extendMainStem(i, -1);
					}
					progress = (int)(100f * j / n)
					j++;
	            	if (progress != oldProgress) {
						pluginHost.updateProgress("Finding main stems...", progress)
	            		oldProgress = progress
	            		// check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
	            	}
				}
			}

			// Output the data into the attribute table.
			featureNum = -1;
			oldProgress = -1;
			for (ShapeFileRecord record : input.records) {
				recNum = record.getRecordNumber()
                points = record.getGeometry().getPoints()
				numPoints = points.length;
				partData = record.getGeometry().getParts()
				numParts = partData.length
				for (part = 0; part < numParts; part++) {
					featureNum++
	                if (linkMag[featureNum] >= minLength) {
	                	Object[] rowData = new Object[5];
	            		rowData[0] = new Double(linkMag[featureNum]);
	            		rowData[1] = new Double(linkMinElev[featureNum]);
	                	rowData[2] = new Double(linkMaxElev[featureNum]);
	                	rowData[3] = new Double(outletNums[featureNum]);
	                	rowData[4] = new Double(numDownstreamNodes[featureNum]);
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
	class IntegerDoublePair implements Comparable<IntegerDoublePair> {
		public int intValue;
		public double doubleValue;

		public IntegerDoublePair(int intValue, double doubleValue) {
			this.intValue = intValue;
			this.doubleValue = doubleValue;
		}

		@Override
        public int compareTo(IntegerDoublePair other) {
        	if (this.doubleValue < other.doubleValue) {
        		return -1;
        	} else if (this.doubleValue > other.doubleValue) {
        		return 1;
        	} else {
        		if (this.intValue < other.intValue) {
        			return -1;
        		} else if (this.intValue > other.intValue) {
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
	def tdf = new PruneVectorStreamNetwork(pluginHost, args, name, descriptiveName)
}
