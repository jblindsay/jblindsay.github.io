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
import whitebox.utilities.StringUtilities
import whitebox.structures.BooleanBitArray2D
import whitebox.structures.NibbleArray2D
import whitebox.structures.DoubleArray2D
import whitebox.structures.IntArray2D
import whitebox.structures.BoundingBox;
import whitebox.structures.KdTree;
import groovy.transform.CompileStatic
import groovy.time.TimeDuration
import groovy.time.TimeCategory

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "TopologicalBreachBurn"
def descriptiveName = "Topological Breach-Burn (Stream Burning)"
def description = "Burns streams into a DEM."
def toolboxes = ["DEMPreprocessing"]

public class TopologicalBreachBurn implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public TopologicalBreachBurn(WhiteboxPluginHost pluginHost, 
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
			sd.addDialogFile("Input DEM raster", "Input DEM Raster:", "open", "Raster Files (*.dep), DEP", true, false)
            sd.addDialogFile("Input streams file", "Input Streams File:", "open", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogFile("Output file", "Output File:", "save", "Raster Files (*.dep), DEP", true, false)
			
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
	KdTree<Integer> pointsTree;

	@CompileStatic
	private void accumulate(int featureNum, int callingNum, int outletNum) {
		if (pluginHost.isRequestForOperationCancelSet()) {
			pluginHost.showFeedback("Operation cancelled")
			return
		}
		double x, y;
		int n, i;
        double sum = linkLengths[featureNum];
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
			for (i = 0; i < results.size(); i++) {
            	n = (int)results.get(i).value;
            	if (!isBeyondEdge[n] && !isInactive[n]) {
            		isInactive[n] = true;
            		outletNums[n] = outletNum;
            		accumulate(n, featureNum, outletNum);
            		sum += linkMag[n];
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
			for (i = 0; i < results.size(); i++) {
            	n = (int)results.get(i).value;
            	if (!isBeyondEdge[n] && !isInactive[n]) {
            		isInactive[n] = true;
            		outletNums[n] = outletNum;
            		accumulate(n, featureNum, outletNum);
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
			Date start = new Date()
			
	  		int progress, oldProgress, col, row, colN, rowN, numPits, r, c
	  		int numSolvedCells = 0
	  		int dir, n, i, j;
	  		double x, y, z, z1, z2, zN, zTest, zN2, lowestNeighbour;
	  		boolean isPit, isEdgeCell, isStream
	  		GridCell gc
	  		double LARGE_NUM = Float.MAX_VALUE
			int numInflowing
  		  	double s, sN
  		  	double searchDist = 0.0000001 // This has to be a small non-zero value and is used in the nearest-neighbour search.
	  		double length;
			int numParts, numPoints, recNum, part, p
			int featureNum = 0;
			int totalNumParts = 0;
			boolean isEdgeLine
			boolean isInterior
			boolean flag = true
			double[][] points;
			int[] partData;
			int outletNum = 1;
			int startingPointInPart, endingPointInPart
			List<KdTree.Entry<Integer>> results;
            double[] entry;
			double x1, y1, x2, y2, xPrime, yPrime, d;
			BoundingBox box;
			int topRow, bottomRow, leftCol, rightCol;
			double rowYCoord, colXCoord;

	  		/*
	  		 * 7  8  1
	  		 * 6  X  2
	  		 * 5  4  3
	  		 */
			int[] dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ]
			int[] dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ]
			int[] backLink = [5, 6, 7, 8, 1, 2, 3, 4]
			double[] outPointer = [0, 1, 2, 4, 8, 16, 32, 64, 128]
			
			if (args.length < 3) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String demFile = args[0]
			String streamsFile = args[1]
			String outputFile = args[2]
			
			// read the input image
			pluginHost.updateProgress("Reading data...", 0);
			WhiteboxRaster dem = new WhiteboxRaster(demFile, "r")
			dem.setForceAllDataInMemory(true);
			double nodata = dem.getNoDataValue()
			int rows = dem.getNumberRows()
			int cols = dem.getNumberColumns()
			int rowsLessOne = rows - 1
			int colsLessOne = cols - 1
			int numCellsTotal = rows * cols
			def paletteName = dem.getPreferredPalette()

			double distMultiplier = 1.0;
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
                    //distMultiplier = 113200 * Math.cos(Math.toRadians(midLat));
                }
			}


			//double minVal = dem.getMinimumValue()
			int elevDigits = (String.valueOf((int)(dem.getMaximumValue()))).length()
			double elevMultiplier = Math.pow(10, 8-elevDigits)
			double SMALL_NUM = 1 / elevMultiplier * 10
			long priority

			ShapeFile input = new ShapeFile(streamsFile)
			ShapeType shapeType = input.getShapeType()
            if (shapeType != ShapeType.POLYLINE) {
            	pluginHost.showFeedback("The input shapefile should be of a POLYLINE ShapeType.")
            	return
            }
            

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
			
	  		// first enter the line end-nodes into a kd-tree
			int numFeatures = input.getNumberOfRecords()
        	int count = 0;
			
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
			pointsTree = new KdTree.SqrEuclid<Integer>(2, new Integer(totalNumParts * 2))

			double[] linkMinElev = new double[totalNumParts];
			for (i = 0; i < totalNumParts; i++) {
				linkMinElev[i] = Double.POSITIVE_INFINITY
			}
			
            // Read the end-nodes into the KD-tree. 
            pluginHost.updateProgress("Building search tree", 0)
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
					isEdgeLine = true

					for (i = startingPointInPart; i < endingPointInPart; i++) {
                    
						row = dem.getRowFromYCoordinate(points[i][1]);
						col = dem.getColumnFromXCoordinate(points[i][0]);
						z = dem.getValue(row, col)
						if (z != nodata) {
						    isEdgeLine = false;
						    if (z < linkMinElev[featureNum]) { linkMinElev[featureNum] = z}
						} else {
							crossesNoData[featureNum] = true;
						}
					}
                    
                    if (isEdgeLine) {
						isBeyondEdge[featureNum] = true;
                    } else {
                    	// calculate the length of this line
	                    length = 0;
	                    for (i = startingPointInPart + 1; i < endingPointInPart; i++) {
	                    	length += Math.sqrt((points[i][0] - points[i - 1][0]) * (points[i][0] - points[i - 1][0]) + (points[i][1] - points[i - 1][1]) * (points[i][1] - points[i - 1][1]))
	                    }
	                    linkLengths[featureNum] = length * distMultiplier;
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

      		boolean isExterior
            int id, idN
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
						exteriorLinks.add(new IntegerDoublePair(i, linkMinElev[i]));
		            } else if (crossesNoData[i]) {
		            	exteriorLinks.add(new IntegerDoublePair(i, linkMinElev[i]));
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
				if (!isInactive[i]) {
					outletNums[i] = outletNum;
					accumulate(i, -1, outletNum);
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

  		  	// now visit each outlet link and propagate it's value to its source
			isInactive = new boolean[totalNumParts];
			j = 0;
			n = exteriorLinks.size() - 1;
			for (IntegerDoublePair link : exteriorLinks) {
				i = link.intValue;
				if (!isInactive[i]) {
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

			// release some of the data from memory
			crossesNoData = null;
            isInactive = null;
            endNodeCoordinates = null;
			linkLengths = null;
			outletNums = new int[totalNumParts];
			pointsTree = null;
			

			// Now rasterize the stream lines, keeping track of collisions.
			IntArray2D linkID = new IntArray2D(rows, cols, -32768)
			DoubleArray2D linkPosition = new DoubleArray2D(rows, cols, -32768.0)
			BooleanBitArray2D linkEndNodes = new BooleanBitArray2D(rows, cols)
			BooleanBitArray2D isCollision = new BooleanBitArray2D(rows, cols)

			IntArray2D linkOutlet = new IntArray2D(rows, cols, -32768)
			
			int numStreamCells = 0;
			List<CollisionCell> collsionsList = new ArrayList<CollisionCell>();
			
			featureNum = 0;

			oldProgress = -1
			for (ShapeFileRecord record : input.records) {
				recNum = record.getRecordNumber()
                points = record.getGeometry().getPoints()
				numPoints = points.length;
				partData = record.getGeometry().getParts()
				numParts = partData.length
				for (part = 0; part < numParts; part++) {
					featureNum++
					if (!isBeyondEdge[featureNum - 1]) {
						box = new BoundingBox();             
						startingPointInPart = partData[part];
	                    if (part < numParts - 1) {
	                        endingPointInPart = partData[part + 1];
	                    } else {
	                        endingPointInPart = numPoints;
	                    }
	
						row = dem.getRowFromYCoordinate(points[startingPointInPart][1])
						col = dem.getColumnFromXCoordinate(points[startingPointInPart][0]);
	                    if (linkID.getValue(row, col) == -32768) {
	                        linkEndNodes.setValue(row, col, true);
	                        linkID.setValue(row, col, featureNum);
	                        linkOutlet.setValue(row, col, outletNums[featureNum - 1])
	                        numStreamCells++;
						}
	                    
	                    row = dem.getRowFromYCoordinate(points[endingPointInPart-1][1])
						col = dem.getColumnFromXCoordinate(points[endingPointInPart-1][0]);
	                    if (linkID.getValue(row, col) == -32768) {
	                        linkEndNodes.setValue(row, col, true);
	                        linkID.setValue(row, col, featureNum);
	                        linkOutlet.setValue(row, col, outletNums[featureNum - 1])
	                        numStreamCells++;
						}
	                    
	                    double[] distances = new double[endingPointInPart - startingPointInPart + 1];
	                    for (i = startingPointInPart + 1; i < endingPointInPart; i++) {
	                    	distances[i] = distances[i - 1] +  Math.sqrt((points[i][0] - points[i - 1][0]) * (points[i][0] - points[i - 1][0]) + (points[i][1] - points[i - 1][1]) * (points[i][1] - points[i - 1][1]))
	                    }
	                    
	                    for (i = startingPointInPart; i < endingPointInPart; i++) {
	                        if (points[i][0] < box.getMinX()) {
	                            box.setMinX(points[i][0]);
	                        }
	                        if (points[i][0] > box.getMaxX()) {
	                            box.setMaxX(points[i][0]);
	                        }
	                        if (points[i][1] < box.getMinY()) {
	                            box.setMinY(points[i][1]);
	                        }
	                        if (points[i][1] > box.getMaxY()) {
	                            box.setMaxY(points[i][1]);
	                        }
	                    }
	                    topRow = dem.getRowFromYCoordinate(box.getMaxY());
	                    bottomRow = dem.getRowFromYCoordinate(box.getMinY());
	                    leftCol = dem.getColumnFromXCoordinate(box.getMinX());
	                    rightCol = dem.getColumnFromXCoordinate(box.getMaxX());
	
						// find each intersection with a row.
	                    for (row = topRow; row <= bottomRow; row++) {
	
	                        rowYCoord = dem.getYCoordinateFromRow(row);
	                        // find the x-coordinates of each of the line segments 
	                        // that intersect this row's y coordinate
	
	                        for (i = startingPointInPart; i < endingPointInPart - 1; i++) {
	                            if (isBetween(rowYCoord, points[i][1], points[i + 1][1])) {
	                                y1 = points[i][1];
	                                y2 = points[i + 1][1];
	                                if (y2 != y1) {
	                                    x1 = points[i][0];
	                                    x2 = points[i + 1][0];
	
	                                    // calculate the intersection point
	                                    xPrime = x1 + (rowYCoord - y1) / (y2 - y1) * (x2 - x1);
	                                    col = dem.getColumnFromXCoordinate(xPrime);
	                                    
										j = linkID.getValue(row, col);

										z = dem.getValue(row, col)
										if (j == -32768 && z != nodata) {
		                                    numStreamCells++;
										} else if ((j != featureNum && j != -32768) && 
											!isCollision.getValue(row, col) &&
											!linkEndNodes.getValue(row, col)) {
//												if (j == -32768 || featureNum == -32768) {
//													println("$j $featureNum")
//												}
											collsionsList.add(new CollisionCell(linkMag[j - 1], linkMag[featureNum - 1]));
											isCollision.setValue(row, col, true);
										}
										
										if ((j == -32768 || linkMag[j - 1] < linkMag[featureNum - 1]) && z != nodata) {
		                                    linkID.setValue(row, col, featureNum);
		                                    linkOutlet.setValue(row, col, outletNums[featureNum - 1])
		                                    d = distances[i + 1] - Math.sqrt((xPrime - x2)*(xPrime - x2) + (rowYCoord - y2)*(rowYCoord - y2))
		                                    linkPosition.setValue(row, col, d);
										}
	                                }
	                            }
	                        }
	                    }
	
	                    // find each intersection with a column.
	                    for (col = leftCol; col <= rightCol; col++) {
	                        colXCoord = dem.getXCoordinateFromColumn(col);
	                        for (i = startingPointInPart; i < endingPointInPart - 1; i++) {
	                            if (isBetween(colXCoord, points[i][0], points[i + 1][0])) {
	                                x1 = points[i][0];
	                                x2 = points[i + 1][0];
	                                if (x1 != x2) {
	                                    y1 = points[i][1];
	                                    y2 = points[i + 1][1];
	
	                                    // calculate the intersection point
	                                    yPrime = y1 + (colXCoord - x1) / (x2 - x1) * (y2 - y1);
	                                    row = dem.getRowFromYCoordinate(yPrime);
	
	                                    j = linkID.getValue(row, col);
	                                    z = dem.getValue(row, col)
										if (j == -32768 && z != nodata) {
		                                    numStreamCells++;
										} else if (j != featureNum && 
										           !isCollision.getValue(row, col) &&
										           !linkEndNodes.getValue(row, col)) {
											collsionsList.add(new CollisionCell(linkMag[j - 1], linkMag[featureNum - 1]));
											isCollision.setValue(row, col, true);
										}
										
										if ((j == -32768 || linkMag[j - 1] < linkMag[featureNum - 1]) && z != nodata) {
		                                    linkID.setValue(row, col, featureNum);
		                                    linkOutlet.setValue(row, col, outletNums[featureNum - 1])
		                                    d = distances[i + 1] - Math.sqrt((yPrime - y2)*(yPrime - y2) + (colXCoord - x2)*(colXCoord - x2))
		                                    linkPosition.setValue(row, col, d);
										}
	                                }
	                            }
	                        }
	                    }
					}
				}
				
				count++
                progress = (int)(100f * count / numFeatures)
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Rasterizing Streams...", progress)
            		oldProgress = progress
            		// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
            	}
            }


			int numLinkCollisions = collsionsList.size();
			if (numLinkCollisions > 5) {
				/*
				 * Finding the threshold TUCL value used to prune the 
				 * network is an optimization problem. In choosing 
				 * a TUCL threshold, we would like to minimize the
				 * percentage of collsions remaining after pruning 
				 * and to minimize the percent of the network that
				 * is pruned. We simply have to find the minimum
				 * value in the combined function of the two 
				 * criteria.
				 */
	           	int numBins = 2000;
	           	int bin;
	           	double maxLinkLength = -1.0;
	           	for (i = 0; i < totalNumParts; i++) {
	           		z = Math.log(linkMag[i] + 1)
	           		if (z > maxLinkLength) { maxLinkLength = z; }
	           	}
	           	maxLinkLength += 0.1;
	           	double binSize = maxLinkLength / numBins;
	           	
	           	// calculate the percent of the network that is pruned
	           	// at a specified TUCL
	           	double[] percentPruned = new double[numBins];
	           	for (i = 0; i < totalNumParts; i++) {
	           		bin = (int)Math.floor(Math.log(linkMag[i] + 1) / binSize);
	           		percentPruned[bin]++;
	           	}
				
				percentPruned[0] = percentPruned[0] / totalNumParts;
	           	for (bin = 1; bin < numBins; bin++) {
	           		percentPruned[bin] = percentPruned[bin-1] + percentPruned[bin] / totalNumParts;
	           	}


	           	// calculate the percent of collsions remaining 
	           	// after pruning to a specified TUCL
	           	double[] collisionsRemaining = new double[numBins];
				for (i = 0; i < numLinkCollisions; i++) {
					bin = (int)Math.floor(Math.log(collsionsList.get(i).tucl + 1) / binSize);
	           		collisionsRemaining[bin]++;
				}
				
				collisionsRemaining[numBins-1] = collisionsRemaining[numBins-1] / numLinkCollisions;
	           	for (bin = numBins-2; bin >= 0; bin--) {
	           		collisionsRemaining[bin] = collisionsRemaining[bin+1] + collisionsRemaining[bin] / numLinkCollisions;
	           	}

	           	// add the two functions and find the min value 
	           	// and the TUCL at the min
	           	double minValue = 99.0;
				double thresholdTUCL;
				int minBin
				double collisionsWeight = 0.65;
				double pruningWeight = 1.0 - collisionsWeight;
				for (bin = 0; bin < numBins; bin++) {
	           		z = percentPruned[bin] * pruningWeight + collisionsRemaining[bin] * collisionsWeight;
	           		if (z < minValue) {
	           			minValue = z;
	           			thresholdTUCL = Math.exp(bin * binSize) - 1;
	           			minBin = bin;
	           		}
				}
				
				// outputs
				String str2 = "";
				str2 += "Total Num. Collisions:\t$numLinkCollisions\n"
				str2 += "Threshold TUCL:\t$thresholdTUCL\n"
				//str2 += "Optimization Function Value at Threshold:\t$minValue\n"
				str2 += "Num. Collisions Remaining:\t${(int)(numLinkCollisions * collisionsRemaining[minBin])} (${(numLinkCollisions * collisionsRemaining[minBin]) / numLinkCollisions * 100.0}%)\n"
				str2 += "% of Network Links Pruned:\t${percentPruned[minBin] * 100}\n"
				pluginHost.returnData(str2);

				// time to reclaim some of the memory
				collsionsList.clear();
				collisionsRemaining = new double[0];
				percentPruned = new double[0];
				
				
				// create the output streams file
	            DBFField[] fields = new DBFField[3];

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
	            fields[2].setName("OUTLET");
	            fields[2].setDataType(DBFField.DBFDataType.NUMERIC);
	            fields[2].setFieldLength(6);
	            fields[2].setDecimalCount(0);

				String outputStreamsFile = outputFile.replace(".dep", "_strms.shp") //streamsFile.replace(".shp", "_pruned.shp")
	      		ShapeFile outputStreams = new ShapeFile(outputStreamsFile, ShapeType.POLYLINE, fields);

				// Output the pruned streams data
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
		                if (linkMag[featureNum] >= thresholdTUCL) {
		                	Object[] rowData = new Object[3];
		            		rowData[0] = new Double(linkMag[featureNum]);
		            		rowData[1] = new Double(linkMinElev[featureNum]);
		                	rowData[2] = new Double(outletNums[featureNum]);
		                	Geometry poly;
		                	poly = new PolyLine(partData, points);
		                	outputStreams.addRecord(poly, rowData);
						}
					}
	
					count++
	                progress = (int)(100f * count / numFeatures)
	            	if (progress != oldProgress) {
						pluginHost.updateProgress("Saving Pruned Streams...", progress)
	            		oldProgress = progress
	            		// check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
	            	}
	            }
	
				outputStreams.write();

				outputStreams = null;
	            
				// display the output streams
				pluginHost.returnData(outputStreamsFile)

				// re-rasterize the network with the pruned value
				linkID = new IntArray2D(rows, cols, -32768)
				linkOutlet = new IntArray2D(rows, cols, -32768)
				linkPosition = new DoubleArray2D(rows, cols, -32768.0)
			
				featureNum = 0;

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
						if (!isBeyondEdge[featureNum - 1] &&
						      linkMag[featureNum - 1] >= thresholdTUCL) {
							box = new BoundingBox();             
							startingPointInPart = partData[part];
		                    if (part < numParts - 1) {
		                        endingPointInPart = partData[part + 1];
		                    } else {
		                        endingPointInPart = numPoints;
		                    }
		
							row = dem.getRowFromYCoordinate(points[startingPointInPart][1])
							col = dem.getColumnFromXCoordinate(points[startingPointInPart][0]);
		                    if (linkID.getValue(row, col) == -32768) {
		                        linkID.setValue(row, col, featureNum);
		                        linkOutlet.setValue(row, col, outletNums[featureNum - 1])
							}
		                    
		                    row = dem.getRowFromYCoordinate(points[endingPointInPart-1][1])
							col = dem.getColumnFromXCoordinate(points[endingPointInPart-1][0]);
		                    if (linkID.getValue(row, col) == -32768) {
		                        linkID.setValue(row, col, featureNum);
		                        linkOutlet.setValue(row, col, outletNums[featureNum - 1])
							}
		                    
		                    double[] distances = new double[endingPointInPart - startingPointInPart + 1];
		                    for (i = startingPointInPart + 1; i < endingPointInPart; i++) {
		                    	distances[i] = distances[i - 1] +  Math.sqrt((points[i][0] - points[i - 1][0]) * (points[i][0] - points[i - 1][0]) + (points[i][1] - points[i - 1][1]) * (points[i][1] - points[i - 1][1]))
		                    }
		                    
		                    for (i = startingPointInPart; i < endingPointInPart; i++) {
		                        if (points[i][0] < box.getMinX()) {
		                            box.setMinX(points[i][0]);
		                        }
		                        if (points[i][0] > box.getMaxX()) {
		                            box.setMaxX(points[i][0]);
		                        }
		                        if (points[i][1] < box.getMinY()) {
		                            box.setMinY(points[i][1]);
		                        }
		                        if (points[i][1] > box.getMaxY()) {
		                            box.setMaxY(points[i][1]);
		                        }
		                    }
		                    topRow = dem.getRowFromYCoordinate(box.getMaxY());
		                    bottomRow = dem.getRowFromYCoordinate(box.getMinY());
		                    leftCol = dem.getColumnFromXCoordinate(box.getMinX());
		                    rightCol = dem.getColumnFromXCoordinate(box.getMaxX());
		
							// find each intersection with a row.
		                    for (row = topRow; row <= bottomRow; row++) {
		
		                        rowYCoord = dem.getYCoordinateFromRow(row);
		                        // find the x-coordinates of each of the line segments 
		                        // that intersect this row's y coordinate
		
		                        for (i = startingPointInPart; i < endingPointInPart - 1; i++) {
		                            if (isBetween(rowYCoord, points[i][1], points[i + 1][1])) {
		                                y1 = points[i][1];
		                                y2 = points[i + 1][1];
		                                if (y2 != y1) {
		                                    x1 = points[i][0];
		                                    x2 = points[i + 1][0];
		
		                                    // calculate the intersection point
		                                    xPrime = x1 + (rowYCoord - y1) / (y2 - y1) * (x2 - x1);
		                                    col = dem.getColumnFromXCoordinate(xPrime);
		                                    
											j = linkID.getValue(row, col);
											if ((j == -32768 || linkMag[j - 1] < linkMag[featureNum - 1]) && dem.getValue(row, col) != nodata) {
			                                    linkID.setValue(row, col, featureNum);
			                                    linkOutlet.setValue(row, col, outletNums[featureNum - 1])
			                                    d = distances[i + 1] - Math.sqrt((xPrime - x2)*(xPrime - x2) + (rowYCoord - y2)*(rowYCoord - y2))
			                                    linkPosition.setValue(row, col, d);
											}
		                                }
		                            }
		                        }
		                    }
		
		                    // find each intersection with a column.
		                    for (col = leftCol; col <= rightCol; col++) {
		                        colXCoord = dem.getXCoordinateFromColumn(col);
		                        for (i = startingPointInPart; i < endingPointInPart - 1; i++) {
		                            if (isBetween(colXCoord, points[i][0], points[i + 1][0])) {
		                                x1 = points[i][0];
		                                x2 = points[i + 1][0];
		                                if (x1 != x2) {
		                                    y1 = points[i][1];
		                                    y2 = points[i + 1][1];
		
		                                    // calculate the intersection point
		                                    yPrime = y1 + (colXCoord - x1) / (x2 - x1) * (y2 - y1);
		                                    row = dem.getRowFromYCoordinate(yPrime);
		
		                                    j = linkID.getValue(row, col);
		                                    if ((j == -32768 || linkMag[j - 1] < linkMag[featureNum - 1]) && dem.getValue(row, col) != nodata) {
			                                    linkID.setValue(row, col, featureNum);
			                                    linkOutlet.setValue(row, col, outletNums[featureNum - 1])
			                                    d = distances[i + 1] - Math.sqrt((yPrime - y2)*(yPrime - y2) + (colXCoord - x2)*(colXCoord - x2))
			                                    linkPosition.setValue(row, col, d);
											}
		                                }
		                            }
		                        }
		                    }
						}
					}
					
					count++
	                progress = (int)(100f * count / numFeatures)
	            	if (progress != oldProgress) {
						pluginHost.updateProgress("Rasterizing Streams...", progress)
	            		oldProgress = progress
	            		// check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
	            	}
	            }

            }


            /*  Scan the raster looking for stream cells
             *  identifying all grid cells that are either the 
             *  maximum or minimum link position value for 
             *  their link ID within their 3 x 3 neighborhood. 
             *  These are link end nodes (either upstream or 
             *  downstream positions). Collisions make it necessary
             *  to search for end-nodes in this way. The 
             *  linkEndNodes array was used above for a different
             *  purpose, i.e. to hold the cells containing feature
             *  end points, which was ued to mark collisions. The
             *  map of end nodes that accounts for the effects of
             *  collisions will be useful for restricing in-stream
             *  flow during the priority-flood op.
             */
            linkEndNodes = new BooleanBitArray2D(rows, cols)
			double position, positionN
			double minNeighbour, maxNeighbour
			oldProgress = -1
			for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					id = linkID.getValue(row, col)
					if (id > 0) { //  it's a stream cell
						position = linkPosition.getValue(row, col)
						minNeighbour = position
						maxNeighbour = position
						for (n = 0; n < 8; n++) {
                    		rowN = row + dY[n];
                    		colN = col + dX[n];
							positionN = linkPosition.getValue(rowN, colN)
							idN = linkID.getValue(rowN, colN)
							if (idN == id) {
								if (positionN < minNeighbour) { minNeighbour = positionN }
								if (positionN > maxNeighbour) { maxNeighbour = positionN }
							}
						}
						if (minNeighbour == position || maxNeighbour == position) {
							linkEndNodes.setValue(row, col, true)
						}
					}
				}
				progress = (int)(100f * row / rowsLessOne)
				if (progress > oldProgress) {
					pluginHost.updateProgress(progress)
					oldProgress = progress

					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}
            
			DoubleArray2D output = new DoubleArray2D(rows, cols, nodata)
			NibbleArray2D flowdir = new NibbleArray2D(rows, cols)
			PriorityQueue<GridCell> queue = new PriorityQueue<GridCell>((2 * rows + 2 * cols) * 2);
			BooleanBitArray2D inQueue = new BooleanBitArray2D(rows, cols)
			
			// initialize the grids
			oldProgress = -1
			for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					z = dem.getValue(row, col)
					output.setValue(row, col, z)
					flowdir.setValue(row, col, 0)
					if (z != nodata) {
						isStream = (linkID.getValue(row, col) > 0)
						isPit = true
						isEdgeCell = false
						for (n = 0; n < 8; n++) {
							zN = dem.getValue(row + dY[n], col + dX[n])
							if (isPit && zN != nodata && zN < z) {
								isPit = false
							} else if (zN == nodata) {
								isEdgeCell = true
							}
						}
						if ((isPit && isEdgeCell) || (isStream && isEdgeCell)) {
							queue.add(new GridCell(row, col, z, 0, isStream))
							inQueue.setValue(row, col, true)
						}
						
					} else {
                        numSolvedCells++
                    }
				}
				progress = (int)(100f * row / rowsLessOne)
				if (progress > oldProgress) {
					pluginHost.updateProgress("Loop 1 of 3", progress)
					oldProgress = progress

					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}

			oldProgress = (int) (100f * numSolvedCells / numCellsTotal);
            pluginHost.updateProgress("Loop 2 of 3: ", oldProgress);
			int flatIndex
			int k
			int linkIDValue, linkIDValueN
			double posDiff
			boolean isAnEndNode, isAnEndNodeN
			
			byte[][] numInflow = new byte[rows][cols]
			int outletNumN;
            while (!queue.isEmpty()) {
                gc = queue.poll();
                row = gc.row;
                col = gc.col;
                outletNum = linkOutlet.getValue(row, col);
                flatIndex = gc.flatIndex;
				linkIDValue = linkID.getValue(row, col);
				z = output.getValue(row, col);
				//isAnEndNode = linkEndNodes.getValue(row, col);
				//if (!(linkIDValue > 0)) { // it's not a stream. only add other non-streams or stream end nodes.
	                for (n = 0; n < 8; n++) {
	                    rowN = row + dY[n];
	                    colN = col + dX[n];
	                    zN = output.getValue(rowN, colN)
	                    if ((zN != nodata) && (!inQueue.getValue(rowN, colN))) {
	                    	linkIDValueN = linkID.getValue(rowN, colN)
	                    	outletNumN = linkOutlet.getValue(rowN, colN)
	                    	isAnEndNodeN = linkEndNodes.getValue(rowN, colN);
	                    	if (linkIDValueN == -32768
	                    	  || (isAnEndNodeN && linkIDValue == -32768) 
	                    	  || (isAnEndNodeN && outletNumN == outletNum)
	                    	  || linkIDValueN == linkIDValue) {
		                	    numSolvedCells++;
		                        flowdir.setValue(rowN, colN, backLink[n])
		                        k = 0
								if (zN == z) {
									k = flatIndex + 1
								}
		                        queue.add(new GridCell(rowN, colN, zN, k, (linkIDValueN > 0)))
								
		                        inQueue.setValue(rowN, colN, true)

		                        numInflow[row][col]++;
	                    	}
	                    }
	                }
//				} else { // it is a stream; add non-streams, link end-nodes, and streams of the same id
//					position = linkPosition.getValue(row, col)
//					int indexOfNextCell = -1
//					double minPosDiff = Integer.MAX_VALUE
//					for (n = 0; n < 8; n++) {
//	                    rowN = row + dY[n];
//	                    colN = col + dX[n];
//	                    zN = output.getValue(rowN, colN)
//	                    if ((zN != nodata) && (!inQueue.getValue(rowN, colN))) {
//	                    	linkIDValueN = linkID.getValue(rowN, colN)
//	                    	if (linkIDValue == linkIDValueN) {
//	                    		positionN = linkPosition.getValue(rowN, colN)
//	                    		posDiff = (positionN - position) * (positionN - position)
//	                    		if (posDiff < minPosDiff) {
//	                    			minPosDiff = posDiff
//	                    			indexOfNextCell = n
//	                    		}
//	                    	} else if (linkIDValueN == -32768 || linkEndNodes.getValue(rowN, colN)) {
//		                	    numSolvedCells++;
//		                        flowdir.setValue(rowN, colN, backLink[n])
//		                        k = 0
//								if (zN == z) {
//									k = flatIndex + 1
//								}
//		                        queue.add(new GridCell(rowN, colN, zN, k, (linkIDValueN > 0)))
//								
//		                        inQueue.setValue(rowN, colN, true)
//	                    	}
//	                    }
//	                }
//
//	                if (indexOfNextCell > -1) {
//	                	n = indexOfNextCell
//	                	rowN = row + dY[n];
//	                    colN = col + dX[n];
//	                	numSolvedCells++;
//                        flowdir.setValue(rowN, colN, backLink[n])
//                        k = 0
//						zN = output.getValue(rowN, colN)
//                        if (zN == z) {
//							k = flatIndex + 1
//						}
//						queue.add(new GridCell(rowN, colN, zN, k, true))
//						
//                        inQueue.setValue(rowN, colN, true)
//	                }
//				}

                progress = (int) (100f * numSolvedCells / numCellsTotal);
                if (progress > oldProgress) {
                    pluginHost.updateProgress("Loop 2 of 3", progress)
                    oldProgress = progress;
                    // check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
                }
            }

            // output the flow direction grid
			WhiteboxRaster outputRaster = new WhiteboxRaster(outputFile, "rw", 
  		  	     demFile, DataType.FLOAT, nodata)
			outputRaster.setPreferredPalette("spectrum.plt")

			oldProgress = -1
			for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					z = output.getValue(row, col)
					if (z != nodata) {
						outputRaster.setValue(row, col, outPointer[flowdir.getValue(row, col)])
					}
				}
				progress = (int)(100f * row / rowsLessOne)
				if (progress > oldProgress) {
					pluginHost.updateProgress("Outputting data...", progress)
					oldProgress = progress

					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}

			outputRaster.addMetadataEntry("Created by the "
	                    + descriptiveName + " tool.")
	        outputRaster.addMetadataEntry("Created on " + new Date())
			Date stop = new Date()
			TimeDuration td = TimeCategory.minus(stop, start)
			outputRaster.addMetadataEntry("Elapsed time: $td")
			outputRaster.close()
	
			// display the output image
			pluginHost.returnData(outputFile)

			// perform the flow accumulation
			WhiteboxRaster outputFA = new WhiteboxRaster(outputFile.replace(".dep", "_FA.dep"), "rw", 
  		  	     demFile, DataType.FLOAT, 1d)
			outputFA.setPreferredPalette("grey.pal");
            outputFA.setDataScale(DataScale.CONTINUOUS);
            outputFA.setZUnits("dimensionless");
            outputFA.setNonlinearity(0.2)
			int rN, cN;
            for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					if (numInflow[row][col] == 0) {
						r = row
						c = col
						numInflow[r][c] = -1
						flag = true
						while (flag) {
							dir = flowdir.getValue(r, c) - 1; //flowDir[r][c]
							if (dir >= 0) {
								rN = r + dY[dir]
								cN = c + dX[dir]
								if (flowdir.getValue(rN, cN) >= 0) {
									z = outputFA.getValue(r, c)
									outputFA.incrementValue(rN, cN, z)
									numInflow[rN][cN] -= 1
									if (numInflow[rN][cN] > 0) {
										flag = false
									} else {
										numInflow[rN][cN] = -1
										r = rN
										c = cN
									}
								} else if (flowdir.getValue(rN, cN) == -1) {
									flag = false
									z = outputFA.getValue(r, c)
									outputFA.incrementValue(rN, cN, z)
								}
							} else {
								flag = false
							}
						}
					}
					if (dem.getValue(row, col) == nodata) {
						outputFA.setValue(row, col, nodata)
					}
				}
				progress = (int)(100f * row / rowsLessOne)
				if (progress != oldProgress) {
					pluginHost.updateProgress("Accumulating flow:", progress)
					oldProgress = progress

					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}

			dem.close()

			outputFA.addMetadataEntry("Created by the "
	                    + descriptiveName + " tool.")
	        outputFA.addMetadataEntry("Created on " + new Date())
			stop = new Date()
			td = TimeCategory.minus(stop, start)
			outputFA.addMetadataEntry("Elapsed time: $td")
			outputFA.close()

			pluginHost.returnData(outputFile.replace(".dep", "_FA.dep"))

				

//            outStream.close()
//            return
            
			// now thin the stream lines
//            int a, i
//            long counter = 1;
//	        int loopNum = 0;
//	        int[][] elements = [ [6, 7, 0, 4, 3, 2], [7, 0, 1, 3, 5], 
//	            [0, 1, 2, 4, 5, 6], [1, 2, 3, 5, 7], 
//	            [2, 3, 4, 6, 7, 0], [3, 4, 5, 7, 1], 
//	            [4, 5, 6, 0, 1, 2], [5, 6, 7, 1, 3] ];
//	        boolean[][] vals = [ [false, false, false, true, true, true], [false, false, false, true, true], 
//	            [false, false, false, true, true, true], [false, false, false, true, true],
//	            [false, false, false, true, true, true], [false, false, false, true, true],
//	            [false, false, false, true, true, true], [false, false, false, true, true] ];
//	        
//	        boolean[] neighbours = new boolean[8];
//	        boolean patternMatch = false;
//
//			while (counter > 0) {
//                loopNum++;
//                pluginHost.updateProgress("Line Thinning (Loop Number " + loopNum + "):", 0);
//                counter = 0;
//                for (row = 0; row < rows; row++) {
//                    for (col = 0; col < cols; col++) {
//                        if (streams.getValue(row, col)) {
//                            // fill the neighbours array
//                            for (i = 0; i < 8; i++) {
//                                neighbours[i] = streams.getValue(row + dY[i], col + dX[i]);
//                            }
//                            
//                            for (a = 0; a < 8; a++) {
//                                // scan through element
//                                patternMatch = true;
//                                for (i = 0; i < elements[a].length; i++) {
//                                    if (neighbours[elements[a][i]] != vals[a][i]) {
//                                        patternMatch = false;
//                                        //break;
//                                    }
//                                }
//                                if (patternMatch) {
//                                    streams.setValue(row, col, false);
//                                    counter++;
//                                }
//                            }
//                        }
//
//                    }
//
//                    progress = (int)(100f * row / rowsLessOne)
//	            	if (progress != oldProgress) {
//						pluginHost.updateProgress("Line Thinning (Loop Number " + loopNum + "):", progress)
//	            		oldProgress = progress
//	            		// check to see if the user has requested a cancellation
//						if (pluginHost.isRequestForOperationCancelSet()) {
//							pluginHost.showFeedback("Operation cancelled")
//							return
//						}
//	            	}
//                }
//            }

//			double[][] output = new double[rows + 2][cols + 2]
//			//double[][] floodOrder = new double[rows + 2][cols + 2]
//			NibbleArray2D flowdir = new NibbleArray2D(rows + 2, cols + 2)
//			PriorityQueue<GridCell> queue = new PriorityQueue<GridCell>((2 * rows + 2 * cols) * 2);
//			BooleanBitArray2D inQueue = new BooleanBitArray2D(rows + 2, cols + 2)
//			
//			// initialize the grids
//			oldProgress = -1
//			for (row = 0; row < rows; row++) {
//				for (col = 0; col < cols; col++) {
//					z = dem.getValue(row, col)
//					output[row + 1][col + 1] = z
//					flowdir.setValue(row + 1, col + 1, 0)
//					if (z != nodata) {
//						isStream = streams.getValue(row, col)
//						isPit = true
//						isEdgeCell = false
//						for (n = 0; n < 8; n++) {
//							zN = dem.getValue(row + dY[n], col + dX[n])
//							if (isPit && zN != nodata && zN < z) {
//								isPit = false
//							} else if (zN == nodata) {
//								isEdgeCell = true
//							}
//						}
//						if ((isPit && isEdgeCell) || (isStream && isEdgeCell)) {
//							queue.add(new GridCell(row + 1, col + 1, z, 0, isStream))
//							inQueue.setValue(row + 1, col + 1, true)
//						}
//						
//					} else {
//                        numSolvedCells++
//                    }
//				}
//				progress = (int)(100f * row / rowsLessOne)
//				if (progress > oldProgress) {
//					pluginHost.updateProgress("Loop 1 of 3", progress)
//					oldProgress = progress
//
//					// check to see if the user has requested a cancellation
//					if (pluginHost.isRequestForOperationCancelSet()) {
//						pluginHost.showFeedback("Operation cancelled")
//						return
//					}
//				}
//			}
//			
//			for (row = 0; row < rows + 2; row++) {
//				output[row][0] = nodata
//				output[row][cols + 1] = nodata
//				flowdir.setValue(row, 0, 0)
//				flowdir.setValue(row, cols + 1, 0)
//			}
//			
//			for (col = 0; col < cols + 2; col++) {
//				output[0][col] = nodata
//				output[rows + 1][col] = nodata
//				flowdir.setValue(0, col, 0)
//				flowdir.setValue(rows + 1, col, 0)
//			}
//
//			oldProgress = (int) (100f * numSolvedCells / numCellsTotal);
//            pluginHost.updateProgress("Loop 2 of 3: ", oldProgress);
//
//			int flatIndex
//			int k
//			//int order = 0
//            while (queue.isEmpty() == false) {
//                gc = queue.poll();
//                row = gc.row;
//                col = gc.col;
//                flatIndex = gc.flatIndex;
//                for (n = 0; n < 8; n++) {
//                    rowN = row + dY[n];
//                    colN = col + dX[n];
//                    zN = output[rowN][colN];
//                    if ((zN != nodata) && (!inQueue.getValue(rowN, colN))) {
//                        numSolvedCells++;
//                        flowdir.setValue(rowN, colN, backLink[n])
//                        k = 0
//						if (zN == output[row][col]) {
//							k = flatIndex + 1
//						}
//                        queue.add(new GridCell(rowN, colN, zN, k, streams.getValue(rowN - 1, colN - 1)))
//						
//                        inQueue.setValue(rowN, colN, true)
//                    }
//                }
//                progress = (int) (100f * numSolvedCells / numCellsTotal);
//                if (progress > oldProgress) {
//                    pluginHost.updateProgress("Loop 2 of 3", progress)
//                    oldProgress = progress;
//                    // check to see if the user has requested a cancellation
//					if (pluginHost.isRequestForOperationCancelSet()) {
//						pluginHost.showFeedback("Operation cancelled")
//						return
//					}
//                }
//            }
//
//            // output the data
//			WhiteboxRaster outputRaster = new WhiteboxRaster(outputFile, "rw", 
//  		  	     demFile, DataType.FLOAT, nodata)
//			outputRaster.setPreferredPalette("spectrum.plt")
//
//			inQueue = new BooleanBitArray2D(rows, cols)
//  		  	oldProgress = -1
//  		  	for (row = 0; row < rows; row++) {
//				for (col = 0; col < cols; col++) {
//					if (streams.getValue(row, col)) {
//						numInflowing = 0			
//		                for (n = 0; n < 8; n++) {
//		                    rowN = row + dY[n];
//		                    colN = col + dX[n];
//		                    if (streams.getValue(rowN, colN) && 
//		                        flowdir.getValue(rowN+1, colN+1) == backLink[n]) { 
//		                        	numInflowing++ 
//		                    }
//		                }
//		                if (numInflowing == 0) {
//		                	r = row
//							c = col
//							zTest = dem.getValue(r, c) - SMALL_NUM
//							flag = true
//							while (flag) {
//								if (zTest > output[r + 1][c + 1] && inQueue.getValue(r, c)) {
//									flag = false // we've already traversed the stream
//								} else {
//									inQueue.setValue(r, c, true)
//									dir = flowdir.getValue(r + 1, c + 1) - 1
//									if (dir >= 0) {
//										// now find the lowest neighbour that isn't the downstream cell
//										lowestNeighbour = LARGE_NUM
//										for (n = 0; n < 8; n++) {
//											if (n != dir) { // we don't expect it to be lower than the cell it flows to
//												zN = output[r + dY[n] + 1][c + dX[n] + 1]
//							                    if (zN != nodata && zN < lowestNeighbour) { 
//							                        lowestNeighbour = zN 
//							                    }
//											}
//						                }
//
//						                zTest = lowestNeighbour - SMALL_NUM
//						                if (zTest < output[r + 1][c + 1]) {
//						                	output[r + 1][c + 1] = zTest
//						                } else {
//											zTest = output[r + 1][c + 1] - SMALL_NUM
//										}
//					                
//										r += dY[dir]
//										c += dX[dir]
//									} else {
//										lowestNeighbour = LARGE_NUM
//										for (n = 0; n < 8; n++) {
//											zN = output[r + dY[n] + 1][c + dX[n] + 1]
//							                if (zN != nodata && zN < lowestNeighbour) { 
//							                    lowestNeighbour = zN 
//							                }
//						                }
//
//						                zTest = lowestNeighbour - SMALL_NUM
//						                if (zTest < output[r + 1][c + 1]) {
//						                	output[r + 1][c + 1] = zTest
//						                }
//
//										flag = false
//									}
//								}
//							}
//							
//		                }
//					}// else {
//					//	outputRaster.setValue(row, col, dem.getValue(row, col))
//					//}
//				}
//  		  		progress = (int)(100f * row / rowsLessOne)
//				if (progress > oldProgress) {
//					pluginHost.updateProgress("Loop 3 of 4", progress)
//					oldProgress = progress
//
//					// check to see if the user has requested a cancellation
//					if (pluginHost.isRequestForOperationCancelSet()) {
//						pluginHost.showFeedback("Operation cancelled")
//						return
//					}
//				}
//			}
//
//			
//
//			oldProgress = -1
//  		  	for (row = 0; row < rows; row++) {
//				for (col = 0; col < cols; col++) {
//					z = output[row + 1][col + 1]
//					outputRaster.setValue(row, col, z)
//					if (streams.getValue(row, col)) {
//						outputRaster.setValue(row, col, 1.0)
//					} else {
//						outputRaster.setValue(row, col, 0.0)
//					}
//				}
//  		  		progress = (int)(100f * row / rowsLessOne)
//				if (progress > oldProgress) {
//					pluginHost.updateProgress("Loop 4 of 4", progress)
//					oldProgress = progress
//
//					// check to see if the user has requested a cancellation
//					if (pluginHost.isRequestForOperationCancelSet()) {
//						pluginHost.showFeedback("Operation cancelled")
//						return
//					}
//				}
//			}

			
			
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
    class GridCell implements Comparable<GridCell> {

        public int row;
        public int col;
        public double z;
        public int flatIndex;
        public boolean streamVal;
        //public long priority

        public GridCell(int row, int col, double z, int flatIndex, boolean streamVal) { //, long priority) { // double Z) {
            this.row = row;
            this.col = col;
            this.z = z;
            this.flatIndex = flatIndex;
            this.streamVal = streamVal;
        }

        @Override
        public int compareTo(GridCell other) {
        	if (this.streamVal && !other.streamVal) {
        		return -1
        	} else if (!this.streamVal && other.streamVal) {
        		return 1
        	} else {
	        	if (this.z > other.z) {
	        		return 1
	        	} else if (this.z < other.z) {
	        		return -1
	        	} else {
	        		if (this.flatIndex > other.flatIndex) {
	        			return 1
	        		} else if (this.flatIndex < other.flatIndex) {
	        			return -1
	        		}
					return 0
	        	}
        	}
        }
    }

    @CompileStatic
    class CollisionCell implements Comparable<CollisionCell> {

        public double tucl;
        
        public CollisionCell(double tucl1, double tucl2) {
        	if (tucl1 < tucl2) {
        		this.tucl = tucl1;
        	} else {
        		this.tucl = tucl2;
        	}
        }

        @Override
        public int compareTo(CollisionCell other) {
        	// sort from highest to lowest
        	if (this.tucl < other.tucl) {
        		return -1;
        	} else if (this.tucl > other.tucl) {
        		return 1;
        	}
	        
	        return 0;
        }
    }

    // Return true if val is between theshold1 and theshold2.
    @CompileStatic
    private static boolean isBetween(double val, double threshold1, double threshold2) {
        if (val == threshold1 || val == threshold2) {
            return true;
        }
        return threshold2 > threshold1 ? val > threshold1 && val < threshold2 : val > threshold2 && val < threshold1;
    }
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def tdf = new TopologicalBreachBurn(pluginHost, args, name, descriptiveName)
}
