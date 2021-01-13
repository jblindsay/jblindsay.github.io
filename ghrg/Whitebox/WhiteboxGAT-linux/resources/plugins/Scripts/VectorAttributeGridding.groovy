import java.awt.event.ActionListener
import java.awt.event.ActionEvent
import java.io.File
import java.util.Date
import java.util.ArrayList
import java.util.HashSet
import java.util.concurrent.Future
import java.util.concurrent.*
import java.util.concurrent.atomic.AtomicInteger
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.ui.plugin_dialog.*
import whitebox.utilities.FileUtilities;
import whitebox.geospatialfiles.VectorLayerInfo
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.geospatialfiles.shapefile.attributes.DBFField.DBFDataType
import whitebox.geospatialfiles.shapefile.ShapeFileRecord
import whitebox.utilities.Topology
import whitebox.structures.KdTree
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "VectorAttributeGridding"
def descriptiveName = "Vector Attribute Gridding"
def description = "Maps the spatial pattern of an attribute associated with vector features."
def toolboxes = ["VectorTools", "RasterCreation"]

public class VectorAttributeGridding implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd
	private String descriptiveName
	private AtomicInteger numSolved = new AtomicInteger(0)
	private int rows
	//private int numSolvedTiles = 0
	
	public VectorAttributeGridding(WhiteboxPluginHost pluginHost, String[] args, String name, String descriptiveName) {

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
			//sd.addDialogFile("Input vector file", "Input Vector File:", "open", "Vector Files (*.shp), SHP", true, false)
            DialogFieldSelector dfs = sd.addDialogFieldSelector("Input file and Value field.", "Input Value Field:", false)															
			DialogDataInput di = sd.addDialogDataInput("Output raster cell size.", "Cell Size (optional):", "", true, true)														
			DialogFile df = sd.addDialogFile("Input base file", "Base Raster File (optional):", "open", "Whitebox Raster Files (*.dep), DEP", true, true)
			sd.addDialogDataInput("Search distance.", "Search Distance:", "", true, false)
            sd.addDialogFile("Output raster file", "Output Raster File:", "saveAs", "Whitebox Raster Files (*.dep), DEP", true, false)

			def lstr2 = { evt -> if (evt.getPropertyName().equals("value")) { 
            		String value = df.getValue()
            		if (!value.isEmpty() && value != null) { 
            			di.setValue("")
            		}
            	} 
            } as PropertyChangeListener
            df.addPropertyChangeListener(lstr2)

            
			// resize the dialog to the standard size and display it
			sd.setSize(800, 400)
			sd.visible = true
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

	// The CompileStatic annotation can be used to significantly
	// improve the performance of a Groovy script to nearly 
	// that of native Java code.
	@CompileStatic
	private void execute(String[] args) {
		try {

			double nodata = -32768.0
			int cols
			
			if (args.length != 5) {
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
			String inputShapefile = inputData[0]
			String fieldName = inputData[1]

			double cellSize = -1.0
			if (args[1] != null && !args[1].isEmpty() && !args[1].toLowerCase().equals("not specified")) {
				cellSize = Double.parseDouble(args[1])
			}
			String inputBaseRaster = args[2]
			if (args[2] != null && !args[2].isEmpty() && !args[2].toLowerCase().equals("not specified")) {
				inputBaseRaster = args[2]
			}
			
			double searchDist = Double.parseDouble(args[3])
			String outputRasterFile = args[4]

			ShapeFile input = new ShapeFile(inputShapefile)

			// see if the input attributes is numerical
            DBFField[] fields = input.getAttributeTable().getAllFields()
        	int fieldNum = input.getAttributeTable().getFieldColumnNumberFromName(fieldName)
			if (fields[fieldNum].getDataType() != DBFDataType.NUMERIC && 
        	     fields[fieldNum].getDataType() != DBFDataType.FLOAT) {
        	    pluginHost.showFeedback("The input attribute must be of numeric type. If you have \n" +
        	                            "categorical attributes, use the Field Calculator to create a new \n" + 
        	                            "Dummy variable from this attribute (i.e. convert it to \n" + 
        	                            "numerical data).")
        	    return
        	}
			
			WhiteboxRaster output
			if ((cellSize > 0) || ((cellSize < 0) & (inputBaseRaster.toLowerCase().contains("not specified")))) {
                if ((cellSize < 0) & (inputBaseRaster.toLowerCase().contains("not specified"))) {
                    cellSize = Math.min((input.getyMax() - input.getyMin()) / 500.0,
                            (input.getxMax() - input.getxMin()) / 500.0);
                }
                double north = input.getyMax() + cellSize / 2.0;
                double south = input.getyMin() - cellSize / 2.0;
                double east = input.getxMax() + cellSize / 2.0;
                double west = input.getxMin() - cellSize / 2.0;
                rows = (int) (Math.ceil((north - south) / cellSize));
                cols = (int) (Math.ceil((east - west) / cellSize));

                // update west and south
                east = west + cols * cellSize;
                south = north - rows * cellSize;

                output = new WhiteboxRaster(outputRasterFile, north, south, east, west,
                        rows, cols, WhiteboxRasterBase.DataScale.CONTINUOUS,
                        WhiteboxRasterBase.DataType.FLOAT, nodata, nodata);
                
            } else {
                output = new WhiteboxRaster(outputRasterFile, "rw",
                        inputBaseRaster, WhiteboxRasterBase.DataType.FLOAT, nodata);
                rows = output.getNumberRows()
                cols = output.getNumberColumns()
            }

            pluginHost.updateProgress("Please wait...", 0)

            int numProcessors = Runtime.getRuntime().availableProcessors()
            int multiplier = 8
            int numBlocks = (int)(numProcessors * multiplier)
            int blockSize = (int)(rows / numBlocks)
            if (blockSize < 1) { blockSize = 1 }
            //numBlocks = rows / blockSize

            ArrayList<DoWorkOnRow> tasks = new ArrayList<>()
            int stRow = 0
            int endingRow = -1
            while (stRow < (rows - 1) || endingRow < (rows - 1)) {
            	endingRow = stRow + blockSize
            	if (endingRow >= rows) { endingRow = rows - 1 }
            	int nRows = endingRow - stRow + 1
				if (stRow < rows) {
					tasks.add(new DoWorkOnRow(input, fieldName, stRow, nRows, 
		    		cols, searchDist, nodata, output.getNorth(), 
		    		output.getWest(), output.getCellSizeY(),
		    		output.getCellSizeX()))
				}
				stRow += blockSize + 1
            }
	    
			ExecutorService executor = Executors.newFixedThreadPool(numProcessors)
	  	    // the only reason for the getExecutorResults method 
	  	    // is that Groovy throws a compilation type mis-match
	  	    // error when compiled statically. I think it's a bug.
	  	    List<Future<RowNumberAndData>> results = getExecutorResults(executor, tasks)
	        executor.shutdown()

			int progress
	        int i = 0
	        int oldProgress = -1
		    for (Future<RowNumberAndData> result : results) {
	    		RowNumberAndData data = result.get()
	    		int startingRow = data.getRow()
	    		double[][] rowData = data.getData()
	    		for (int r in 0..<rowData.length) {
	    			//output.setRowValues(startingRow + r, rowData[r])
	    			for (int c in 0..<rowData[r].length) {
	    				output.setValue(startingRow + r, c, rowData[r][c])
	    			}
	    			i++
	    		}
	        	
	        	// update progress bar
				progress = (int)(100f * i / rows)
				if (progress > oldProgress) {
					pluginHost.updateProgress("Progress", progress)
					oldProgress = progress
				}
				// check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
		    }
		
			output.setPreferredPalette("spectrum.plt")
			output.addMetadataEntry("Created by the "
                    + descriptiveName + " tool.")
	        output.addMetadataEntry("Created on " + new Date())
			output.close()
			
			pluginHost.returnData(outputRasterFile)
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

	public List<Future<RowNumberAndData>> getExecutorResults(ExecutorService executor, ArrayList<DoWorkOnRow> tasks) {
    	List<Future<RowNumberAndData>> results = executor.invokeAll(tasks)
		return results
    }

    boolean cancelledMessageGiven = false

    @CompileStatic
    class DoWorkOnRow implements Callable<RowNumberAndData> {
		private int startingRow
		private int nRows
		private int cols
		private ShapeFile input
		private double nodata
		private double searchDist
		private double north
		private double west
		private double cellSizeX
		private double cellSizeY
		private String fieldName
		
	    DoWorkOnRow(ShapeFile input, String fieldName, int startingRow, int nRows, 
	    	int cols, double searchDist, double nodata,
	       double north, double west, double cellSizeY, 
	       double cellSizeX) {
        	this.input = input
        	this.startingRow = startingRow
        	this.nRows = nRows
        	this.cols = cols
        	this.searchDist = searchDist
            this.nodata = nodata
            this.north = north
            this.west = west
            this.cellSizeX = cellSizeX
            this.cellSizeY = cellSizeY
            this.fieldName = fieldName
       	}
        	
        @Override
	    public RowNumberAndData call() {
	    	int solved = 0
	    	// check to see if the user has requested a cancellation
			if (pluginHost.isRequestForOperationCancelSet()) {
				if (!cancelledMessageGiven) { 
					pluginHost.showFeedback("Operation cancelled")
					cancelledMessageGiven = true
				}
				return
			}

			double halfCellSizeX = cellSizeX / 2.0d
			double halfCellSizeY = cellSizeY / 2.0d
			double area = Math.PI * searchDist * searchDist / 1000000
			
	    	double[][] retData = new double[nRows][cols]

	       	KdTree tree = input.getKdTree()

	       	int numFeatures = input.getNumberOfRecords()
	       	AttributeTable table = input.getAttributeTable()
	       	double[] fieldData = new double[numFeatures]
	       	for (int rec = 0; rec < numFeatures; rec++) {
	       		double val = (double)table.getValue(rec, fieldName)
	       		fieldData[rec] = val
	       	}
			
			List<KdTree.Entry<Integer>> results
			double x, y
			for (int row in 0..<nRows) {
				y = north - halfCellSizeY - (startingRow + row) * cellSizeY
				for (int col in 0..<cols) {
					x = west + halfCellSizeX + col * cellSizeX;
			        double[] entry = new double[2]
			        entry[0] = x
			        entry[1] = y
			        
			        results = tree.neighborsWithinRange(entry, searchDist)
			        if (results.size() > 0) {
				        HashSet<Integer> hs = new HashSet() 
				        for (int j = 0; j < results.size(); j++) {
				        	KdTree.Entry result = results.get(j)
							if (result.distance > 0) {
								int k = result.value
								hs.add(k)
							}
				        }
				        double total = 0
				        hs.each() { it -> 
							int val = (int)it
							total += fieldData[val]
						}
			        	retData[row][col] = (double)(total / hs.size())
			        } else {
			        	retData[row][col] = nodata
			        }
				}
				if (pluginHost.isRequestForOperationCancelSet()) {
					if (!cancelledMessageGiven) { 
						pluginHost.showFeedback("Operation cancelled")
						cancelledMessageGiven = true
					}
					return
				}

				solved = numSolved.incrementAndGet()
			}
  			
	        def ret = new RowNumberAndData(startingRow, retData)
	        
	        int progress = (int) (100f * solved / rows)
			pluginHost.updateProgress("Interpolated $solved rows:", progress)
	        
    	    return ret
        }               
    }

	@CompileStatic
    class RowNumberAndData {
		private int startingRow
		private int nRows
		private double[][] data
		private int nCols = 0
		
    	RowNumberAndData(int startingRow, double[][] data) {
    		this.startingRow = startingRow
    		this.data = data
    		this.nRows = data.length
    		this.nCols = data[0].length
    	}

    	public int getRow() {
    		return startingRow
    	}

    	public void setRow(int startingRow) {
    		this.startingRow = startingRow
    	}

    	public int getNumRows() {
    		return nRows
    	}

    	public void setNumRows(int nRows) {
    		this.nRows = nRows
    	}

    	public double[][] getData() {
    		return data
    	}

		public void setData(double[][] data) {
			this.data = data
			this.nRows = data.length
    		this.nCols = data[0].length
		}
		
    	public double getDataAt(int row, int col) {
    		if (row >= 0 && row < nRows && col >= 0 && col < nCols) {
    			return data[row][col]
    		} else {
    			return null
    		}
    	}
    }
	
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def myClass = new VectorAttributeGridding(pluginHost, args, name, descriptiveName)
}
