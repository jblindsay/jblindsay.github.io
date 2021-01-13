// imports
var Runnable = Java.type('java.lang.Runnable');
var Thread = Java.type('java.lang.Thread');
var WhiteboxRaster = Java.type('whitebox.geospatialfiles.WhiteboxRaster');
var DataType = Java.type('whitebox.geospatialfiles.WhiteboxRasterBase.DataType');

function execute() {
	try {
		var inputDEMFile = "/Users/johnlindsay/Documents/Data/SouthernOnt/DEM_erased.dep"
		var inputStreamsFile = "/Users/johnlindsay/Documents/Data/SouthernOnt/streams.dep"
		var outputFile = "/Users/johnlindsay/Documents/Data/SouthernOnt/tmp1.dep"
		
		var dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ];
		var dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ];
					
		var dem = new WhiteboxRaster(inputDEMFile, "r");
		var rows = dem.getNumberRows();
		var rowsLessOne = rows - 1;
		var columns = dem.getNumberColumns();
		var demNodata = dem.getNoDataValue();
		var maxZ = dem.getMaximumValue();
		
		var streams = new WhiteboxRaster(inputStreamsFile, "r");
		if (streams.getNumberRows() != rows || streams.getNumberColumns() != columns) {
			pluginHost.showFeedback("Input files must have the same dimensions");
			return;
		}
		var streamNodata = dem.getNoDataValue();
		
		var output = new WhiteboxRaster(outputFile, "rw", inputDEMFile, DataType.FLOAT, demNodata);
		output.setPreferredPalette(dem.getPreferredPalette());

		print("Starting...");
		
		var progress, oldProgress = -1;
		var z, zN, streamVal;
		var i;
		for (row = 0; row < rows; row++) {
		    for (col = 0; col < columns; col++) {
		        z = dem.getValue(row, col);
		        if (z !== demNodata) {
		            streamVal = streams.getValue(row, col);
		            if (streamVal > 0 && streamVal != streamNodata) {
		            	output.setValue(row, col, z);
		            } else {
		            	isStreamNeighbour = false;
		            	for (i = 0; i < 8; i++) {
							zN = streams.getValue(row + dY[i], col + dX[i]);
							if (zN > 0 && zN != streamNodata) {
								output.setValue(row, col, maxZ);
								break;
							}
		            	}
										
		            }
		        }
		    }
		    progress = Math.round(row * 100.0 / rowsLessOne);
		    if (progress !== oldProgress) {
		        pluginHost.updateProgress(progress);
		        oldProgress = progress;
		      	print(progress, "%");
		        // check to see if the user has requested a cancellation
		        if (pluginHost.isRequestForOperationCancelSet()) {
		            pluginHost.showFeedback("Operation cancelled");
		            return;
		        }
		    }
		}
		
		dem.close();
		streams.close();
		//output.addMetadataEntry("Created by the " + descriptiveName + " tool.");
		output.addMetadataEntry("Created on " + new Date());
		output.close();

		print("Done!");
		
		// display the output image
		pluginHost.returnData(outputFile);
	} catch (err) {
        pluginHost.showFeedback("An error has occurred:\n" + err);
        pluginHost.logException("Error in " + descriptiveName, err);
    } finally {
        // reset the progress bar
        pluginHost.updateProgress("Progress", 0);
    }
}

var r = new Runnable({
    run: function () {
        execute(args);
    }
});
var t = new Thread(r);
t.start();