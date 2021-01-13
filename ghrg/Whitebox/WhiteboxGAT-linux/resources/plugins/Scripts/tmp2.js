// imports
var Runnable = Java.type('java.lang.Runnable');
var Thread = Java.type('java.lang.Thread');
var WhiteboxRaster = Java.type('whitebox.geospatialfiles.WhiteboxRaster');
var DataType = Java.type('whitebox.geospatialfiles.WhiteboxRasterBase.DataType');

function execute() {
	try {
		var inputDEMFile1 = "/Users/johnlindsay/Documents/Data/AlbertaLidar/tmp4.dep"
		var inputDEMFile2 = "/Users/johnlindsay/Documents/Data/AlbertaLidar/DEM_final_no_gaps.dep"
		var outputFile = "/Users/johnlindsay/Documents/Data/AlbertaLidar/DEM_10m_no_gaps.dep"
		
		var dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ];
		var dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ];
					
		var dem = new WhiteboxRaster(inputDEMFile1, "r");
		var rows = dem.getNumberRows();
		var rowsLessOne = rows - 1;
		var columns = dem.getNumberColumns();
		var demNodata = dem.getNoDataValue();
		var maxZ = dem.getMaximumValue();
		
		var dem2 = new WhiteboxRaster(inputDEMFile2, "r");
		var dem2Nodata = dem2.getNoDataValue();
		
		var output = new WhiteboxRaster(outputFile, "rw", inputDEMFile1, DataType.FLOAT, demNodata);
		output.setPreferredPalette(dem.getPreferredPalette());

		print("Starting...");
		
		var progress, oldProgress = -1;
		var z, zN, z2, x, y, row2, col2;
		var i;
		for (row = 0; row < rows; row++) {
		    for (col = 0; col < columns; col++) {
		        z = dem.getValue(row, col);
		        if (z === demNodata) {
		        	x = dem.getXCoordinateFromColumn(col);
		        	y = dem.getYCoordinateFromRow(row);
		        	col2 = dem2.getColumnFromXCoordinate(x);
		        	row2 = dem2.getRowFromYCoordinate(y);
		            z2 = dem2.getValue(row2, col2);
		            if (z2 !== dem2Nodata) {
		            	output.setValue(row, col, z2);
		            } //else {
//		            	isStreamNeighbour = false;
//		            	for (i = 0; i < 8; i++) {
//							zN = streams.getValue(row + dY[i], col + dX[i]);
//							if (zN > 0 && zN != streamNodata) {
//								output.setValue(row, col, maxZ);
//								break;
//							}
//		            	}
//										
//		            }
		        } else {
		        	output.setValue(row, col, z);
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
		dem2.close();
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
