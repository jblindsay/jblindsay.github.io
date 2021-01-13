var WhiteboxRaster = Java.type('whitebox.geospatialfiles.WhiteboxRaster');
var DataType = Java.type('whitebox.geospatialfiles.WhiteboxRasterBase.DataType');

var descriptiveName = "Fill DEM voids with other elevation data";

var file1 = "/Users/johnlindsay/Documents/Data/AlbertaLidar/tmp4.dep";
var file2 = "/Users/johnlindsay/Documents/Data/AlbertaLidar/Alberta/OldLidar.dep";
var file3 = "/Users/johnlindsay/Documents/Data/AlbertaLidar/DEM_final_no_gaps.dep";

var input1 = new WhiteboxRaster(file1, "r");
var rows = input1.getNumberRows();
var cols = input1.getNumberColumns();
var nodata = input1.getNoDataValue();

var input2 = new WhiteboxRaster(file2, "r");
var nodata2 = input2.getNoDataValue();

var output = new WhiteboxRaster(file3, "rw", file1, DataType.FLOAT, nodata);
output.setNoDataValue(nodata);
output.setPreferredPalette("high_relief.pal");

var z1, z2, x, y
var r, c
var progress, oldProgress = -1;
for (row = 0; row < rows; row++) {
	y = input1.getYCoordinateFromRow(row);
	for (col = 0; col < cols; col++) {
		var z1 = input1.getValue(row, col);
		if (z1 === nodata) {
			x = input1.getXCoordinateFromColumn(col);
			r = input2.getRowFromYCoordinate(y);
			c = input2.getColumnFromXCoordinate(x);
			var z2 = input2.getValue(r, c);
			if (z2 !== nodata2) {
				output.setValue(row, col, z2);
			} else {
				output.setValue(row, col, nodata);
			}
		} else {
			output.setValue(row, col, z1);
		}
	}
	progress = row * 100.0 / (rows - 1);
	if (progress !== oldProgress) {
		pluginHost.updateProgress(progress);
		oldProgress = progress;
		// check to see if the user has requested a cancellation
//		if (pluginHost.isRequestForOperationCancelSet()) {
//            pluginHost.showFeedback("Operation cancelled");
//			break;
//        }
	}
}

// reset the progress bar
pluginHost.updateProgress(-1);

input1.close();
input2.close();
output.addMetadataEntry("Created by the " + descriptiveName + " tool.");
output.addMetadataEntry("Created on " + new Date());
output.close();

// display the output image
pluginHost.returnData(file3);
