var wd = pluginHost.getWorkingDirectory();
	
function filter(lasFileName, minElev) {
	var inputFile = wd + lasFileName;
	var outputFile = wd + filteredName;
	var searchDistance = "5.0";
	var minNeighbours = "5";
	var maxSlope = "35.0";
	var minOGPHgt = "0.0";
	var displayOutput = "false";
	var args = [inputFile, outputFile, searchDistance, minNeighbours, maxSlope, minOGPHgt, minElev, displayOutput];
	pluginHost.runPlugin("IsolateGroundPoints", args, false, true);
}

function interpolate(filteredName, demName) {
	var inputData = wd + filteredName;"";
	var useZValues="true";
	var outputFile = wd + demName;
	var gridRes = "1.5";
	var baseFile = "not specified";
	var weightType = "shepard's";
	var exponent = "2.0";
	var nodalFunc = "constant";
	var maxDistance = "3.0";
	var minNumNeighbours = "5";
	var useQuadSearch = "false";
	var args = [inputData, useZValues, outputFile, 
	gridRes, baseFile, weightType, exponent, nodalFunc, 
	maxDistance, minNumNeighbours, useQuadSearch];
	pluginHost.runPlugin("InterpolationIDW", args, false, false)
}

function mosaic(filesList, fileName) {
	var inputFiles = ""
	for (i = 0; i < filesList.length; i++) {
		if (i < filesList.length - 1) {
			inputFiles += wd + filesList[i] + ";";
		} else {
			inputFiles += wd + filesList[i];
		}
	}
	var outputFile = wd + fileName;
	var resamplingMethod = "nearest neighbour"; 
	var args = [inputFiles, outputFile, resamplingMethod];
	pluginHost.runPlugin("Mosaic", args, false, false);
}

function hillshade(fileName, outFile) {
	var demFile = wd + fileName;
	var outputFile = wd + outFile;
	var azimuth = "315";
	var altitude = "30.0";
	var zFactor = "1.0";
	var args = [demFile, outputFile, azimuth, altitude, zFactor];
	pluginHost.runPlugin("Hillshade", args, false, false);
}

var lasFiles = ["444000_4824000.las", "444000_4825000.las", "444000_4826000.las", "444000_4827000.las",
                "445000_4824000.las", "445000_4825000.las", "445000_4826000.las", "445000_4827000.las",
                "446000_4824000.las", "446000_4825000.las", "446000_4826000.las", "446000_4827000.las",
                "447000_4824000.las", "447000_4825000.las", "447000_4826000.las", "447000_4827000.las",
                "448000_4824000.las", "448000_4825000.las", "448000_4826000.las", "448000_4827000.las"];

var minElevs = ["180.0", "200.0", "200.0", "200.0", "180.0", "200.0", "200.0", "200.0",
                "180.0", "200.0", "200.0", "200.0", "180.0", "200.0", "200.0", "200.0",
                "180.0", "200.0", "200.0", "200.0"];

var filteredFiles = ["tile1 filtered.shp", "tile2 filtered.shp",
                     "tile3 filtered.shp", "tile4 filtered.shp",
                     "tile5 filtered.shp", "tile6 filtered.shp",
                     "tile7 filtered.shp", "tile8 filtered.shp",
                     "tile9 filtered.shp", "tile10 filtered.shp",
                     "tile11 filtered.shp", "tile12 filtered.shp",
                     "tile13 filtered.shp", "tile14 filtered.shp",
                     "tile15 filtered.shp", "tile16 filtered.shp",
                     "tile17 filtered.shp", "tile18 filtered.shp",
                     "tile19 filtered.shp", "tile20 filtered.shp"];
                     
var demFiles = ["tile1.dep", "tile2.dep", "tile3.dep", "tile4.dep",
                "tile5.dep", "tile6.dep", "tile7.dep", "tile8.dep",
                "tile9.dep", "tile10.dep", "tile11.dep", "tile12.dep",
                "tile13.dep", "tile14.dep", "tile15.dep", "tile16.dep",
                "tile17.dep", "tile18.dep", "tile19.dep", "tile20.dep"];

var tileActive = [true, true, true, true,
                  true, true, true, true,
                  true, true, true, true,
                  true, true, true, true,
                  true, true, true, true]

var filterBool = false;
var intepolateBool = true;

for (i = 0; i < lasFiles.length; i++) {
	if (tileActive[i]) {
		if (filterBool) {
			filterAndInterpolate(lasFiles[i], minElevs[i]);
		}
		if (intepolateBool) {
			interpolate(filteredFiles[i], demFiles[i]);
		}
		if (pluginHost.isRequestForOperationCancelSet()) {
			break;
		}
	}
}

mosaic(demFiles, "DEM IDW.dep");
hillshade("DEM IDW.dep", "Hillshade IDW.dep");

pluginHost.showFeedback("Complete!");