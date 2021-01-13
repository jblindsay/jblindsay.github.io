/*
 * This script will use the Run Plugin In Parallel tool to concurrently 
 * process a workflow in which dozens of shapefiles containing LiDAR 
 * points are processed for identifying their ground points.
 */
String pluginName = "IsolateGroundPoints"
String fieldName = "Z"
String suffix = "grd"
String searchDistance = "0.5"
String minPoints = "10"
String slopeThreshold = "25"
String displayOutput = "false"

def tag = "25deg_0_5radius"

/* The lidar shapefiles each have a name like 'CVC_all_RowY_ColX.shp'
 *  and are all in the same directory. The code below is used to 
 *  retreive each of these shapefiles and add a line in the string 
 *  that is fed as an argument to the RunPluginInParallel script.
 */

StringBuilder sb = new StringBuilder()

int numFiles = 0
for (int row = 1; row < 50; row++) {
	for (int col = 1; col < 50; col++) {
		String fileName = pluginHost.getWorkingDirectory() + "CVC_all_Row${row}_Col${col}.shp"
		File file = new File(fileName)
		if (file.exists()) {
			String line = "${fileName},${fieldName},${suffix},${searchDistance},${minPoints},${slopeThreshold},${displayOutput}\n"
			sb.append(line)
			numFiles++
		}
	}
}

def suppressReturns = "true" 
String[] args = [pluginName, sb.toString(), suppressReturns] 
pluginHost.runPlugin("RunPluginInParallel", args, false, true)

// merge the points
sb = new StringBuilder()
for (int row = 1; row < 40; row++) {
	for (int col = 1; col < 40; col++) {
		String fileName = pluginHost.getWorkingDirectory() + "CVC_all_Row${row}_Col${col} ${suffix}.shp"
		File file = new File(fileName)
		if (file.exists()) {
			sb.append(fileName).append(";")
		}
	}
}
def mergedPoints = pluginHost.getWorkingDirectory() + "tmp1.shp"
def inputFiles = (sb.toString()).substring(0, sb.toString().length() - 1)
args = [inputFiles, mergedPoints]
pluginHost.runPlugin("MergePointsFiles", args, false, true) 

// clip the points
def clipFile = pluginHost.getWorkingDirectory() + "CVC mask.shp" 
def clippedPointsFile = pluginHost.getWorkingDirectory() + "ground points clipped ${tag}.shp"
args = [mergedPoints, clipFile, clippedPointsFile]
pluginHost.runPlugin("Clip", args, false, true)

// interpolate the grid
def inputData = clippedPointsFile + ";Z" 
def useZ = "false"
def demFile = pluginHost.getWorkingDirectory() + "ground surface ${tag}.dep" 
def gridRes = "0.05" 
def baseFileHeader = ""
def weightType = "shepard"
def weight = "2.0" 
def nodalFunc = "constant"
def searchDist = "0.15"
def numNeighbours = "5"
args = [inputData, useZ, demFile, gridRes, baseFileHeader, weightType, weight, nodalFunc, searchDist, numNeighbours, "false"] 
pluginHost.runPlugin("InterpolationIDW", args, false, true) 

// clip the DEM to the mask
def maintainDimensions = "false" 
args = [demFile, clipFile, maintainDimensions] 
pluginHost.runPlugin("ClipRasterToPolygon", args, false)

 
pluginHost.showFeedback("Done")
