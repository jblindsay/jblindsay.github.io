/* Created by John Lindsay, 2015 
 * This script will take three input images, scale 
 * their absolute values to a 'cutoff' value and then
 * combine them into a single RGB image.
 */
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.WhiteboxRasterBase.DataScale

try {
double cutoff = 2.58

//def redImageName = "/Users/johnlindsay/Documents/Research/Multi-scale Topographic Position paper/data/Appalachians/UTM/large mag3.dep"
//def greenImageName = "/Users/johnlindsay/Documents/Research/Multi-scale Topographic Position paper/data/Appalachians/UTM/med mag3.dep"
//def blueImageName = "/Users/johnlindsay/Documents/Research/Multi-scale Topographic Position paper/data/Appalachians/UTM/small mag3.dep"
//def hillshadeImageName = "/Users/johnlindsay/Documents/Research/Multi-scale Topographic Position paper/data/Appalachians/UTM/hillshade.dep"
//def outputImageName = "/Users/johnlindsay/Documents/Research/Multi-scale Topographic Position paper/data/Appalachians/UTM/colour_comp3.dep"

//def redImageName = "/Users/johnlindsay/Documents/Data/West_Coast_SRTM/broad mag.dep"
//def greenImageName = "/Users/johnlindsay/Documents/Data/West_Coast_SRTM/meso mag.dep"
//def blueImageName = "/Users/johnlindsay/Documents/Data/West_Coast_SRTM/local mag.dep"
//def outputImageName = "/Users/johnlindsay/Documents/Data/West_Coast_SRTM/colour comp.dep"

//def redImageName = "/Users/johnlindsay/Documents/Data/SouthernOnt/DEVmax broad.dep"
//def greenImageName = "/Users/johnlindsay/Documents/Data/SouthernOnt/DEVmax meso.dep"
//def blueImageName = "/Users/johnlindsay/Documents/Data/SouthernOnt/DEVmax local.dep"
//def outputImageName = "/Users/johnlindsay/Documents/Data/SouthernOnt/colour comp.dep"

def redImageName = "/Users/johnlindsay/Documents/Research/Multi-scale Topographic Position paper/data/Hawaii/tmp5.dep"
def greenImageName = "/Users/johnlindsay/Documents/Research/Multi-scale Topographic Position paper/data/Hawaii/tmp7.dep"
def blueImageName = "/Users/johnlindsay/Documents/Research/Multi-scale Topographic Position paper/data/Hawaii/tmp6.dep"
def outputImageName = "/Users/johnlindsay/Documents/Research/Multi-scale Topographic Position paper/data/Hawaii/MTPCC1.dep"

def redraster = new WhiteboxRaster(redImageName, "r")
def rows = redraster.getNumberRows()
def cols = redraster.getNumberColumns()
def nodata = redraster.getNoDataValue()

def greenraster = new WhiteboxRaster(greenImageName, "r")
def blueraster = new WhiteboxRaster(blueImageName, "r")
//def hillshaderaster = new WhiteboxRaster(hillshadeImageName, "r")

def outputraster = new WhiteboxRaster(outputImageName, "rw", redImageName, DataType.FLOAT, nodata)
outputraster.setPreferredPalette("rgb.pal");
outputraster.setDataScale(DataScale.RGB);

double z, hs
int progress, r, g, b
def oldprogress = -1
for (int row = 0; row < rows; row++) {
	for (int col = 0; col < cols; col++) {
		redVal = redraster.getValue(row, col)
        greenVal = greenraster.getValue(row, col)
        blueVal = blueraster.getValue(row, col)
        if ((redVal != nodata) && (greenVal != nodata) && (blueVal != nodata)) {
        	redVal = Math.abs(redVal)
        	greenVal = Math.abs(greenVal)
        	blueVal = Math.abs(blueVal)
        	if (redVal > cutoff) { redVal = cutoff }
        	if (greenVal > cutoff) { greenVal = cutoff }
        	if (blueVal > cutoff) { blueVal = cutoff }

        	//hs = hillshaderaster.getValue(row, col) / 150.0
        	//if (hs > 1) { hs = 1; }
        	
            r = (int)(redVal / cutoff * 255) //* hs;
            if (r < 0) {
                r = 0;
            }
            if (r > 255) {
                r = 255;
            }
            g = (int)(greenVal / cutoff * 255) //* hs;
            if (g < 0) {
                g = 0;
            }
            if (g > 255) {
                g = 255;
            }
            b = (int)(blueVal / cutoff * 255) //* hs;
            if (b < 0) {
                b = 0;
            }
            if (b > 255) {
                b = 255;
            }
            z = (double) ((255 << 24) | (b << 16) | (g << 8) | r);
            outputraster.setValue(row, col, z);
        }
	}
	progress = (int)(100.0 * row / (rows - 1))
	if (progress != oldprogress) {
		oldprogress = progress
		pluginHost.updateProgress(progress)
		if (pluginHost.isRequestForOperationCancelSet()) {
			pluginHost.showFeedback("Operation cancelled")
			return
		}
	}
}
redraster.close()
greenraster.close()
blueraster.close()
//hillshaderaster.close()
outputraster.addMetadataEntry("Created by a custom script.");
outputraster.addMetadataEntry("Created on " + new Date());
outputraster.close()

pluginHost.returnData(outputImageName)

pluginHost.updateProgress(-1)
} catch (Exception e) {
	println("Error!")
	println(e)
}
