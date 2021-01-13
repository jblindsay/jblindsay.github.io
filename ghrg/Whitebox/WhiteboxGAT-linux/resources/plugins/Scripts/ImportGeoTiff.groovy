import java.awt.image.BufferedImage;
import java.io.File;
import java.util.List;
import java.nio.ByteBuffer;
import java.awt.image.ColorModel;
import org.apache.commons.imaging.ImageInfo;
import org.apache.commons.imaging.Imaging;
import org.apache.commons.imaging.common.IImageMetadata;
import org.apache.commons.imaging.formats.tiff.taginfos.*
import org.apache.commons.imaging.formats.tiff.fieldtypes.*
import org.apache.commons.imaging.formats.tiff.constants.TiffDirectoryType

def fileName = "/Users/johnlindsay/Documents/Data/ned10m38122b5.tif/ned10m38122b5.tif"
//def fileName = "/Users/johnlindsay/Documents/Data/GeoTiff/Sample64BitFloatingPointPix451x337.tif"

def file = new File(fileName)
final TiffImageMetadata metadata = Imaging.getMetadata(file);
//println metadata
metadata.getAllFields().each() { it -> 
	println "Name: ${it.getTagName()}"
	println "Tag: ${it.getTag()}"
	println "FieldType Name: ${it.getFieldType().getName()}"
	println "FieldType Type: ${it.getFieldType().getType()}"
	println "FieldType Size: ${it.getFieldType().getSize()}"
	println "Value: ${it.getValue()}"
	println ""
}


TagInfoShortOrLong colInfo = new TagInfoShortOrLong("ImageWidth", 256, 1, TiffDirectoryType.TIFF_DIRECTORY_ROOT)
int cols = metadata.findField(colInfo).getValue();
println("Columns: $cols")

TagInfoShortOrLong rowInfo = new TagInfoShortOrLong("ImageLength", 257, 1, TiffDirectoryType.TIFF_DIRECTORY_ROOT)
int rows = metadata.findField(rowInfo).getValue();
println("Rows: $rows")

TagInfoShort sfInfo = new TagInfoShort("SampleFormat", 339, 1, TiffDirectoryType.TIFF_DIRECTORY_ROOT)
int sampleFormat = metadata.findField(sfInfo).getValue();
println("Sample Format: $sampleFormat")

TagInfoShort piInfo = new TagInfoShort("PhotometricInterpretation", 262, 1, TiffDirectoryType.TIFF_DIRECTORY_ROOT)
int photometricInterpretation = metadata.findField(piInfo).getValue();
println("Photometric Interpretation: $photometricInterpretation")


TagInfoDouble tiPixScale = new TagInfoDouble("ModelPixelScaleTag", 33550, 8, TiffDirectoryType.EXIF_DIRECTORY_UNKNOWN)
TagInfoDouble tiTiePoint = new TagInfoDouble("ModelTiepointTag", 33922, 8, TiffDirectoryType.EXIF_DIRECTORY_UNKNOWN)

double[] modelPixelScale = metadata.findField(tiPixScale).getValue();
double[] modelTiepoint = metadata.findField(tiTiePoint).getValue()

double north = modelTiepoint[4] + modelTiepoint[1] * modelPixelScale[1];
double south = modelTiepoint[4] - (rows - modelTiepoint[1]) * modelPixelScale[1];
double east = modelTiepoint[3] + (cols - modelTiepoint[0]) * modelPixelScale[0];
double west = modelTiepoint[3] - modelTiepoint[0] * modelPixelScale[0];
            
println("North: $north")
println("South: $south")
println("East: $east")
println("West: $west")

TagInfoAscii nodataInfo = new TagInfoAscii("GDALNoData", 42113, 1, TiffDirectoryType.EXIF_DIRECTORY_UNKNOWN)
double nodata = Double.parseDouble(metadata.findField(nodataInfo).getValue());
println("NoData: $nodata")



final ImageInfo imageInfo = Imaging.getImageInfo(file);
//println imageInfo
final BufferedImage image = Imaging.getBufferedImage(file);

int r, g, b, a;
Integer z
double val
z = image.getRGB(501, 400);
val = (double)z.floatValue()
println(val)
b = (int) z & 0xFF;
g = ((int) z >> 8) & 0xFF;
r = ((int) z >> 16) & 0xFF;
a = ((int) z >> 24) & 0xFF;

println("$z $r $g $b $a")
//for (int row = 0; row < rows; row++) {
//    for (int col = 0; col < cols; col++) {
//        z = image.getRGB(col, row);
//        r = (int) z & 0xFF;
//        g = ((int) z >> 8) & 0xFF;
//        b = ((int) z >> 16) & 0xFF;
//        output.setValue(row, col, (double) ((255 << 24) | (b << 16) | (g << 8) | r));
//    }
//}
                        
//ColorModel cm = image.getColorModel();
//println cm.getPixelSize()
//if (image.getColorModel().getPixelSize() == 24) {
//    dataType = "float";
//    dataScale = "rgb";
//}
