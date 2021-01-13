import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.attributes.DBFField
import whitebox.geospatialfiles.shapefile.attributes.*

def inputFile = pluginHost.getWorkingDirectory() + "ground points.shp"
def outputFile = pluginHost.getWorkingDirectory() + "ground points only.shp"

ShapeFile shape = new ShapeFile(inputFile)

AttributeTable table = shape.getAttributeTable()
DBFField[] fields = table.getAllFields()

ShapeFile output = new ShapeFile(outputFile, shape.getShapeType(), fields)

int numRecords = shape.getNumberOfRecords()
String fieldName = "CLASS"
int i = 0
Object[] rowData
ShapeFileRecord record
Geometry geom
int progress
int oldProgress = -1
for (int rec in 0..<numRecords) {
	double val = (double)table.getValue(rec, fieldName)
	if (val == 1.0d) { 
		record = shape.getRecord(rec)
		geom = record.getGeometry()
		rowData = table.getRecord(rec)
		output.addRecord(geom, rowData)
		i++
	}

	progress = (int)(100f * rec / (numRecords - 1))
	if (progress != oldProgress) {
		pluginHost.updateProgress(progress)
		oldProgress = progress
	}
}

output.write()

println("${i} records found")
pluginHost.returnData(outputFile)
pluginHost.updateProgress(0)
