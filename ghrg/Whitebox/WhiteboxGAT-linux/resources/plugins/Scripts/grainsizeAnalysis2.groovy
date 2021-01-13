import whitebox.utilities.FileUtilities

/* Jackie, the script will output the data as a tab-delimited table
 *  in the output window below. You should be able to copy it 
 *  (command a to select all then command c to copy) and paste it 
 *  directly into a blank Excel sheet.
 */
 
/* Jackie, change the following 'directory variable to 
 *  a directory path on your system that contains all 
 *  of the various files you want to process. They must 
 *  be in a csv format. If you have too many data files 
 *  to be able to export them all from Excel in csv 
 *  format in a reasonable amount of time, let me know
 *  and I'll see about writing a script that will do 
 *  it automatically. This is just a bit more involved.
 */
//def directory = "/Users/jcockburn/data/grainsize files/"
def directory = "/Users/johnlindsay/Dropbox/BF2 CSV" // my test directory...it contained only the one file you sent me.

// change these next two as needed
def fileExtension = "csv"
def delimiter = "\t"

def dataFiles = FileUtilities.findAllFilesWithExtension(directory, fileExtension, true)
int numFiles = dataFiles.size()

// create an array to store the file names
String[] fileNames = new String[numFiles]

// first figure out how many grain size classes there are
int numClasses = 0
def startClassRead = false
new File(dataFiles.get(0)).eachLine { line ->
	if (startClassRead) {
		def cellData = line.split(delimiter)
		if (!cellData[0].trim().isEmpty()) {
			numClasses++
		}
	}
	if (line.contains("um${delimiter}Volume")) {
		startClassRead = true
	}
}

/* Now create the output data array and initialize the 
 *  first column with the class sizes.
 */
String[][] outData = new String[numClasses][numFiles + 1]
int fileNum = 0
for (String str : dataFiles) {
	startClassRead = false
	int classNum = 0
	new File(str).eachLine { line ->
		def cellData = line.split(delimiter)
		if (line.contains("File name:")) {
			fileNames[fileNum] = cellData[1].trim()
		}
		if (startClassRead) {
			if (!cellData[0].trim().isEmpty()) {
				if (fileNum == 0) {
					outData[classNum][0] = cellData[0].trim()
				}
				if (cellData.size() == 2) {
					outData[classNum][fileNum + 1] = cellData[1].trim()
				} else {
					outData[classNum][fileNum + 1] = "0"
				}
				classNum++
			}
		}
		if (line.contains("um${delimiter}Volume")) {
			startClassRead = true
		}
	}
	fileNum++
}

// create and print the header row
def str = "Size Listing um"
for (String fileName : fileNames) {
	str += "\t" + fileName
}
println(str)

// now print the rows of data tab delimited
for (int i = 0; i < numClasses; i++) {
	str = outData[i][0]
	for (int j = 1; j < outData[i].length; j++) {
		str += "\t" + outData[i][j]
	}
	println(str)
}
