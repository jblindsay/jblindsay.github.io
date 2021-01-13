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
def directory = "/Users/johnlindsay/Documents/JackiesData/"
//def directory = "/Users/johnlindsay/Dropbox/BF2 CSV/" // my test directory...it contained only the one file you sent me.

// change these next two as needed
def fileExtension = '$ls'
def delimiter = "\t"

def dataFiles = FileUtilities.findAllFilesWithExtension(directory, fileExtension, true)
int numFiles = dataFiles.size()
String[][] outData = new String[numFiles][26]
int fileNum = 0
for (String str : dataFiles) {
	def obscuration1read = false 
	new File(str).eachLine { line -> 
		def cellData = line.split(delimiter)
		if (line.contains("File name:")) {
			outData[fileNum][0] = cellData[1].trim()
		} else if (line.contains("File ID:")) {
			outData[fileNum][1] = cellData[1].trim()
		} else if (line.contains("Sample ID:")) {
			outData[fileNum][2] = cellData[1].trim()
		} else if (line.contains("Operator:")) {
			outData[fileNum][3] = cellData[1].trim()
		} else if (line.contains("Run number:")) {
			outData[fileNum][4] = cellData[1].trim()
		} else if (line.contains("Start time:")) {
			outData[fileNum][5] = cellData[1].trim()
		} else if (line.contains("Obscuration:") && !obscuration1read) {
			outData[fileNum][6] = cellData[1].trim()
			obscuration1read = true
		} else if (line.contains("PIDS Obscur:")) {
			outData[fileNum][7] = cellData[1].trim()
		} else if (line.contains("Obscuration:") && obscuration1read) {
			outData[fileNum][8] = cellData[1].trim()
		} else if (line.contains("From")) {
			outData[fileNum][9] = cellData[1].trim()
		} else if (line.contains("To")) {
			outData[fileNum][10] = cellData[1].trim()
		} else if (line.contains("Volume")) {
			outData[fileNum][11] = cellData[1].trim()
		} else if (line.contains("Mean:")) {
			outData[fileNum][12] = cellData[1].trim()
		} else if (line.contains("Median:")) {
			outData[fileNum][13] = cellData[1].trim()
		} else if (line.contains("Mean/Median ratio:")) {
			outData[fileNum][14] = cellData[1].trim()
		} else if (line.contains("Mode:")) {
			outData[fileNum][15] = cellData[1].trim()
		} else if (line.contains("S.D.:")) {
			outData[fileNum][16] = cellData[1].trim()
		} else if (line.contains("Variance:")) {
			outData[fileNum][17] = cellData[1].trim()
		} else if (line.contains("C.V.:")) {
			outData[fileNum][18] = cellData[1].trim()
		} else if (line.contains("Skewness:")) {
			outData[fileNum][19] = cellData[1].trim()
		} else if (line.contains("Kurtosis:")) {
			outData[fileNum][20] = cellData[1].trim()
		} else if (line.contains("10")) {
			outData[fileNum][21] = cellData[1].trim()
		} else if (line.contains("25")) {
			outData[fileNum][22] = cellData[1].trim()
		} else if (line.contains("50")) {
			outData[fileNum][23] = cellData[1].trim()
		} else if (line.contains("75")) {
			outData[fileNum][24] = cellData[1].trim()
		} else if (line.contains("90")) {
			outData[fileNum][25] = cellData[1].trim()
		}
	}
	fileNum++
}

// create and print the header row
def str = "File name\tFile ID\tSample ID\tOperator\tRun number\tStart time\tObscuration\tPIDS Obscur\tObscuration\tFrom\tTo\tVolume\tMean\tMedian\tMean/Median ratio\tMode\tS.D.\tVariance\tC.V.\tSkewness\tKurtosis\t10\t25\t50\t75\t90"
println(str)

// now print the rows of data as tab delimited data
for (int i = 0; i < numFiles; i++) {
	str = outData[i][0]
	for (int j = 1; j < outData[i].length; j++) {
		str += "\t" + outData[i][j]
	}
	println(str)
}
