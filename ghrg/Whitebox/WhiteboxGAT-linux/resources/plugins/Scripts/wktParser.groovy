import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
//import com.vividsolutions.jts.io.WKTReader

def wd = pluginHost.getWorkingDirectory()
//def prjFileName = wd + "StatsCanada_2006_CartoProv.prj"
//def prjFileName = "/Users/johnlindsay/Documents/Data/TM_WORLD_BORDERS-0/TM_WORLD_BORDERS-0.2.prj"
//def prjFileName = wd + "StatsCan_2006_Ont.prj"
def prjFileName = "/Users/johnlindsay/Documents/Research/StreamBurning/Data/tmp1.prj"

//WKTReader wktReader = new WKTReader()

def prjFile = new File(prjFileName)

def getTabString = { i -> 
	StringBuilder sb = new StringBuilder()
	for (j in 0..<i) {
		sb.append("\t")
	}
	return sb.toString()
}

def parseWKT
parseWKT = { str, order ->
	//println str
	int c = str.indexOf(",")
	int s = str.indexOf("[")
	String orderStr = getTabString(order)
	if (c >= 0 && c < s) {
		if (c > 0) {
			def word = str.substring(0, c).trim()
			def content = str.substring(c + 1).trim()
			println "${orderStr}${word}"
			parseWKT(content, order)
		} else {
			content = str.substring(1)
			parseWKT(content, order)
		}
	} else if (s < c && s >= 0) {
		def word = str.substring(0, s).trim()
		// find the closing bracket
		char[] ch = str.getChars()
		int k = 0
		int e = -1
		for (int m in 0..<ch.length) {
			if (ch[m] == '[') {
				k++
			} else if (ch[m] == ']') {
				k--
				if (k == 0) {
					e = m
					break
				}
			}
		}
		if (e >= 0) {
			def content = str.substring(s + 1, e).trim()
			println "${orderStr}${word}"
			if (!content.contains("[")) {
				println "${orderStr}\t${content}"
			}
			parseWKT(content, order + 1)
			if (e < str.length() - 1) {
				def remainder = str.substring(e + 1)
				parseWKT(remainder, order)
			}
		}
	}
}

prjFile.eachLine() { str -> 
	//println str
	parseWKT(str, 0)
}
