def executabledir = "/Users/johnlindsay/Projects/whitebox/bin/"
def executablestr = "./gospatial"
def arg1 = " -cwd=\"/Users/johnlindsay/Documents/Research/FastBreaching/data/\""
def arg2 = " -run=\"filldepressions\""
def arg3 = " -args=\"quebec DEM.dep;tmp11.dep;true\""
def cmd = executablestr + arg1 + arg2 + arg3
 
def executeOnShell(String command, File workingDir) {
  //println command
  def process = new ProcessBuilder(addShellPrefix(command))
                                    .directory(workingDir)
                                    .redirectErrorStream(true) 
                                    .start()
  process.inputStream.eachLine {
  		str = ((String)it).trim()
  		if (!str.isEmpty()) {
  			if (str.contains("%")) {
  				def strArray = str.split(" ")
  				String label = str.replace(strArray[strArray.length-1], "")
  				int progress = Integer.parseInt(strArray[strArray.length-1].replace("%", "").trim())
				pluginHost.updateProgress(label, progress)	  				
	  		} else {
	  			println(str)
	  		}
  		}
  }
  process.waitFor();
  
  return process.exitValue()
}
 
def addShellPrefix(String command) {
  commandArray = new String[3]
  commandArray[0] = "sh"
  commandArray[1] = "-c"
  commandArray[2] = command
  return commandArray
}

executeOnShell(cmd, new File(executabledir))

pluginHost.returnData("/Users/johnlindsay/Documents/Research/FastBreaching/data/tmp11.dep")
pluginHost.updateProgress(-1)
