// imports
var File = Java.type('java.io.File');
var ProcessBuilder = Java.type('java.lang.ProcessBuilder');


var executabledir = "/Users/johnlindsay/Projects/whitebox/bin/";
var executablestr = "./gospatial";
var arg1 = " -cwd=\"/Users/johnlindsay/Documents/Research/FastBreaching/data/\"";
var arg2 = " -run=\"filldepressions\"";
var arg3 = " -args=\"quebec DEM.dep;tmp11.dep;true\"";
var cmd = executablestr + arg1 + arg2 + arg3;
 
function executeOnShell(command, workingDir) {
  //println command
  var process = new ProcessBuilder(addShellPrefix(command))
                                    .directory(workingDir)
                                    .redirectErrorStream(true)
                                    .start()
//  process.inputStream.forEach(function(str) {
//  		if (!str.trim()) {
//  			if (str.contains("%")) {
//  				var strArray = str.split(" ")
//  				var label = str.replace(strArray[strArray.length-1], "")
//  				var progress = parseInt(strArray[strArray.length-1].replace("%", "").trim())
//				pluginHost.updateProgress(label, progress)	  				
//	  		} else {
//	  			print(str)
//	  		}
//  		}
//  });
  process.waitFor();
  
  return process.exitValue();
}

function addShellPrefix(command) {
  var commandArray = ["sh", "-c", command];
  return commandArray;
}

executeOnShell(cmd, new File(executabledir));

pluginHost.returnData("/Users/johnlindsay/Documents/Research/FastBreaching/data/tmp11.dep");
pluginHost.updateProgress(-1);
