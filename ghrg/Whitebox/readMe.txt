Welcome to Whitebox GAT!

Whitebox is a Java program. As such, its executable file is a .jar file rather than the usual .exe (on Windows, at least). Double-clicking the .jar file that you see in this folder (e.g. WhiteboxGAT_2_0_3.jar) will launch the actual program JAR, contained in the 'lib' folder (WhiteboxGIS.jar) using the default parameters that are set in the 'appstart.properties' file, also contained within this folder. You're system may have to be set-up to run .jar files, or permission may have to be granted to the file to launch, depending on your system set-up.

Importantly, if your computer has more memory (RAM) than the specified 1GB within this properties file, you should increase this default setting to allow Whitebox to better take advantage of the system resources. This can be done by changing the line:

app.vm.options=-Xmx1g

to something like:

app.vm.options=-Xmx8g

In the example above, I have reset the default to assume the java virtual machine (JVM) on which Whitebox runs can use 8GB of memory. If this value is higher than the available system memory, Whitebox will likely not launch and you will need to lower this setting to an appropriate level. Notice that it will likely still work if you double-click the WhiteboxGIS.jar within the 'lib' folder.

If you double-click the Whitebox .jar file and nothing happens, after having set the parameters above, it is likely because the file association for jar files has been hijacked. This can be easily resolved on most version of Windows, but unfortunately the required software has been removed on Windows 7 (typical of MS!). Try downloading Jarfix 2.0.0 and follow the instructions on this site:

http://johann.loefflmann.net/en/software/jarfix/index.html
