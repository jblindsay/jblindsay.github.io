from math import atan2
from math import pi
from math import degrees

def convertToAzimuth(x1, y1, x2, y2):
	# Calculate the clockwise angle between points relative to north
	theta = -atan2(y2 - y1, x2 - x1) + pi / 2
	# The above statement is in the range -180 to 180,
	# convert this so it is 0 to 360 
	if theta < 0:
		theta += 2 * pi
	# Convert radians to degrees
	azimuth = degrees(theta)
	return azimuth

x1 = -81328.998084106
y1 = 7474929.8690234
x2 = 4125765.0381464
y2 = 7474929.8690234
print(convertToAzimuth(x1, y1, x2, y2))