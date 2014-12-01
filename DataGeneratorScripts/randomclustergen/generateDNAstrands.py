import sys
import csv
import numpy
import getopt
import math
import difflib

# Either clusters/points per cluster or total number of points must be specified
def usage():
    print '$> python generateDNAstrands.py <required args> [optional args]\n' + \
        '\t-o <file>\tFilename for the output of the dna strands\n' + \
        '\t-s <#>\t\tSize of DNA strand\n' + \
        '\t-n [#]\t\tTotal number of points\n' + \
        '\t-c [#]\t\tNumber of clusters to generate\n' + \
        '\t-p [#]\t\tNumber of points per cluster\n' 

def stringDistance(p1, p2):
    '''
    Takes two DNA strands and computes the distance between them.
    '''
    result = len(p1)         
    edits = 0
    difference = difflib.ndiff(p1, p2)
    for i,s in enumerate(difference):
        if s[0]=='-' or s[0]=='+':
          edits = edits + 1
    edits = edits/2
    return result - edits


def tooClose(point, points, minDist):
    '''
    Computes the edit distance between the point and all points
    in the list, and if any points in the list are closer than minDist,
    this method returns true.
    '''
    for pair in points:
        if euclideanDistance(point, pair) < minDist:
                return True

    return False

def handleArgs(args):
    # set up return values
    sizeStrand = -1
    numClusters = -1
    numPoints = -1
    totPoints = -1
    output = None

    try:
        optlist, args = getopt.getopt(args[1:], 'o:s:n:c:p:')
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    for key, val in optlist:
        # first, the required arguments
        if key == '-o':
            output = val
        elif key == '-s':
            sizeStrand = int(val)
        # now, the optional arguments
        elif key == '-n':
            totPoints = int(val)
        elif key == '-c':
            numClusters = int(val)
        elif key == '-p':
            numPoints = int(val)

    # check required arguments were inputted  
    if numClusters < 0 or numPoints < 0 or \
            totPoints < 0 or \
            output is None:
        usage()
        sys.exit()
    return (sizeStrand, numClusters, numPoints, \
			totPoints, output)

# probably a redundant method
def drawOrigin(maxValue):
    return numpy.random.uniform(0, maxValue, 2)

# get a random base
def getRandomBase():
	r = numpy.random.rand()
	if r < 0.25:
		return 'a'
	elif r >= 0.25 and r < 0.5:
		return 'g'
	elif r >= 0.5 and r < 0.75:
		return 'c'
	else:
		return 't'


# start by reading the command line
sizeStrand, \
numClusters, \
numPoints, \
totPoints, \
output = handleArgs(sys.argv)

writer = csv.writer(open(output, "w"))

# if total points is given => generate strands randomly all over the place
if totPoints > 0:
		for i in range(totPoints):
			s = ""
			# generate a random strand
			for j in range(sizeStrand):
				s += getRandomBase()
			writer.writerow(s)	

	
# else, generate points according to clusters - also store underlying 
# cluster and centroid info for future verification
else :
		# step 1: generate each 2D centroid
		centroids_radii = []
		minDistance = 0
		for i in range(0, numClusters):
			centroid_radius = drawOrigin(maxValue)
			# is it far enough from the others?
			while (tooClose(centroid_radius, centroids_radii, minDistance)):
				centroid_radius = drawOrigin(maxValue)
			centroids_radii.append(centroid_radius)

		# step 2: generate the points for each centroid
		points = []
		minClusterVar = 0
		maxClusterVar = 0.5
		for i in range(0, numClusters):
			# compute the variance for this cluster
			variance = numpy.random.uniform(minClusterVar, maxClusterVar)
			cluster = centroids_radii[i]
			for j in range(0, numPoints):
				# generate a 2D point with specified variance
				# point is normally-distributed around centroids[i]
				x, y = numpy.random.normal(cluster, variance)
				# write the points out
				writer.writerow([x, y])
