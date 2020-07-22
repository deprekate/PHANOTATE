#############################################################################
# Full Imports
from __future__ import division
import math
import random
import sys
from decimal import Decimal

class KMeans():
	#def __init__(self, n_clusters=2, init='random', n_init=10, max_iter=300, tol=Decimal("0.0001") ):
	def __init__(self, n_clusters=2, init='random', n_init=10, max_iter=300, tol=0.0001 ):
		self.n_clusters = n_clusters
		self.init = init
		self.n_init = n_init
		self.max_iter = max_iter
		self.tol = tol
		self.cluster_centers_ = []
		self.labels_ = []
		self.withinss_ = []

	def fit(self, X):
		points = []
		for p in X:
			points.append(Point(p, 'na'))
		self.best_clusters = iterative_kmeans(
					points,
					self.n_clusters,
					self.tol,
					self.n_init
					)
		for c in self.best_clusters:
			self.cluster_centers_.append(c.centroid)
			self.withinss_.append( c.getTotalDistance() / c.length() )
		for p in points:
			for i, c in enumerate(self.best_clusters):
				if c.has_point(p):
					self.labels_.append(i)
		return self


#############################################################################
# K-means Methods

def iterative_kmeans(points, num_clusters, cutoff, iteration_count):
	"""
	K-means isn't guaranteed to get the best answer the first time. It might
	get stuck in a "local minimum."

	Here we run kmeans() *iteration_count* times to increase the chance of
	getting a good answer.

	Returns the best set of clusters found.
	"""
	candidate_clusters = []
	errors = []
	for _ in range(iteration_count):
		clusters = kmeans(points, num_clusters, cutoff)
		error = calculateError(clusters)
		candidate_clusters.append(clusters)
		errors.append(error)

	highest_error = max(errors)
	lowest_error = min(errors)
	ind_of_lowest_error = errors.index(lowest_error)
	best_clusters = candidate_clusters[ind_of_lowest_error]
	return best_clusters

def kmeans(points, k, cutoff, initial_centroids=False):
	# Pick out k random points to use as our initial centroids
	if(not initial_centroids):
		initial_centroids = random.sample(points, k)

	# Create k clusters using those centroids
	# Note: Cluster takes lists, so we wrap each point in a list here.
	clusters = [Cluster([p]) for p in initial_centroids]

	# Loop through the dataset until the clusters stabilize
	loopCounter = 0
	while True:
		# Create a list of lists to hold the points in each cluster
		lists = [[] for _ in clusters]
		clusterCount = len(clusters)

		# Start counting loops
		loopCounter += 1
		# For every point in the dataset ...
		for p in points:
			# Get the distance between that point and the centroid of the first
			# cluster.
			smallest_distance = getDistance(p, clusters[0].centroid)

			# Set the cluster this point belongs to
			clusterIndex = 0

			# For the remainder of the clusters ...
			for i in range(1, clusterCount):
				# calculate the distance of that point to each other cluster's
				# centroid.
				distance = getDistance(p, clusters[i].centroid)
				# If it's closer to that cluster's centroid update what we
				# think the smallest distance is
				if distance < smallest_distance:
					smallest_distance = distance
					clusterIndex = i
			# After finding the cluster the smallest distance away
			# set the point to belong to that cluster
			lists[clusterIndex].append(p)

		# Set our biggest_shift to zero for this iteration
		biggest_shift = 0 #Decimal(0)

		# For each cluster ...
		for i in range(clusterCount):
			# Calculate how far the centroid moved in this iteration
			shift = clusters[i].update(lists[i])
			# Keep track of the largest move from all cluster centroid updates
			biggest_shift = max(biggest_shift, shift)

		# Remove empty clusters
		clusters = [c for c in clusters if len(c.points) != 0]

		# If the centroids have stopped moving much, say we're done!
		if biggest_shift < cutoff:
			#print "Converged after %s iterations" % loopCounter
			break
	return clusters


#############################################################################
# Classes

class Point(object):
	'''
	A point in n dimensional space
	'''
	def __init__(self, coords, name):
		'''
		coords - A list of values, one per dimension
		'''

		self.coords = coords
		self.n = len(coords)
		self.name = name

	def __repr__(self):
		return str(self.coords)

class Cluster(object):
	'''
	A set of points and their centroid
	'''

	def __init__(self, points):
		'''
		points - A list of point objects
		'''

		if len(points) == 0:
			raise Exception("ERROR: empty cluster")

		# The points that belong to this cluster
		self.points = points

		# The dimensionality of the points in this cluster
		self.n = points[0].n

		# Assert that all points are of the same dimensionality
		for p in points:
			if p.n != self.n:
				raise Exception("ERROR: inconsistent dimensions")

		# Set up the initial centroid (this is usually based off one point)
		self.centroid = self.calculateCentroid()

	def __repr__(self):
		'''
		String representation of this object
		'''
		return str(self.points)

	def has_point(self, point):
		'''
		Check if cluster has point
		'''
		return point in self.points

	def length(self):
		return len(self.points)

	def update(self, points):
		'''
		Returns the distance between the previous centroid and the new after
		recalculating and storing the new centroid.

		Note: Initially we expect centroids to shift around a lot and then
		gradually settle down.
		'''
		old_centroid = self.centroid
		self.points = points
		# Return early if we have no points, this cluster will get
		# cleaned up (removed) in the outer loop.
		if len(self.points) == 0:
			return 0

		self.centroid = self.calculateCentroid()
		shift = getDistance(old_centroid, self.centroid)
		return shift

	def calculateCentroid(self):
		'''
		Finds a virtual center point for a group of n-dimensional points
		'''
		numPoints = len(self.points)
		# Get a list of all coordinates in this cluster
		coords = [p.coords for p in self.points]
		# Reformat that so all x's are together, all y'z etc.
		unzipped = zip(*coords)
		# Calculate the mean for each dimension
		centroid_coords = [sum(dList)/numPoints for dList in unzipped]

		return Point(centroid_coords, "CENTROID")

	def getTotalDistance(self):
		'''
		Return the sum of all squared Euclidean distances between each point in 
		the cluster and the cluster's centroid.
		'''
		sumOfDistances = 0 #Decimal(0)
		for p in self.points:
			sumOfDistances += getDistance(p, self.centroid)

		return sumOfDistances
	
	def getMeanDistance(self):
		'''
		Return the mean of all squared Euclidean distances between each point in 
		the cluster and the cluster's centroid.
		'''
		Distances = []
		for p in self.points:
			Distances.append(getDistance(p, self.centroid))

		return sum(Distances) / len(Distances)

#############################################################################
# Helper Methods

def getDistance(a, b):
	'''
	Squared Euclidean distance between two n-dimensional points.
	https://en.wikipedia.org/wiki/Euclidean_distance#n_dimensions
	Note: This can be very slow and does not scale well
	'''
	if a.n != b.n:
		raise Exception("ERROR: non comparable points")

	accumulatedDifference = 0 #Decimal(0)
	for i in range(a.n):
		squareDifference = pow((a.coords[i]-b.coords[i]), 2)
		accumulatedDifference += squareDifference

	return accumulatedDifference

def makeRandomPoint(n, lower, upper):
	'''
	Returns a Point object with n dimensions and values between lower and
	upper in each of those dimensions
	'''
	p = Point([random.uniform(lower, upper) for _ in range(n)])
	return p

def calculateError(clusters):
	'''
	Return the average squared distance between each point and its cluster
	centroid.

	This is also known as the "distortion cost."
	'''
	accumulatedDistances = 0
	num_points = 0
	for cluster in clusters:
		num_points += len(cluster.points)
		accumulatedDistances += cluster.getTotalDistance()

	error = accumulatedDistances / num_points
	return error
