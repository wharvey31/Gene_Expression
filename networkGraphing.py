#!/usr/bin/env python

# This script is intended to compare relative resource distribution of multiple organisms given KEGG Pathway maps and count files
	# Hierarchy of connections: Genes -> Pathways -> SubNetworks -> Networks

import argparse
import pandas as pd
import json
import math
from matplotlib import pyplot as plt
import re

parser = argparse.ArgumentParser() 

# Two arguments accepted, -d accepts a directory file with locations of count and pathway files
	# -g is for only displaying specific levels of the KEGG organism map
parser.add_argument("--graphmode", "-g", type=int, required=False, default=99)
parser.add_argument("--directory", "-d", type=str, required=True)

args = parser.parse_args()

# This also performs minor text parsing to remove variance between datasets
	# KEGG pathway maps add organism specific tags if genes are contained in a pathway, network, etc. 
		# This step removes those specific tags so organisms can be compared
def labelClean(label, joinString):
	uniqueID = str(label).split(" ", 1)
	# Unique tags are contained within brackets, so a bracket test is performed
	if uniqueID[1][-1] == "]":
		cleanedLabel = label.split(" ")
		cleanedLabel.pop()
		cleanedLabel = joinString.join(cleanedLabel)
	# If there is no specific tag, the labels are returned with required replacements
	else:
		cleanedLabel = uniqueID[1].replace(' ', joinString)
	return str(cleanedLabel)

def dictMaker(traceBack, networkTrackDict, subNetworkTrackDict, pathwayTrackDict, locus):
	# traceBack contains the path to the gene, the following steps parse and extract that path
	networkLabel = labelClean(traceBack[1], "\n")
	subNetworkLabel = labelClean(traceBack[2], " ")
	pathwayLabel = labelClean(traceBack[3], " ")
	# Checks for existence of specific networks/subNetworks/pathways and generates a list of subNetworks/pathways/genes contained within
	if str(networkLabel) in networkTrackDict:
		networkTrackDict[networkLabel].append(subNetworkLabel)
	else:
		networkTrackDict[networkLabel] = [subNetworkLabel]
	if str(subNetworkLabel) in subNetworkTrackDict:
		subNetworkTrackDict[subNetworkLabel].append(pathwayLabel)
	else:
		subNetworkTrackDict[subNetworkLabel] = [pathwayLabel]
	if str(pathwayLabel) in pathwayTrackDict:
		pathwayTrackDict[pathwayLabel].append(locus)
	else:
		pathwayTrackDict[pathwayLabel] = [locus]
	return networkTrackDict, subNetworkTrackDict, pathwayTrackDict


def addition(networkDict, subNetworkDict, pathwayDict, refDict, network, subNetwork, path, timepoint, tag, expression):
	# Tracks total expression at this timepoint
	expression += float(refDict[tag][timepoint])
	# Adds matrix counts to the appropriate dictionary
	if network in networkDict:
		networkDict[network] += float(refDict[tag][timepoint])
	else:
		networkDict[network] = float(refDict[tag][timepoint])
	if subNetwork in subNetworkDict:
		subNetworkDict[subNetwork] += float(refDict[tag][timepoint])
	else:
		subNetworkDict[subNetwork] = float(refDict[tag][timepoint])
	if path in pathwayDict:
		pathwayDict[path] += float(refDict[tag][timepoint])
	else:
		pathwayDict[path] = float(refDict[tag][timepoint])
	return expression, networkDict, subNetworkDict, pathwayDict

#Graphing function
def displayGraphs(countDict = None, timepoint = 0, figSize = (10,6), labSize = 5, yPlotLabel = "Relative Read Distribution as Percent of Total",
						xPlotLabel = "Pathway", rotate_legend = False, bottom_adj = None):
	labelList = []
	countList = []
	finalDict = {}
	# Adds network/subNEtwork/pathway labels to a list for future sorting 
	for key in countDict:		
		labelList.append(key)
	for i in range(len(organismList)):
		countList.append([])
		for j in range(len(labelList)):
			countList[i].append(countDict[labelList[j]][i])
	for i in range(len(organismList)):
		# List comprehension to ensure sorted order of results
		countList[i] = [x for j, x in sorted(zip(labelList, countList[i]))]
		finalDict[organismList[i]] = countList[i]
	labelList.sort()
	# Final dict contains the proper data structure for pandas graphing
	df = pd.DataFrame(finalDict, index=labelList)
	if rotate_legend == True:
		# Rotates x-axis labels
		ax = df.plot.bar(rot=90, style='ggplot',figsize=figSize)
	else: 
		ax = df.plot.bar(rot=0, style='ggplot',figsize=figSize)
	# Anchors legend to outer edge of graph
	plt.legend(bbox_to_anchor=(1., 1.), loc=2, borderaxespad=0., title="Organism")
	ax.xaxis.set_tick_params(labelsize=labSize)
	plt.subplots_adjust(bottom=bottom_adj)
	plt.title("Resource Distribution at "+str(timepointList[timepoint]))
	plt.ylabel(yPlotLabel)
	plt.xlabel(xPlotLabel, labelpad=20)
	plt.show()
	
# Recursive method for searching through pathway map given dictionary generated from JSON and an NCBI locus ID
def all_keys(search_dict, key_id):
	# Calls itself (recursion)
	def _all_keys(search_dict, key_id, keys=None):
		# Establishes key list
		if not keys:
			keys = []
		# Searches through specific level of dictionary
		for i in search_dict:
			# Looks for a leaf which matches the NCBI locus ID
			if re.match(str(key_id),str(search_dict[i])):
				# Returns the current list concatenated with the 'name' label for desired leaf
				return keys+[i]
			if isinstance(search_dict[i], dict):
				# If the dictionary is nested, grab all possible dictionaries and add them to a key search list
				potential_keys = _all_keys(search_dict[i], key_id, keys + [search_dict['name']])
				if 'name' in potential_keys:
					keys = potential_keys
					break
			if isinstance(search_dict[i], list):
				# If the dictionary reveals a list, iterate over the list
				for item in search_dict[i]:
					potential_keys = _all_keys(item, key_id, keys + [search_dict['name']])
					if 'name' in potential_keys:
						keys = potential_keys
						break
		return keys
	return _all_keys(search_dict, key_id)[:-1]

# Final function to check if a given dictionary contains values for at least one KEGG category
def cleanDict(checkDict):
	for i in checkDict.keys():
		if float(sum(checkDict[i])) == 0.0:
			del checkDict[i]
	return checkDict

# Creates lists to hold the respective json and count files along with additional basic information
jsonFileList = []
countFileList = []
organismList = []
geneList = []
tagDict = {}
minTimepoints = 999
# Standard List of sampling times
timepointList = ["Early Log Phase", "Middle Log Phase", "Late Log Phase", "Starving Phase", "Late Starving Phase"]
data = []

# Opens the directory file and parses through it accordingly
with open(args.directory) as d:
	for line in d:
		line = line.rstrip()
		line = line.split("\t")
		# JSON files should always be in the first position
		jsonFileList.append(line[0])
		# Count files in the second
		countFileList.append(line[1])

# Reads in gene count for each file
for i in range(len(countFileList)):
	geneList.append([])
	with open(countFileList[i]) as f:
		next(f)
		for line in f:
			line = line.rstrip()
			line = line.split("\t")
			# Takes first column of matrix as locus IDs
			tag = line.pop(0)
			# adds count values to dictionary holding all locus IDs
			tagDict[tag] = line
			geneList[i].append(tag)
		time = len(line)
		# Checks for minimum dimensionality of count data
		if time < minTimepoints:
			minTimepoints = time

# Creates unique file names based off of KEGG organism IDs and cycles through JSONs
for file in jsonFileList:
	filename = file.split("/")
	filename = filename[-1]
	filename = filename[0:3]
	organismList.append(filename)
	# Opens JSON
	with open(file) as f:
		data.append(json.load(f))

# Removes Brite Hierarchies categories (redundant/uninformative)
for listing in data:
	for category in listing['children']:
		if str(category['name']) == '09180 Brite Hierarchies':
			listing['children'].pop(listing['children'].index(category))

# Basic holding structure for respective levels of KEGG pathway
networkTrackDict = {}
subNetworkTrackDict = {}
pathwayTrackDict = {}

for organism in range(len(geneList)):
	print "Establishing pathway map for "+str(organismList[organism])+". . ."
	for locus in geneList[organism]:
		# Calls recursive function to connect genes to pathways, etc.
		traceBack = all_keys(data[organism], locus)
		# Checks to see if a hit is found. Hits will always be of full length
		if len(traceBack) == 4:
			networkTrackDict, subNetworkTrackDict, pathwayTrackDict = dictMaker(traceBack, networkTrackDict, subNetworkTrackDict, pathwayTrackDict, locus)
		else:
			continue

# Cycles through timepoints and organisms to add counts per pathway
for timepoint in range(minTimepoints):
	# Initialize a combine mapping dictionary for wach level which is used to speed up graphing functionality
	combinedNetworkDict = {}
	for network in networkTrackDict:
		combinedNetworkDict[network] = [0.0] * len(organismList)
	combinedSubNetworkDict = {}
	for subNetwork in subNetworkTrackDict:
		combinedSubNetworkDict[subNetwork] = [0.0] * len(organismList)
	combinedPathwayDict = {}
	for pathway in pathwayTrackDict:
		combinedPathwayDict[pathway] = [0.0] * len(organismList)
	for organism in range(len(geneList)):
		# Tracking for total expression per organism per timepoint
		expression_total = 0
		# Establishes individual dictionaries containing counts for each organism and KEGG network/subNetwork/pathway
		networkCountDict = {}
		subNetworkCountDict = {}
		pathwayCountDict = {}
		# Cycles through genes identified in count sheet
		for locus in geneList[organism]:
		# Checks to see if those genes are identified in any pathways
			for path in pathwayTrackDict:
				# If they are found, then check which category the pathway falls into
				if str(locus) in pathwayTrackDict[path]:
					for subNetwork in subNetworkTrackDict:
						# If it is found in a subNetwork, continue to nook for network
						if str(path) in subNetworkTrackDict[subNetwork]:
							for network in networkTrackDict:
								# Checks for final match to a network
								if str(subNetwork) in networkTrackDict[str(network)]:
									# Sends network name, pathway name, timpoint, and gene locus to addition function 
										# returns an updated expression total for the timepoint, network dictionary, subNetwork dictionary, and pathway dictionary
									expression_total, networkCountDict, subNetworkCountDict, pathwayCountDict = \
										addition(networkCountDict, subNetworkCountDict, pathwayCountDict, tagDict,
											str(network), str(subNetwork), str(path), timepoint, locus, expression_total)			
		# Adds counts to a shared dictionary with values being n-dimensional lists where n is number of organisms examined
		for key in networkCountDict:
			combinedNetworkDict[key][organism] = networkCountDict[key]/expression_total
		for key in subNetworkCountDict:
			combinedSubNetworkDict[key][organism] = subNetworkCountDict[key]/expression_total
		for key in pathwayCountDict:
			combinedPathwayDict[key][organism] = pathwayCountDict[key]/expression_total
	# graphing function for networks only
	if args.graphmode == 1:
		print "Cleaning graph output for "+str(timepointList[timepoint]+". . .")
		combinedNetworkDict = cleanDict(combinedNetworkDict)
		displayGraphs(countDict=combinedNetworkDict, timepoint=timepoint, figSize=(10,6))
	# graphing function for subnetworks only
	elif args.graphmode == 2:
		print "Cleaning graph output for "+str(timepointList[timepoint]+". . .")
		combinedSubNetworkDict = cleanDict(combinedSubNetworkDict)
		displayGraphs(countDict=combinedSubNetworkDict, timepoint=timepoint, figSize=(15,6), rotate_legend=True, bottom_adj=0.31)
	# graphing function for pathways only
	elif args.graphmode == 3:
		print "Cleaning graph output for "+str(timepointList[timepoint]+". . .")
		combinedPathwayDict  = cleanDict(combinedPathwayDict)
		displayGraphs(countDict=combinedPathwayDict, timepoint=timepoint, figSize=(15,6), rotate_legend=True, bottom_adj=0.40)
	# graphing function for networks, subnetworks, and pathways (default)
	else:
		print "Cleaning graph output for "+str(timepointList[timepoint]+". . .")
		combinedNetworkDict = cleanDict(combinedNetworkDict)
		displayGraphs(countDict=combinedNetworkDict, timepoint=timepoint, figSize=(10,6))
		combinedSubNetworkDict = cleanDict(combinedSubNetworkDict)
		displayGraphs(countDict=combinedSubNetworkDict, timepoint=timepoint, figSize=(15,6), rotate_legend=True, bottom_adj=0.31)
		combinedPathwayDict  = cleanDict(combinedPathwayDict)
		displayGraphs(countDict=combinedPathwayDict, timepoint=timepoint, figSize=(15,6), rotate_legend=True, bottom_adj=0.40)