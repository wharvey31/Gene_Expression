#!/usr/bin/env python

import argparse
import pandas as pd
import json
import math
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser() 

parser.add_argument("--proportionalCount", "-p", type=int, required=False, default=0)
parser.add_argument("--graphmode", "-g", type=int, required=False, default=99)
parser.add_argument("--directory", "-d", type=str, required=True)

args = parser.parse_args()

# Addition mechanism for tracking Overall Gene Network Categories (OGNC) and Pathway Counts
def addition(networkDict, subNetworkDict, pathwayDict, refDict, network, subnet, pathway, timepoint, tag, organism):
	expression = float(refDict[tag][timepoint])
	if args.proportionalCount == 1:
		if float(refDict[tag][timepoint]) > 0.0:
			# Values are percentage of total genes expressed at every timepoint
			networkDict[str(network)] += 1.0/float(totalGenes[organismList[organism]])
			pathwayDict[str(pathway)] += 1.0/float(totalGenes[organismList[organism]])
			subNetworkDict[str(subnet)] += 1.0/float(totalGenes[organismList[organism]])
	else:
		try:
			networkDict[str(network)] += float(refDict[tag][timepoint])
			pathwayDict[str(pathway)] += float(refDict[tag][timepoint])
			subNetworkDict[str(subnet)] += float(refDict[tag][timepoint])
		except:
			networkDict[str(network)] += 0.0
			pathwayDict[str(pathway)] += 0.0
			subNetworkDict[str(subnet)] += 0.0
	return expression

#Graphing function
def displayGraphs(countDict = None, timepoint = 0, figSize = (10,6), labSize = 5, yPlotLabel = "Total Read Count",
		 				xPlotLabel = "Pathway/Network", rotate_legend= False, bottom_adj = None):
	labelList = []
	countList = []
	finalDict = {}
	# Adds Pathway or Network
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
	df = pd.DataFrame(finalDict, index=labelList)
	if rotate_legend== True:
		# Sets up labels for condensed legend
		ax = df.plot.bar(rot=90, style='ggplot',figsize=figSize)
	else: 
		ax = df.plot.bar(rot=0, style='ggplot',figsize=figSize)
	plt.legend(bbox_to_anchor=(1., 1.), loc=2, borderaxespad=0., title="Organism")
	ax.xaxis.set_tick_params(labelsize=labSize)
	plt.subplots_adjust(bottom=bottom_adj)
	plt.title("Resource Distribution at "+str(timepointList[timepoint]))
	plt.ylabel(yPlotLabel)
	plt.xlabel(xPlotLabel, labelpad=20)
	plt.show()
	

def mergeDict(key, initDict, combDict, total):
	if key in combDict:
		combDict[key].append(float(initDict[key])/float(total))
		if (float(combDict[key][0]) + float(combDict[key][1])) == 0:
			del combDict[key]
	else:
		combDict[key] = [float(initDict[key])/(total)]
	return combDict


jsonFileList = []
countFileList = []
totalGenes = {}

with open(args.directory) as d:
	for line in d:
		line = line.rstrip()
		line = line.split("\t")
		jsonFileList.append(line[0])
		countFileList.append(line[1])

# Establishes dictionaries for each 
networkTrackDict = {}
subNetworkTrackDict = {}
pathwayTrackDict = {}
organismList = []
geneList = []
tagDict = {}

# Standard List of sampling times
timepointList = ["Early Log Phase", "Middle Log Phase", "Late Log Phase", "Starving Phase", "Late Starving Phase"]

minTimepoints = 999

# Reads in gene count for each file
for file in range(len(countFileList)):
	geneList.append([])
	with open(countFileList[file]) as f:
		next(f)
		for line in f:
			line = line.rstrip()
			line = line.split("\t")
			tag = line.pop(0)
			tagDict[tag] = line
			geneList[file].append(tag)
		time = len(line)
		# Checks for minimum dimensionality of count data
		if time < minTimepoints:
			minTimepoints = time

# Creates unique file names based off of KEGG organism IDs and cycles through JSONs
for file in range(len(jsonFileList)):
	filename = jsonFileList[file].split("/")
	filename = filename[-1]
	filename = filename[0:3]
	organismList.append(filename)
	# Opens JSON
	with open(jsonFileList[file]) as f:
		data = json.load(f)
		networks = data['children']
		for j in range(len(networks)):
			networkCat = networks[j]['name'].split(" ",1)
			networkLabel = str(networkCat[1]).replace(' ', '\n')
			if networkLabel[-1] == "]":
				networkLabel = networkLabel.split("\n")
				networkLabel.pop()
				networkLabel = "\n".join(networkLabel)
			if str(networkLabel) in networkTrackDict:
				pass
			else:
				networkTrackDict[str(networkLabel)] = []
			for k in range(len(networks[j]['children'])):
				# Adds dictionary which adds pathways to Networks
				subNetworkCat = networks[j]['children'][k]['name'].split(" ",1)
				subNetworkLabel = str(subNetworkCat[1])
				if subNetworkLabel[-1] == "]":
					subNetworkLabel = subNetworkLabel.split(" ")
					subNetworkLabel.pop()
					subNetworkLabel = " ".join(subNetworkLabel)
				if str(subNetworkLabel) in networkTrackDict[networkLabel]:
					pass
				else:
					networkTrackDict[str(networkLabel)].append(str(subNetworkLabel))
				if str(subNetworkLabel) in subNetworkTrackDict:
					pass
				else:
					subNetworkTrackDict[str(subNetworkLabel)] = []
				for l in range(len(networks[j]['children'][k]['children'])):
					pathwayCat =  networks[j]['children'][k]['children'][l]['name'].split(" ",1)
					pathwayLabel = str(pathwayCat[1])
					if pathwayLabel[-1] == "]":
						pathwayLabel = pathwayLabel.split(" ")
						pathwayLabel.pop()
						pathwayLabel = " ".join(pathwayLabel)
					if str(pathwayLabel) in subNetworkTrackDict[subNetworkLabel]:
						pass
					else:
						subNetworkTrackDict[subNetworkLabel].append(str(pathwayLabel))
					if str(pathwayLabel) in pathwayTrackDict:
						pass
					else:
						pathwayTrackDict[str(pathwayLabel)] = []
					try:
						for m in range(len(networks[j]['children'][k]['children'][l]['children'])):
							try:
								# Links genes - specifically NCBI loci - to pathways
								gene = networks[j]['children'][k]['children'][l]['children'][m]['name'].split(" ", 1)
								if str(gene[0]) in pathwayTrackDict[str(pathwayLabel)]:
									pass
								else:
									if gene[0] in geneList[file]:
										pathwayTrackDict[str(pathwayLabel)].append(str(gene[0]))
							except:
								continue
					except:
						gene = networks[j]['children'][k]['children'][l]['name'].split(" ", 1)
						if str(pathwayLabel) in pathwayTrackDict:
							if str(gene[0]) in pathwayTrackDict[str(pathwayLabel)]:
								pass
							else:
								if gene[0] in geneList[file]:
									pathwayTrackDict[str(pathwayLabel)].append(str(gene[0]))
						else:
							if gene[0] in geneList[file]:
								pathwayTrackDict[str(pathwayLabel)] = [gene[0]]

# Removes Brite Hierarchies categories (redundant/uninformative)
networkTrackDict.pop('Brite\nHierarchies',None)

# Total gene counts based off of total expression profile
totalGenes = {}
for i in range(len(organismList)):
	totalGenes[organismList[i]] = len(geneList[i])

# Cycles through timepoints and organisms to add counts per pathway
for timepoint in range(0,minTimepoints):
	combinedNetworkDict = {}
	combinedPathwayDict = {}
	combinedSubNetworkDict = {}
	organismTracker = 0
	for organism in geneList:
		expression_total = 0
		# Sets counts to 0 for all categories along each timepoint
		pathwayCountDict = {}
		networkCountDict = {}
		subNetworkCountDict = {}
		for key in pathwayTrackDict:
			pathwayCountDict[key] = 0
		for key in networkTrackDict:
			networkCountDict[key] = 0
		for key in subNetworkTrackDict:
			subNetworkCountDict[key] = 0
		# Cycles through genes identified in count sheet
		for locus in organism:
		# Checks to see if those genes are identified in any pathways
			for path in pathwayTrackDict:
				# If they are found, then check which category the pathway falls into
				if str(locus) in pathwayTrackDict[path]:
					for subNetwork in subNetworkTrackDict:
						if str(path) in subNetworkTrackDict[subNetwork]:
							for network in networkTrackDict:
								if str(subNetwork) in networkTrackDict[str(network)]:
									# Sends OGNC name, pathway name, timpoint, and gene locus to addition function									
									expression_total += addition(networkCountDict, subNetworkCountDict, pathwayCountDict, tagDict, str(network), str(subNetwork), str(path), timepoint, locus, organismTracker)
		# Adds counts to a shared dictionary with values being n-dimensional lists
		# Adds counts to a shared dictionary with values being n-dimensional lists
		for key in networkCountDict:
			combinedNetworkDict = mergeDict(key, networkCountDict, combinedNetworkDict, expression_total)
		for key in subNetworkCountDict:
			combinedSubNetworkDict = mergeDict(key, subNetworkCountDict, combinedSubNetworkDict, expression_total)
		for key in pathwayCountDict:
			combinedPathwayDict = mergeDict(key, pathwayCountDict, combinedPathwayDict, expression_total)
		organismTracker += 1
	# graphing function for networks only
	if args.graphmode == 1:
		displayGraphs(countDict=combinedNetworkDict, timepoint=timepoint, figSize=(10,6))
	# graphing function for pathways only
	elif args.graphmode == 2:
		displayGraphs(countDict=combinedPathwayDict, timepoint=timepoint, figSize=(15,6), rotate_legend=True, bottom_adj=0.37)
	# graphing function for subnetworks only
	elif args.graphmode == 3:
		displayGraphs(countDict=combinedSubNetworkDict, timepoint=timepoint, figSize=(15,6), rotate_legend=True, bottom_adj=0.31)
	# graphing function for networks, subnetworks, and pathways (default)
	else:
		displayGraphs(countDict=combinedNetworkDict, timepoint=timepoint, figSize=(10,6))
		displayGraphs(countDict=combinedSubNetworkDict, timepoint=timepoint, figSize=(15,6), rotate_legend=True, bottom_adj=0.31)
		displayGraphs(countDict=combinedPathwayDict, timepoint=timepoint, figSize=(15,6), rotate_legend=True, bottom_adj=0.37)
