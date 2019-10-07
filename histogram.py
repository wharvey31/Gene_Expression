#!/usr/bin/env python

import argparse
import pandas as pd
import json
import math
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser() 

parser.add_argument("--directory", "-d", type=str, required=True, default=0)
parser.add_argument("--bins", "-b", type=int, required=False, default=100)

args = parser.parse_args()

dfHolder = []
plotTitle = []
cutoff = 0.9

with open(args.directory) as f:
	for line in f:
		line = line.rstrip()
		title = line.split("/")
		plotTitle.append(title[1])
		with open(line) as g:
			dfHolder.append(pd.read_table(line, names=["locus", "ELP", "MLP", "LLP", "LSP", "SP"]))


for i in range(len(dfHolder)):
	dfHolder[i] = dfHolder[i][dfHolder[i].ELP < dfHolder[i].ELP.quantile(cutoff)]
	dfHolder[i] = dfHolder[i][dfHolder[i].MLP < dfHolder[i].MLP.quantile(cutoff)]
	dfHolder[i] = dfHolder[i][dfHolder[i].LLP < dfHolder[i].LLP.quantile(cutoff)]
	dfHolder[i] = dfHolder[i][dfHolder[i].LSP < dfHolder[i].LSP.quantile(cutoff)]
	dfHolder[i] = dfHolder[i][dfHolder[i].SP < dfHolder[i].SP.quantile(cutoff)]
	dfHolder[i] = dfHolder[i][["locus", "ELP", "MLP", "LLP", "LSP", "SP"]]
	dfHolder[i].hist(bins=args.bins)
	plt.subplots_adjust(hspace=0.5)
	plt.suptitle("Histogram for "+str(plotTitle[i])+" - "+str(args.directory))
	plt.show()