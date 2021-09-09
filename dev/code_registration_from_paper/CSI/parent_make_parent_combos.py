#!/usr/bin/env python
import sys
import pandas as pd
import itertools

# set parameters
path = sys.argv[1]
maxNumParents = int(sys.argv[2])  # the maximum number of genes to consider in
                                  # combination as potential regulators of
                                  # "target".
target = sys.argv[3]  # id of target gene


# load gene expression data  prepared by "parent_prepare_data.R"
df = pd.read_csv(path+'/express_data.csv')
names = list(df.columns.values)


# get all gene ids got expression data for - these are the candidate
# regulators of "target".
try:
    names.remove('sample_id')
    names.remove('group')
except ValueError:
    pass
# remove target gene from candidate parent combos
parents = names
parents = [p for p in parents if target not in p]


# generate all combinations of possible of regulators of "target"
combos = []
for i in range(1, maxNumParents+1):
    for combo in itertools.combinations(parents, i):
        combos.append(','.join(combo))

# return output to bash by printing
print(' '.join(combos))
