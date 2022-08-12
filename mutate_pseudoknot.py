#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

# ./mutate_pseudoknot.py hcov_229e_pk.fasta 7-15 41-49

import sys
from itertools import product
import subprocess, shlex

import numpy as np
from scipy.stats import norm
from scipy.stats.mstats import zscore


command = f'pKiss --mode=shapes'.split()
nucleotides = ['A','C','G','U']

sequenceFile = sys.argv[1]
leftside = sys.argv[2] # * format: xxxxx-yyyyy
rightside = sys.argv[3] # * format: xxxxx-yyyyy


def predict_pk_shape(arg):
  p = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
  pkiss_output = p.communicate(input=f"{arg}".encode())[0]
  return pkiss_output.decode().split("\n")[2].split()


header = ""
seq = ""

with open(sequenceFile) as inputStream:
  for line in inputStream:
    if line.endswith("#"):
      continue
    if line.startswith(">"):  # parse header
      header = line.strip().split(" ")[0][1:]  # remove '>' and ' '
    else:
      seq += line.strip()
    

try:
  leftside = list(map(int, leftside.split('-')))
  rightside = list(map(int, rightside.split('-')))
except Exception:
  print(leftside, rightside)
  sys.exit(1)

left_stem = seq[leftside[0]:leftside[1]]
right_stem = seq[rightside[0]:rightside[1]]
assert len(left_stem) == len(right_stem)

WT_SHAPE = predict_pk_shape(seq)
mutatedSequences = product(nucleotides, repeat = 2*len(left_stem))

similar_shapes = [WT_SHAPE]

for mutant in mutatedSequences:
  left = ''.join(mutant[:len(mutant)//2])
  right = ''.join(mutant[len(mutant)//2:])
  mutated_pk = seq[:leftside[0]] + left  + seq[leftside[1]:rightside[0]]  + right + seq[rightside[1]:]
  if mutated_pk == seq: continue
  assert len(mutated_pk) == len(seq)
  mutantPK = predict_pk_shape(mutated_pk)
  if mutantPK[-1] == WT_SHAPE[-1]:
    similar_shapes.append(mutantPK)

mfes = list(map(lambda x: x[0], similar_shapes))
a = np.array(mfes).astype(float)
z = zscore(a)
pvalues = norm.cdf(-1*abs(z)) * 2

WT_z = z[0]
WT_pvalue = pvalues[0]

print(seq)
print("\n".join(WT_SHAPE))
print(WT_z)
print(WT_pvalue)

print(len(mfes))