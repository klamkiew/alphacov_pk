#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

# ./mutate_pk_parallel.py hcov_229e_pk.fasta 7-15 41-49

import sys
from itertools import product, islice
import subprocess, shlex

#from multiprocessing import Pool
import multiprocessing as mp
import queue

import numpy as np
from scipy.stats import norm
from scipy.stats.mstats import zscore



command = f'pKiss --mode=shapes'.split()
nucleotides = ['A','C','G','U']
bpTypes = {
  'A' : ['U'],
  'C' : ['G'],
  'G' : ['C','U'],
  'U' : ['A','G']
}

sequenceFile = sys.argv[1]
leftside = sys.argv[2] # * format: xxxxx-yyyyy
rightside = sys.argv[3] # * format: xxxxx-yyyyy

output = sequenceFile.split('.')[0]

def predict_pk_shape(arg):
  """
  """
  p = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
  pkiss_output = p.communicate(input=f"{arg}".encode())[0]
  return pkiss_output.decode().split("\n")[2].split()


def prepare_pk(idx, mutant):
  """
  """
  left = ''.join(mutant[:len(mutant)//2])
  right = ''.join(mutant[len(mutant)//2:])
  mutated_pk = seq[:leftside[0]] + left  + seq[leftside[1]:rightside[0]]  + right + seq[rightside[1]:]
  if mutated_pk == seq: 
    return
  assert len(mutated_pk) == len(seq)
  mutantPK = predict_pk_shape(mutated_pk)
  if mutantPK[-1] == WT_SHAPE[-1]:
    with lock:
      outputStream.write(f">mutated_sequence_{idx}\n{mutated_pk}\n{' '.join(mutantPK)}\n")
      outputStream.flush()
      similar_shapes.append(float(mutantPK[0]))

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
similar_shapes = [float(WT_SHAPE[0])]

####################################################################
mutatedSequences = product(nucleotides, repeat = 2*len(left_stem))
#mutatedSequences = zip(mutatedSequences, range(4**(2*len(left_stem))))
manager = mp.Manager() # overpaid a**hole

similar_shapes = manager.list(similar_shapes)
lock = mp.Lock()
outputStream = open(f"{output}_background.out", 'w')
p = mp.Pool(8)
p.starmap(prepare_pk, enumerate(mutatedSequences))

outputStream.close()

mfes = similar_shapes
a = np.array(mfes).astype(float)
z = zscore(a)
pvalues = norm.cdf(-1*abs(z)) * 2

WT_z = z[0]
WT_pvalue = pvalues[0]

with open(f"{output}_eval.out", 'w') as outputStream:
  outputStream.write(f"{seq}\n")
  outputStream.write(f"{' '.join(WT_SHAPE)}\n")
  outputStream.write(f"{WT_z} {WT_pvalue}\n")
