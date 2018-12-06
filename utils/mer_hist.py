import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i','--infile',type=int,action='store',help='input file')
args = parser.parse_args()

with open(f'{args.infile}','r') as file:
	rawtext = file.readlines()

hist = np.array([[float(x) for x in line.split()] for line in rawtext])
plt.bar(hist[:,0],hist[:,1])
	
