#!/usr/bin/env python

# Author: Jason Buenrostro, Stanford University 
# modified from plotV_vC.py (Alicia)

# Will make a V-plot from bed regions

##### IMPORT MODULES #####
# import necessary for python
import os
import sys
import numpy as np
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Pool
import argparse
import random

#### OPTIONS ####
# read options from the command line
parser = argparse.ArgumentParser(description="Description of your program")
parser.add_argument("-a", help="<Reads> Accepts sorted BAM file")
parser.add_argument("-b", help="<Bed> Accepts bed file")
parser.add_argument("-o", help="OutputFile")
parser.add_argument("-e", default="1000", help="number of bases to extend to each side, default=1000")
parser.add_argument("-p", default="center", help="options:center,ends, default=center")
parser.add_argument("-c", default="20", help="number of threads to use, default=20")
parser.add_argument("-s", default='4', help="column in which strand information is located (1 being first), default=4")
parser.add_argument("-u", action="store_true", default=False, help="Print uncompressed output file")
parser.add_argument("-v", action="store_true", default=False, help="Print profile around bed file")
parser.add_argument("-i", action="store_true", default=False, help="Print insert sizes across intervals")
parser.add_argument("--window", default='20', help="window size for plotting")
args = parser.parse_args()

# return usage information if no arguments given
if not any(vars(args).values()):
    parser.print_help()
    sys.exit()

##### DEFINE FUNCTIONS #####
# assign mat
def asn_mat(val, mat, s_int, e_int, t, i, weight):
    if float(val) >= s_int and float(val) < e_int - 1 and t < rows:  # -1 correct?
        base = val - s_int
        if len(p1_ints[0]) == 3:
            mat[t][int(base)] += weight
        elif p1_ints[i][int(args.s) - 1] == "-":
            mat[t][len(mat[0]) - int(base) - 1] += weight
        else:
            mat[t][int(base)] += weight
    return mat

#compute Vplot Matrix for a particular chunk
def sub_Mat(start):
    # initialize data matrix
    mat = np.zeros([rows, cols])
    # loop through the intervals and get relevant info
    bamfile = pysam.AlignmentFile(args.a, "rb")
    end = min(start + chunksize, len(p1_ints))
    for i in range(start, end):
        # get interval as num
        center = int(p1_ints[i][1]) + (int(p1_ints[i][2]) - int(p1_ints[i][1])) / 2
        s_int = center - int(args.e)
        e_int = center + int(args.e)
        # loop through rds
        for p2_rds in bamfile.fetch(p1_ints[i][0], max(0, s_int - 2000), e_int + 2000):
            # check mapping quality
            if p2_rds.mapping_quality < 30:
                continue
            # get read positions
            if p2_rds.is_reverse:
                continue
            else:
                l_pos = p2_rds.reference_start + 4
                # calculate center point
                ilen = abs(p2_rds.template_length) - 9
                r_pos = l_pos + ilen
                c_pos = l_pos + ilen / 2
                if ilen % 2 == 1 and args.p == 'center':
                    mat = asn_mat(c_pos, mat, s_int, e_int, ilen, i, 0.5)
                    mat = asn_mat(c_pos + 1, mat, s_int, e_int, ilen, i, 0.5)
                elif ilen % 2 != 1 and args.p == 'center':
                    mat = asn_mat(c_pos, mat, s_int, e_int, ilen, i, 1)
                # save ends or read centers to v-plot
                elif args.p == 'ends':
                    mat = asn_mat(l_pos, mat, s_int, e_int, ilen, i, 1)
                    mat = asn_mat(r_pos, mat, s_int, e_int, ilen, i, 1)
                else:
                    sys.exit('Error, check parameters')
    return mat

##### INPUTS AND OUTPUTS #####
# Get intervals
p1_ints = np.loadtxt(args.b, dtype='str')

# Determine number of rows and columns for the matrix
rows = 1000
cols = int(args.e) * 2

# Split bedfile into chunks
maxi = len(p1_ints)
chunksize = maxi // int(args.c)
chunks = maxi // chunksize
starts = list(range(0, maxi, chunksize))

def sub_Mat(start):
    # Initialize data matrix
    mat = np.zeros([rows, cols])
    # Open and read BAM file
    with pysam.AlignmentFile(args.a, "rb") as bamfile:
        end = min(start + chunksize, len(p1_ints))
        for i in range(start, end):
            # Get interval as num
            center = int(p1_ints[i][1]) + (int(p1_ints[i][2]) - int(p1_ints[i][1])) / 2
            s_int = center - int(args.e)
            e_int = center + int(args.e)
            # Loop through rds
            for p2_rds in bamfile.fetch(p1_ints[i][0], max(0, s_int - 2000), e_int + 2000):
                # Check mapping quality
                if p2_rds.mapping_quality < 30:
                    continue
                # Get read positions
                if p2_rds.is_reverse:
                    continue
                else:
                    l_pos = p2_rds.reference_start + 4
                    # Calculate center point
                    ilen = abs(p2_rds.template_length) - 9
                    r_pos = l_pos + ilen
                    c_pos = l_pos + ilen / 2
                    if ilen % 2 == 1 and args.p == 'center':
                        mat = asn_mat(c_pos, mat, s_int, e_int, ilen, i, 0.5)
                        mat = asn_mat(c_pos + 1, mat, s_int, e_int, ilen, i, 0.5)
                    elif ilen % 2 != 1 and args.p == 'center':
                        mat = asn_mat(c_pos, mat, s_int, e_int, ilen, i, 1)
                    # Save ends or read centers to v-plot
                    elif args.p == 'ends':
                        mat = asn_mat(l_pos, mat, s_int, e_int, ilen, i, 1)
                        mat = asn_mat(r_pos, mat, s_int, e_int, ilen, i, 1)
                    else:
                        sys.exit('Error, check parameters')
    return mat

if __name__ == "__main__":
    # Parallel processed computation of matrix for each chunk
    with Pool(processes=int(args.c)) as pool:
#    pool = Pool(processes=int(args.c))
        sub_mats = pool.map(sub_Mat, starts, 1)

    # Sum up matrices for each chunk into matrix for all
    mat = np.zeros([rows, cols])
    for i in range(len(starts)):
        mat = mat + sub_mats[i]

    # Get column vector
    if args.v:
        mat = np.sum(mat, 0)
    if args.i:
        mat = np.sum(mat, 1)

    # Save matrix
    if not args.o:
        n1 = os.path.basename(args.a)
        n2 = os.path.basename(args.b)
        if args.v:
            args.o = n1 + '.' + n2 + '.vect'
        elif args.i:
            args.o = n1 + '.' + n2 + '.iSize'
        else:
            args.o = n1 + '.' + n2 + '.vplot'
    if args.u:
        np.savetxt(args.o, mat, delimiter='\t', fmt='%s')
    else:
        np.save(args.o, mat)

    # Plot
    fig = plt.figure(figsize=(8.0, 5.0))
    xran = min(500, int(args.e))
    yran = min(500, rows)
    if args.v:
        plt.plot(mat / np.mean(mat[1:200]), 'k.')
        plt.plot(np.convolve(mat, np.ones(int(args.window)), 'same') / int(args.window) / np.mean(mat[1:200]), 'r')
        plt.xlabel('Position relative to center')
        plt.ylabel('Insertions')
    elif args.i:
        plt.plot(mat[0:990])
        plt.xlabel('Insert size (bp)')
        plt.ylabel('Count')
    else:
        plt.imshow(mat[0:yran, (int(args.e) - xran):(int(args.e) + xran + 1)], origin='lower', aspect='equal', extent=[-xran, xran, 1, yran + 1])
        plt.xlabel('Position relative to center')
        plt.ylabel('Insert size')

    # Save figure
    fig.savefig(args.o + '.png')
    plt.close(fig)
