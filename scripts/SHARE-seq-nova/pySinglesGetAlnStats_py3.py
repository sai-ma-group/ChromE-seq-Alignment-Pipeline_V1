#!/usr/bin/env python

# Author: Sai Ma, MSSM
# This will collect QC stats

##### IMPORT MODULES #####
# import necessary for python
import os
import re
import sys
import numpy as np
import subprocess
from argparse import ArgumentParser
import gzip
import pysam
import scipy.io as sio

##### DEFINE FUNCTIONS #####

#### OPTIONS ####
# read options from command line
parser = ArgumentParser()
parser.add_argument("-i", help="<SampleDataTable>")
parser.add_argument("-d", help="<Directory>")
parser.add_argument("-o", help="OutputFile")
args = parser.parse_args()

# return usage information if no argvs given
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()

##### INPUTS AND OUTPUTS #####
# get table to append to
table = np.loadtxt(args.i, dtype='str')

# initialize out
if not args.o:
    outName = os.path.basename(args.i).split('.')[0] + '.merged.xls'
else:
    outName = args.o

# dups log
cmd1 = 'find ' + args.d + ' -name "*.dups.log" | sort -V > ReadDups'
os.system(cmd1)
dupsFiles = np.loadtxt('ReadDups', dtype='str')


# get list of align stats data
cmd2 = 'find ' + args.d + ' -name "*.stats.log" | sort -V > ReadStats'
os.system(cmd2)
statFiles = np.loadtxt('ReadStats', dtype='str')

# get iSizes
cmd3 = 'find ' + args.d + ' -name "*.hist_data.log" | sort -V > ReadiSizes'
os.system(cmd3)
iSizeFiles = np.loadtxt('ReadiSizes', dtype='str')

# get enrich info
cmd4 = 'find ' + args.d + ' -name "*.RefSeqTSS" | sort -V > ReadTSS'
os.system(cmd4)
tssFiles = np.loadtxt('ReadTSS', dtype='str')

# get alignment info
cmd5 = 'find ' + args.d + ' -name "*.align.log" | sort -V > ReadAlignInfo'
os.system(cmd5)
alignFiles = np.loadtxt('ReadAlignInfo', dtype='str')

## set outputs
# duplicate file
frags = ['Frags']
alnRate = ['Align_Rate']
libSize = ['Lib_Size']
dupRate = ['Dup_Rate']
filtFrags = ['Filt_Frags']
alnFragsNoM = ['Aligned_noMT']

# align stats
chrMFrags = ['MT_frags']
chrMrate = ['MT_rate']
totalFiltRds = ['Final_frags']

# isize info
mediSize = ['Median_iSize']
meaniSize = ['Mean_iSize']
stdiSize = ['std_iSize']
widthiSize = ['80perc_Width_iSize']

# enrich info
tssEnrich = ['TSS_rate']

# isize matrix
iSizeMat = []

# pass vector
vector = [0]

# loop and cat data
for i in range(0, len(table) - 1):
    # name
    if len(table[i + 1][0]) == 1:
        name = table[i + 1]
    else:
        name = table[i + 1][0]

    # find
    if name in dupsFiles[i].split('/'):
        j = i
        print(name + ' ' + dupsFiles[j])
    else:
        for j in range(0, len(dupsFiles)):
            if name in dupsFiles[j].split('/'):
                print(name + ' ' + dupsFiles[j])
                break
            
    # check
    if name not in dupsFiles[j]:
        sys.exit('Check file: ' + dupsFiles[j])
    if name not in statFiles[j]:
        sys.exit('Check file: ' + statFiles[j])
    if name not in iSizeFiles[j]:
        sys.exit('Check file: ' + iSizeFiles[j])
    if name not in tssFiles[j]:
        sys.exit('Check file: ' + tssFiles[j])
    if name not in alignFiles[j]:
        sys.exit('Check file: ' + alignFiles[j])

    # skip if any problems
    try:
        iSizeHeaders = np.loadtxt(args.d + '/' + iSizeFiles[j], dtype='str', delimiter=',')[0].split('\t')
    except:
        continue
    
    # pass filter
    vector.append(i + 1)

    # get alignment info
    alignData = np.loadtxt(args.d + '/' + alignFiles[j], dtype='str', delimiter=',')
    for k in range(0, len(alignData)):
        if 'overall alignment rate' in alignData[k]:
            idx = k
    alnRate.append(float(alignData[idx].split('%')[0]) / 100.)
    for k in range(0, len(alignData)):
        if 'reads; of these:' in alignData[k]:
            idx = k
    frags.append(int(alignData[idx].split(' reads;')[0]))

    # get dups
    dupH = np.loadtxt(args.d + '/' + dupsFiles[j], dtype='str', delimiter=',')[0].split(' ')[1:]
#    print(dupH)
    try:
        dupV = np.array(np.loadtxt(args.d + '/' + dupsFiles[j], dtype='str', delimiter=',')[1].split(' ')[1:], float)
    except:
        dupV = np.append(np.array(np.loadtxt(args.d + '/' + dupsFiles[j], dtype='str', delimiter=',')[1].split(' ')[1:-1], float), 0)
    dupData = dict(zip(dupH, dupV))
    libSize.append(int(dupData['ESTIMATED_LIBRARY_SIZE']))
    dupRate.append(round(dupData['PERCENT_DUPLICATION'], 4))
    filtFrags.append(int(dupV[1] - dupV[4]))

    # get align stats
    statData = np.loadtxt(args.d + '/' + statFiles[j], dtype='str', delimiter='\t')
    idx = np.where(statData[:, 0] == 'chrM')[0][0]
    alnFragsNoM.append(np.sum(np.array(statData[1:idx, 2], int)) / 2)
    chrMFrags.append(int(statData[idx][2]) / 2)
    chrMrate.append(round(chrMFrags[-1] / float(alnFragsNoM[-1] + chrMFrags[-1]), 4))
    idx = np.where(statData[:, 0] == 'Chromosome')[0][1] + 1
    totalFiltRds.append(np.sum(np.array(statData[idx:-1][:, 2], int)) / 2)

    # get iSize
    iSizeHeaders = np.loadtxt(args.d + '/' + iSizeFiles[j], dtype='str', delimiter=',')[0].split('\t')
    iSizeVals = np.loadtxt(args.d + '/' + iSizeFiles[j], dtype='str', delimiter=',')[1].split('\t')
    iSizeData = dict(zip(iSizeHeaders, iSizeVals))
    mediSize.append(iSizeData['MEDIAN_INSERT_SIZE'])
    meaniSize.append(round(float(iSizeData['MEAN_INSERT_SIZE']), 2))
    try:
        stdiSize.append(round(float(iSizeData['STANDARD_DEVIATION']), 2))
    except:
        stdiSize.append(0)
    widthiSize.append(iSizeData['WIDTH_OF_80_PERCENT'])

    # raw counts
    iSizeData = np.loadtxt(args.d + '/' + iSizeFiles[j], dtype='str', delimiter=',')[3:]
    outIsize = np.zeros(1000)
    for l in range(0, len(iSizeData)):
        k = iSizeData[l]
        try:
            outIsize[int(k.split('\t')[0]) - 1] = int(k.split('\t')[1])
        except:
            continue
    iSizeMat.append(outIsize)
    # print(iSizeMat)

    # get enrichment stats
    # tssEnrich.append(round(np.sum(np.loadtxt(options.d+'/'+tssFiles[j]))/(totalFiltRds[-1]*2),4)) ## *2 bec TSS is in reads not frags
    mat = np.loadtxt(args.d + '/' + tssFiles[j])
    win = 20
    #print(round(max(np.convolve(mat, np.ones(win), 'same'))/win / np.mean(mat[1:200]), 2))
    #print(tssEnrich)
    tssEnrich.append(round(max(np.convolve(mat, np.ones(win), 'same'))/win / np.mean(mat[1:200]), 2))
    #tssEnrich.append(round(max(np.convolve(mat, np.ones(win), 'same'))/win / np.mean(mat[1:200]), 2)))
    # tssEnrich.append(round(max(np.convolve(mat, np.ones(win), 'same') / win / np.mean(mat[1:200]), 2)))
    
# append to table
data = np.vstack((frags, alnRate, libSize, dupRate, filtFrags, alnFragsNoM, chrMFrags, chrMrate, totalFiltRds, mediSize, meaniSize, stdiSize, widthiSize, tssEnrich)).T
final = np.column_stack((table[vector], data))
np.savetxt(outName, final, fmt='%s', delimiter='\t')

# save iSize matrix
sio.savemat(outName + '.iSize.txt', {'iSize': np.array(iSizeMat), 'annotation': np.array(final, dtype=object)}, oned_as='column', do_compression='True')

# sign out
print("Created by Sai Ma.")
print("Completed.")
