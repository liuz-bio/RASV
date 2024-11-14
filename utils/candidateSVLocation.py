import gzip
import json
import numpy as np
import pandas as pd
import copy
import time
import os
from sklearn.cluster import DBSCAN

#L:0 R:1
def appendRecord(chrom, query_start, query_end, is_reverse, line9, allDict, queryDict):
    for lx in line9:
        lx = lx.strip(",").split(',')
        if lx[2] =='1':
            if chrom in allDict:
                try:
                    allDict[chrom]['G'].append(int(lx[1]))
                    allDict[chrom]['Q'].append(int(query_start))
                    allDict[chrom]['RL'].append(0)#L
                    queryDict[query_start] = [query_end, is_reverse]
                except KeyError:
                    allDict[chrom]['G'] = [int(lx[1])]
                    allDict[chrom]['Q'] = [int(query_start)]
                    allDict[chrom]['RL'] = [0]#L
                    queryDict[query_start] = [query_end, is_reverse]
            else:
                allDict[chrom] = {'G':[int(lx[1])]}
                allDict[chrom]['Q'] = [int(query_start)]
                allDict[chrom]['RL']  = [0]#L
                queryDict[query_start] = [query_end, is_reverse]
        elif lx[2] == '4':
            if chrom in allDict:
                try:
                    allDict[chrom]['G'].append(int(lx[0]))
                    allDict[chrom]['Q'].append(int(query_end))
                    allDict[chrom]['RL'].append(1)#R
                    queryDict[query_end] = [query_start, is_reverse]
                except KeyError:
                    allDict[chrom]['G'] = [int(lx[0])]
                    allDict[chrom]['Q'] = [int(query_end)]
                    allDict[chrom]['RL'] = [1]#R
                    queryDict[query_end] = [query_start, is_reverse]
            else:
                allDict[chrom] = {'G':[int(lx[0])]}
                allDict[chrom]['Q'] = [int(query_end)]
                allDict[chrom]['RL'] = [1]#R
                queryDict[query_end] = [query_start, is_reverse]
    return allDict, queryDict

#allDict: {'chr1': {'genome': [2033372], 'query': [9046], 'RL': ['L']}, 'chr4': {'genome': [71592083], 'query': [58287], 'RL': ['L']}}
#allDict: {'chr1': {'genome': [2033372], 'query': [1167], 'RL': ['R']}, 'chr4': {'genome': [71592081], 'query': [3446], 'RL': ['R']}}
#DEL 0, INS 1, INV 2, DUP 3, BND 4
def dealBND(allDict, queryDict, records, query_length, reverseLists, outBnds):
    breakpointDeviation = 100
    reverseDeviation = 100
    chroms = list(allDict.keys())
    if len(chroms) < 1:
        return []

    for i, CHR1 in enumerate(chroms):
        for CHR2 in chroms[i:]:
            query_CHR1_size = len(allDict[CHR1]['Q'])
            query_CHR2_size = len(allDict[CHR2]['Q'])
            for i1 in range(query_CHR1_size):
                for i2 in range(i1+1, query_CHR2_size):
                    #print("CHR1==CHR2:",CHR1, CHR2, allDict[CHR1]['RL'][i1], allDict[CHR2]['RL'][i2], allDict[CHR1]['G'][i1], allDict[CHR2]['G'][i2], allDict[CHR1]['Q'][i1], allDict[CHR2]['Q'][i2])
                    CHR1_RL_i1 =  allDict[CHR1]['RL'][i1]
                    CHR2_RL_i2 = allDict[CHR2]['RL'][i2]
                    CHR1_G_i1 = allDict[CHR1]['G'][i1]
                    CHR2_G_i2 = allDict[CHR2]['G'][i2]
                    CHR1_Q_i1 = allDict[CHR1]['Q'][i1]
                    CHR2_Q_i2 = allDict[CHR2]['Q'][i2]
                    if CHR1 == CHR2:
                        if CHR1_RL_i1 == 0:#L
                            if CHR2_RL_i2 == 0:#L
                                if abs(CHR1_G_i1 - CHR2_G_i2) > breakpointDeviation:
                                    if abs(CHR1_Q_i1 + CHR2_Q_i2 - query_length) <= reverseDeviation:
                                        a = [CHR1_G_i1, CHR2_G_i2]
                                        a.sort()
                                        records[2].setdefault(CHR1, []).append([a[0], a[1], 2])
                            else:
                                if CHR1_G_i1 > CHR2_G_i2:
                                    if abs(CHR1_G_i1 - CHR2_G_i2) > breakpointDeviation:
                                        if abs(CHR1_Q_i1 + CHR2_Q_i2 - query_length) <= reverseDeviation:
                                            a = [CHR1_G_i1, CHR2_G_i2]
                                            a.sort()
                                            records[0].setdefault(CHR1, []).append([a[0], a[1], 0])
                                    elif abs(CHR1_Q_i1 - CHR2_Q_i2) > breakpointDeviation:
                                        a = [CHR1_G_i1, CHR2_G_i2]
                                        a.sort()
                                        records[1].setdefault(CHR1, []).append([a[0], a[1], 1])
                                else:
                                    #print('R>LDUP:', CHR1, CHR2, CHR1_G_i1, CHR2_G_i2, CHR1_Q_i1, CHR2_Q_i2, CHR1_RL_i1, CHR2_RL_i2, "query_length:", query_length)
                                    if abs(CHR1_Q_i1 - CHR2_Q_i2) <= breakpointDeviation:
                                        if abs(CHR1_G_i1 - CHR2_G_i2) > breakpointDeviation:
                                            a = [CHR1_G_i1, CHR2_G_i2]
                                            a.sort()
                                            records[3].setdefault(CHR1, []).append([a[0], a[1], 3])
                                    elif abs(queryDict[CHR1_Q_i1][0] - query_length) <= 100 and abs(queryDict[CHR2_Q_i2][0]) <= 100:
                                        a = [CHR1_G_i1, CHR2_G_i2]
                                        a.sort()
                                        records[3].setdefault(CHR1, []).append([a[0], a[1], 3])
                                    elif len(set(reverseLists)) == 2:
                                        tmp = queryDict[CHR1_Q_i1]
                                        reverseList = reverseLists.copy()
                                        reverseList.remove(tmp[1])
                                        if tmp[0] == CHR2_Q_i2 and len(set(reverseList) & {tmp[1]}) == 0:
                                            if abs(abs(CHR1_G_i1 - CHR2_G_i2) - abs(CHR1_Q_i1 - CHR2_Q_i2)) <= 50:
                                                a = [CHR1_G_i1, CHR2_G_i2]
                                                a.sort()
                                                records[2].setdefault(CHR1, []).append([a[0], a[1], 2])
                        elif CHR2_RL_i2 == 0:#L
                            if CHR1_G_i1 > CHR2_G_i2:
                                #print('R>LDUP:', CHR1, CHR2, CHR1_G_i1, CHR2_G_i2, CHR1_Q_i1, CHR2_Q_i2, CHR1_RL_i1, CHR2_RL_i2, "query_length:", query_length)
                                if abs(CHR1_Q_i1 - CHR2_Q_i2) <= breakpointDeviation:
                                    if abs(CHR1_G_i1 - CHR2_G_i2) > breakpointDeviation:
                                        a = [CHR1_G_i1, CHR2_G_i2]
                                        a.sort()
                                        records[3].setdefault(CHR1, []).append([a[0], a[1], 3])
                                elif abs(queryDict[CHR2_Q_i2][0] - query_length) <= 100 and abs(queryDict[CHR1_Q_i1][0]) <= 100:
                                    a = [CHR1_G_i1, CHR2_G_i2]
                                    a.sort()
                                    records[3].setdefault(CHR1, []).append([a[0], a[1], 3])
                                elif len(set(reverseLists)) == 2:
                                    tmp = queryDict[CHR1_Q_i1]
                                    reverseList = reverseLists.copy()
                                    reverseList.remove(tmp[1])
                                    if tmp[0] == CHR2_Q_i2 and len(set(reverseList) & {tmp[1]}) == 0:
                                        if abs(abs(CHR1_G_i1 - CHR2_G_i2) - abs(CHR1_Q_i1 - CHR2_Q_i2)) <= 100:
                                            a = [CHR1_G_i1, CHR2_G_i2]
                                            a.sort()
                                            records[2].setdefault(CHR1, []).append([a[0], a[1], 2])
                            else:
                                if abs(CHR1_G_i1 - CHR2_G_i2) > breakpointDeviation:
                                    if abs(CHR1_Q_i1 + CHR2_Q_i2 - query_length) <= reverseDeviation:
                                        a = [CHR1_G_i1, CHR2_G_i2]
                                        a.sort()
                                        records[0].setdefault(CHR1, []).append([a[0], a[1], 0])
                                elif abs(CHR1_Q_i1 - CHR2_Q_i2) > breakpointDeviation:
                                    a = [CHR1_G_i1, CHR2_G_i2]
                                    a.sort()
                                    records[1].setdefault(CHR1, []).append([a[0], a[1], 1])
                        else:
                            if abs(CHR1_G_i1 - CHR2_G_i2) > breakpointDeviation:
                                if abs(CHR1_Q_i1 + CHR2_Q_i2 - query_length) <= reverseDeviation:
                                    a = [CHR1_G_i1, CHR2_G_i2]
                                    a.sort()
                                    records[2].setdefault(CHR1, []).append([a[0], a[1], 2])
                    elif len(set([CHR1_RL_i1, CHR2_RL_i2])) == 2:
                        idx = '|'.join([CHR1, CHR2])
                        if CHR1_RL_i1 == 1 and CHR1_G_i1 <= CHR2_G_i2: #R
                            records[4].setdefault(idx, []).append([CHR1_G_i1, CHR2_G_i2, 4, CHR1_RL_i1, CHR2_RL_i2])
                            outBnds.write('\t'.join(map(str, [idx,CHR1_G_i1, CHR2_G_i2, 4,CHR1_Q_i1,CHR2_Q_i2,CHR1_RL_i1,CHR2_RL_i2]))+'\n')
                        elif CHR1_G_i1 >= CHR2_G_i2:
                            records[4].setdefault(idx, []).append([CHR1_G_i1, CHR2_G_i2, 4, CHR1_RL_i1, CHR2_RL_i2])
                            outBnds.write('\t'.join(map(str, [idx,CHR1_G_i1, CHR2_G_i2, 4,CHR1_Q_i1,CHR2_Q_i2,CHR1_RL_i1,CHR2_RL_i2]))+'\n')
                    elif abs(CHR1_Q_i1 + CHR2_Q_i2 - query_length) <= 200:
                        idx = '|'.join([CHR1, CHR2])
                        records[4].setdefault(idx, []).append([CHR1_G_i1, CHR2_G_i2, 4, CHR1_RL_i1, CHR2_RL_i2])
                        outBnds.write('\t'.join(map(str, [idx,CHR1_G_i1, CHR2_G_i2, 4,CHR1_Q_i1,CHR2_Q_i2,CHR1_RL_i1,CHR2_RL_i2]))+'\n')

def splitSigs(chrom, svType, record, dbscan):
    groups = []
    step = 700
    record_array = np.array(record, dtype=int)
    record_array = record_array[np.argsort(record_array[:, 0])]
    record_size = len(record_array)

    if record_size >1000:
        start_index = 0
        while True:
            head = record_array[start_index:step, :]
            if len(head) == 0:
                break
            else:
                stop_index = start_index+step
                if stop_index >= record_size:
                    break
                else:
                    while True:
                        if np.linalg.norm(head[-1] - record_array[stop_index,:]) >1500:
                            start_index = stop_index
                            break
                        else:
                            head = np.append(head, record_array[stop_index,:])
                            stop_index += 1

            if svType == "BND":
                head = head.reshape(-1, 5)
                clusterBND(chrom, svType, head, dbscan, groups)
            else:
                head = head.reshape(-1, 3)
                cluster(chrom, svType, head, dbscan, groups)

    elif svType == "BND":
        clusterBND(chrom, svType, record_array, dbscan, groups)
    else:
        cluster(chrom, svType, record_array, dbscan, groups)

    return groups

def cluster(chrom, svType, record, dbscan, groups):
    clusters = dbscan.fit_predict(record[:,:2])
    labels = np.unique(clusters)
    for label in labels:
        cluster_points = record[clusters == label]
        if label == -1:
            for i in cluster_points:
                groups.append([chrom, i[0], i[1], svType])
    else:
        pos1, pos2 = np.mean(cluster_points[:, :2], axis=0).astype(int)
        groups.append([chrom, pos1, pos2, svType])


def clusterBND(chrom, svType, record, dbscan, groups):
    clusters = dbscan.fit_predict(record[:,:2])
    labels = np.unique(clusters)
    groups = []
    for label in labels:
        cluster_points = record[clusters == label]
        if label == -1:
            for i in cluster_points:
                groups.append([chrom, i[0], i[1], svType, i[3], i[4]])
        else:
            pos1, pos2 = np.mean(cluster_points[:,:2], axis=0).astype(int) 
            rl1 = pd.Series(cluster_points[:,3]).value_counts()
            rl2 = pd.Series(cluster_points[:,4]).value_counts()
            groups.append([chrom, pos1, pos2, svType, rl1.idxmax(), rl2.idxmax()])

def candidateSVLocation(bed_gz_file, outPath, sample):
    startTime = time.time()
    #outDUPs = open('outDUPs.txt','w')
    outBnds = open(os.path.join(outPath, sample+".outBnds.txt"),'w')
    with gzip.open(bed_gz_file, 'rt') as f:
        while True:
            head = next(f).split('\t')
            if len(head) >5:
                break
        allDict = {}
        queryDict = {}
        reverseList = []
        #records = {'DEL':{},'INS':{},'INV':{},'DUP':{},'BND':{}}
        records = [{},{},{},{},{}]
        svTypeDict = {0:'DEL',1:'INS',2:'INV',3:'DUP',4:'BND'}
        dbscan = DBSCAN(eps=500, min_samples=2)
        for line in f:
            line = line.strip().split('\t')
            if head[9] != line[9]:
                if len(allDict) >0:
                    for chrom in allDict.keys():
                        sorted_index = sorted(range(len(allDict[chrom]['Q'])), key=lambda k: allDict[chrom]['Q'][k])
                        allDict[chrom]['G'] = [allDict[chrom]['G'][i] for i in sorted_index]
                        allDict[chrom]['Q'] = [allDict[chrom]['Q'][i] for i in sorted_index]
                        allDict[chrom]['RL'] = [allDict[chrom]['RL'][i] for i in sorted_index]
                    query_length = int(head[5])
                    #print("allDict:", allDict)
                    dealBND(allDict, queryDict, records, query_length, reverseList, outBnds)
                allDict = {}
                queryDict = {}
                reverseList = []
                head = line
                continue
            elif int(line[6]) >= 20 :
                if len(allDict) <1:
                    reverseList.append(int(head[7]))
                    reverseList.append(int(line[7]))
                    allDict, queryDict = appendRecord(head[0], int(head[3]), int(head[4]), int(head[7]), head[10:], allDict, queryDict)
                    allDict, queryDict = appendRecord(line[0], int(line[3]), int(line[4]), int(line[7]), line[10:], allDict, queryDict)
                else:
                    allDict, queryDict = appendRecord(line[0], int(line[3]), int(line[4]), int(line[7]), line[10:], allDict, queryDict)
                    reverseList.append(int(line[7]))

        print("Time:", time.time()-startTime)
        #out_raw = open('out_raw.tab','w')
        with open(os.path.join(outPath, sample+".out.tab"), 'w') as out:
            for i in range(len(records)):
                svType = svTypeDict[i]
                with  open(os.path.join(outPath, sample+"."+svType+".sigs"),'w') as outsvType:
                    for chrom in records[i].keys():
                        records[i][chrom] = splitSigs(chrom, svType, records[i][chrom], dbscan)
                        for ii in  records[i][chrom]:
                            out.write('\t'.join(map(str, ii))+'\n')
                            outsvType.write('\t'.join(map(str, ii))+'\n')
   
    outBnds.close()
