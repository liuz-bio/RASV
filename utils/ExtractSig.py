import sys
sys.path.append("..")
import multiprocessing
from utils.LogExceptions import LogExceptions
from utils.candidateSVLocation import candidateSVLocation
import pysam
import os
current_path = os.path.abspath(os.path.dirname(__file__))
current_path = current_path.replace("utils", "")
print("当前路径为：%s ExtractSignal" % current_path)
from loguru import logger
import ctypes
import json
import time
import subprocess
import numpy as np


class ExtractSig:
    def __init__(self, samples, bamFiles, mutProces, outPath, chromLengths, minSignal, bam_cram, winSize, is_call=True):
        self.sampleToBams = dict(zip(samples, bamFiles))
        self.mutProces = mutProces
        self.minSignal = minSignal
        self.winSize = winSize
        self.outPath = outPath
        self.bam_cram = bam_cram
        self.is_call = is_call
        self.chromLengthGroups = self.splitChrom(chromLengths)

    def splitChrom(self, chromLengths):
        chromLengthGroups = {}
        for chrom, length in chromLengths.items():
            length = int(length)
            chromLengthGroups[chrom] = []
            if length < self.winSize:
                chromLengthGroups[chrom].append([1, length])
            else:
                multipleTimes = length // self.winSize
                for i in range(multipleTimes):
                    chromLengthGroups[chrom].append([i * self.winSize + 1, (i + 1) * self.winSize])
                if length % self.winSize > 0:
                    chromLengthGroups[chrom].append([multipleTimes * self.winSize + 1, length])
        return chromLengthGroups

    def mergeSignals(self, allFiles):
        task_pool = {}
        for sample, chrData in allFiles.items():
            bedName = os.path.join(self.outPath, sample + '.bed')
            bedNames = '  '.join([ ii for i in chrData.values() for ii in i])
            task1 = "cat %s |env LC_ALL=C sort -S 200M -k1V -k2n -k3n --unique --compress-program=pzstd  --parallel=5|bgzip --threads 5  >%s.gz"%(bedNames ,bedName)
            task2 = "cat %s |env LC_ALL=C sort -S 200M -k10V  --unique --compress-program=pzstd  --parallel=5|bgzip --threads 5  >%s.readName.gz"%(bedNames ,bedName)
            p1 = subprocess.Popen(task1, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            if self.is_call:
                p2 = subprocess.Popen(task2, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)
            else:
                p2 = p1
            task_pool[bedName] = [p1, p2,  [ ii for i in chrData.values() for ii in i]]

        bedNames = list(task_pool.keys())
        while len(task_pool) != 0:
            for bedName in bedNames:
                try:
                    intask_Popen1 = task_pool[bedName][0]
                    intask_Popen2 = task_pool[bedName][1]
                except KeyError:
                    continue
                if intask_Popen1.poll() != None and intask_Popen2.poll() != None:
                    # print(time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime()) + ' '+str(intask_Popen.pid)+': '+' finish...')
                    logger.info(time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime()) + ' '+str(intask_Popen1.pid)+': '+' finish...')
                    pysam.tabix_index(bedName + ".gz", force=True, seq_col=0, start_col=1, end_col=2, preset="bed")
                    for i in task_pool[bedName][2]:
                        try:
                            os.remove(i)
                        except FileNotFoundError:
                            continue
                    del task_pool[bedName]
                    if self.is_call:
                        candidateSVLocation(bedName + ".readName.gz", self.outPath, sample)
                    #os.remove(bedName)

    def run(self):
        t1 = time.time()
        func = Func()
        multiprocessing.log_to_stderr()
        process_pool = multiprocessing.Pool(self.mutProces)
        allFiles = {}
        for sample, bamfile in self.sampleToBams.items():
            for chrom, lengths in self.chromLengthGroups.items():
                for length in lengths:
                    start, stop = length[0], length[1]
                    tmpBed = str(
                        os.path.join(self.outPath, '_'.join(map(str, [sample, chrom, start, stop])) + ".bed"))
                    if sample in allFiles:
                        if chrom in allFiles[sample]:
                            allFiles[sample][chrom].append(tmpBed)
                        else:
                            allFiles[sample][chrom] = [tmpBed]
                    else:
                        allFiles[sample] = {chrom: [tmpBed]}
                    #self.extractSignals(sample, bamfile, chrom, length[0], length[1], tmpBed)
                    process_pool.apply_async(func=LogExceptions(self.extractSignals),
                                             args=(sample, bamfile, chrom, length[0], length[1], tmpBed, ),
                                             callback=func.call_back, error_callback=func.err_call_back)
        process_pool.close()
        process_pool.join()
        t2 = time.time()
        print(t2 - t1)
        self.mergeSignals(allFiles)
        print(time.time() - t2)

    def dealReads(self, read):
        signals = self.lib.DealWith(read.cigarstring.encode(),
                               ctypes.c_int(read.reference_start),
                               read.query_sequence.encode(),
                               self.c_min_signal,
                              ).decode()
        #print(read.get_tag("NM"), "##################")
        is_reverse = '1' if read.is_reverse else '0'
        is_secondart = '1' if read.is_secondary else '0'
        return  [read.reference_name,
                read.reference_start,
                read.reference_end,
                read.query_alignment_start,
                read.query_alignment_end,
                read.query_length,
                read.mapping_quality,
                is_reverse,
                is_secondart,
                read.query_name.split('|')[0],
                signals]

    def low_quality(self, read):
        return read.mapping_quality >= 20

    def a(self, read):
        return 1

    def extractSignals(self, sample, bamfile, chrom, start, stop, tmpBed):
        lib = ctypes.CDLL('%s/golang/cigar.so'%current_path)
        lib.DealWith.argtypes = [ctypes.c_char_p, 
                                 ctypes.POINTER(ctypes.c_int),
                                 ctypes.c_char_p,
                                 ctypes.POINTER(ctypes.c_int),                                 
                                 ]

        lib.DealWith.restype = ctypes.c_char_p

        if self.bam_cram:
            samfile = pysam.AlignmentFile(bamfile, "rb")
        else:
            samfile = pysam.AlignmentFile(bamfile, "rc")

        allreads = samfile.fetch(contig=chrom, start=start, stop=stop)
        
        c_min_signal = ctypes.c_int(self.minSignal)

        #lines = map(self.dealReads, filter(self.low_quality, allreads))
        lines = []
        start_time = time.time()
        aa = 0
        for read in allreads:
            mapping_quality = read.mapping_quality
           
            if mapping_quality < 20:
                continue

            #a, b = read.get_cigar_stats()
            #ad = abs(a[1]- b[1])  < 30
            #ai = abs(a[2]- b[2])  < 30
            #ab = abs(a[4]-b[4]) < 30

            is_reverse = '1' if read.is_reverse else '0'
            is_secondary = '1'  if read.is_secondary else '0'
            query_name = read.query_name.split('|')[0]

            #if ab & ai & ab :
            #    signals = ""
            #else:
            signals = lib.DealWith(read.cigarstring.encode(),
                                   ctypes.c_int(read.reference_start),
                                   read.query_sequence.encode(),
                                   c_min_signal,
                                  ).decode().strip("\t")


            lines.append([
                chrom,
                read.reference_start,
                read.reference_end,
                read.query_alignment_start,
                read.query_alignment_end,
                read.query_length,
                mapping_quality,
                is_reverse,
                is_secondary,
                query_name,
                signals
            ])
            aa +=1


        start_time2 = time.time()
        if aa > 0:
            OutHandle = open(tmpBed, 'w')
            OutHandle.write('\n'.join(['\t'.join(map(str, i)) for i in lines])+"\n")
            OutHandle.close()
        #print("time:",time.time()-start_time, aa, time.time()-start_time2)

class Func(object):
    def __init__(self):
        # 利用匿名函数模拟一个不可序列化象
        # 更常见的错误写法是，在这里初始化一个数据库的长链接
        self.num = lambda: None

    def work(self, num=None):
        self.num = num
        return self.num

    @staticmethod
    def call_back(res):
        #print(f'Hello,World! {res}')
        pass

    @staticmethod
    def err_call_back(err):
        print(f'出错啦~ error：{str(err)}')

'''
if __name__ == "__main__":
    sampleToBams = {
        "Father": "/home/lz/work_space/Database/VISOR_test/Chr1_two_new/T2/FY0/VISOR_LASeR_Father_DEL.DUP.INS.INV/sim.srt.bam",
        "Mother": "/home/lz/work_space/Database/VISOR_test/Chr1_two_new/T2/FY0/VISOR_LASeR_Mother_DEL.DUP.INS.INV/sim.srt.bam",
        "Son": "/home/lz/work_space/Database/VISOR_test/Chr1_two_new/T2/FY0/VISOR_LASeR_Son_DEL.DUP.INS.INV/sim.srt.bam"}
    samples = ["Father", "Mother", "Son"]
    bamFiles = ["/home/lz/work_space/Database/VISOR_test/Chr1_two_new/T2/FY0/VISOR_LASeR_Father_DEL.DUP.INS.INV/sim.srt.bam",
                "/home/lz/work_space/Database/VISOR_test/Chr1_two_new/T2/FY0/VISOR_LASeR_Mother_DEL.DUP.INS.INV/sim.srt.bam",
                "/home/lz/work_space/Database/VISOR_test/Chr1_two_new/T2/FY0/VISOR_LASeR_Son_DEL.DUP.INS.INV/sim.srt.bam"]
    genome = "/home/lz/work_space/Database/VISOR_test/Genome_old/Genome.fa"
    chromLengths = dict([i.strip().split('\t')[:2] for i in open(genome + ".fai").read().strip().split('\n')])
    mutProces = 90
    minSignal = 30
    winSize = 10000000
    outPath = "./tmp"
    #samples, bamFiles, mutProces, outPath, chromLengths, minSignal, winSize=10000000
    ext = ExtractSignals(samples, bamFiles, mutProces, outPath, chromLengths,  minSignal, winSize)
    ext.run()
'''
