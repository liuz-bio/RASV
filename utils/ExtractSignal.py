import sys
sys.path.append("..")
import multiprocessing
from utils.LogExceptions import LogExceptions
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

class ExtractSignals:
    def __init__(self, samples, bamFiles, mutProces, outPath, chromLengths, minSignal, bam_cram, winSize=10000000):
        self.sampleToBams = dict(zip(samples, bamFiles))
        self.mutProces = mutProces
        self.minSignal = minSignal
        self.winSize = winSize
        self.outPath = outPath
        self.bam_cram = bam_cram
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
            out = open(bedName, 'w')
            # out =  open(bedName, 'wt', compresslevel=9)
            bedLines = []
            for chrom, tmpBeds in chrData.items():
                chromList = []  # sorted_data = sorted(data_list, key=lambda x: (x[0], x[1]))
                for tmpBed in tmpBeds:
                    if os.path.exists(tmpBed):
                        # with gzip.open(tmpBed, 'rt') as inputs:
                        with open(tmpBed, 'rt') as inputs:
                            for lx in inputs:
                                chromList.append(chrom + '\t' + lx)
                                # chromList.append([int(lx.split("\t")[0]), chrom+'\t'+lx])
                        os.remove(tmpBed)
                bedLines.append(''.join(list(set(chromList))))
            out.write(''.join(bedLines))
            out.close()
            #task = "LC_ALL=C sort -S 5%% -k1V -k2n -k3n %s |bgzip >%s.gz" % (bedName, bedName)
            task = "LC_ALL=C sort -S 5%% -k1V -k2n -k3n --parallel = 20 --compress-program=gzip %s |bgzip >%s.gz" % (bedName, bedName)
            p = subprocess.Popen(task, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            task_pool[bedName] = p

        bedNames = list(task_pool.keys())
        while len(task_pool) != 0:
            for bedName in bedNames:
                try:
                    intask_Popen = task_pool[bedName]
                except KeyError:
                    continue
                if intask_Popen.poll() != None:
                    # print(time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime()) + ' '+str(intask_Popen.pid)+': '+' finish...')
                    logger.info(time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime()) + ' '+str(intask_Popen.pid)+': '+' finish...')
                    pysam.tabix_index(bedName + ".gz", force=True, seq_col=0, start_col=1, end_col=2, preset="bed")
                    del task_pool[bedName]
                    os.remove(bedName)

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

    def extractSignals(self, sample, bamfile, chrom, start, stop, tmpBed):
        lib = ctypes.CDLL('%s/golang/extractSignal.so'%current_path)
        lib.DealWith.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p,
                                 ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p,
                                 ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p,
                                 ctypes.POINTER(ctypes.c_int)]
        lib.DealWith.restype = ctypes.c_char_p

        if self.bam_cram:
            samfile = pysam.AlignmentFile(bamfile, "rb")
        else:
            samfile = pysam.AlignmentFile(bamfile, "rc")

        allreads = samfile.fetch(contig=chrom, start=start, stop=stop)
        
        c_min_signal = ctypes.c_int(self.minSignal)
        cigs, alignStarts, alignStops = [], [], []
        query_names, map_quality, query_seqs = [], [], []
        query_is_reverse, query_is_secondary = [], []
        queryAlignStarts, queryAlignStops =[], []
        queryLengths = []

        for read in allreads:
            #if read.flag > 2000:
            #    continue
            map_quality.append(read.mapping_quality+500)
            query_names.append(read.query_name.split('|')[0])
            cigs.append(read.cigarstring)
            alignStarts.append(read.reference_start)
            alignStops.append(read.reference_end)
            queryAlignStarts.append(read.query_alignment_start)
            queryAlignStops.append(read.query_alignment_end)
            queryLengths.append(read.query_length)
            query_is_reverse.append('1' if read.is_reverse else '0')
            query_is_secondary.append('1' if read.is_secondary else '0')
            query_seqs.append(read.query_sequence)

        if len(alignStarts) > 0:
            lib.DealWith(
                         ','.join(cigs).encode(), 
                         tmpBed.encode(),
                         '|'.join(query_names).encode(), 
                         ','.join(map(str, alignStarts)).encode(),
                         ','.join(map(str, alignStops)).encode(), 
                         ','.join(map(str, queryAlignStarts)).encode(),
                         ','.join(map(str, queryAlignStops)).encode(),
                         ','.join(map(str, queryLengths)).encode(),
                         ','.join(map(str, map_quality)).encode(),
                         ','.join(query_is_reverse).encode(),
                         ','.join(query_is_secondary).encode(),
                         ','.join(query_seqs).encode(),
                         c_min_signal,
                         ).decode().split(',')

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
