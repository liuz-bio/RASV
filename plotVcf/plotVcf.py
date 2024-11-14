import os
import sys
current_path = os.path.abspath(os.path.dirname(__file__))
current_path = current_path.replace("plotVcf", "")
import pysam
import ctypes
from utils.LogExceptions import LogExceptions
from utils.readVcf import readVcf
from utils.GeneGtfSnp import GeneGtf
from utils.GeneGtfSnp import ClinvarRead
from utils.ExtractSignal import ExtractSignals
import multiprocessing
from multiprocessing import Pool
import time
from loguru import logger


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


class plotVcf:
    def __init__(self, recoders, flanks, samples, inputFiles, outPath, imgFormat, dpi, imgwith, imgheight, mutProces, minSignal, minQuality, winSize):
        self.samples = samples
        self.recoders = recoders
        self.flanks = flanks
        self.legendFlag = True
        self.outPath = outPath
        self.minSignal = minSignal
        self.winSize = winSize
        self.minQuality = minQuality
        self.mutProces = mutProces
        self.imgFormat = imgFormat
        self.supportGene = False
        self.supportClinvar = False
        self.dpi = dpi
        self.imgwith = imgwith
        self.imgheight = imgheight
        self.gtfName = "/home/lz/work_space/Database/VISOR_test/Chr1_two_new/compare2_500/extract/plotSV/Homo_sapiens.GRCh38.110.chr.sorted.gtf.gz"
        self.clinvarName = "/home/lz/work_space/Database/VISOR_test/Chr1_two_new/compare2_500/extract/plotSV/clinvar.vcf.gz"
        self.creatBedGz(inputFiles)

    def creatBedGz(self, inputFiles):
        bedGzFormate = []
        bamCramFormate = []
        for i in inputFiles:
            i = i.split('.')
            bamCramFormate.append(i[-1].lower())
            bedGzFormate.append('.'.join(i[-2:]).lower())
        bamCram = list(set(bamCramFormate))
        print("bedGzFormate:",bedGzFormate, inputFiles)
        print("bamCram:", bamCram)
        os.makedirs(self.outPath, exist_ok=True)
        if list(set(bedGzFormate)) == ["bed.gz"]:
            print("bed.gz")
            self.bedGzFiles = inputFiles
            self.run()            
        elif bamCram == ["bam"]:
            print("bam")
            self.bam_cram = True
            self.signals(inputFiles)  
            print([os.path.join(self.outPath, ii+".bed.gz") for ii in self.samples])
            self.bedGzFiles = [os.path.join(self.outPath, ii+".bed.gz") for ii in self.samples]
            self.run()
        elif bamCram == ["cram"]:
            print("cram")
            sys.exit()
            self.bam_cram = False
            self.signals(inputFiles)
            self.bedGzFiles = [os.path.join(self.outPath, ii+".bed.gz") for ii in self.samples]
            self.run()
        else:
           logger.error('-i/--input, The input file type is incorrect. Please make sure that the input file is one of the three types of bam/cram/bed.') 

    def signals(self, bamFiles):
        samFile =  pysam.AlignmentFile(bamFiles[0], 'rb')
        chromLengths = {}
        for i in range(1000000):
            try:
                chrom = samFile.header.get_reference_name(i)
                length = samFile.header.get_reference_length(chrom)
                chromLengths[chrom] = length
            except ValueError:
                break

        os.makedirs(self.outPath, exist_ok=True)
        print("#samples, bamFiles, mutProces, outPath, chromLengths, minSignal, winSize=10000000")
        ext = ExtractSignals(
                             self.samples,
                             bamFiles,
                             self.mutProces,
                             self.outPath,
                             chromLengths,
                             self.minSignal,
                             self.bam_cram,
                             self.winSize,
                            )
        ext.run()

    def caculationfStartEnd(self,start, end):
        fStart, fStop = [], []

        for flank in self.flanks:
            fStart_tmp = start - flank
            fStop.append(str(end + flank))
            if fStart_tmp <= 0:
                fStart_tmp = 1
            fStart.append(str(fStart_tmp))

        return fStart, fStop

    def dealBedSV(self, chrom, start, stop, svLen, svType, genotypes):
        t1 = time.time()
        lib = ctypes.CDLL(os.path.join(current_path, "golang", "plotVcf.so"))
        print("current_path:",current_path)
        lib.DealWith.argtypes = [
                                ctypes.c_char_p,
                                ctypes.c_char_p,
                                ctypes.c_char_p,
                                ctypes.c_char_p,
                                ctypes.c_char_p,
                                ctypes.c_char_p,
                                ctypes.POINTER(ctypes.c_int),
                                ctypes.POINTER(ctypes.c_int),
                                ctypes.POINTER(ctypes.c_int),
                                ctypes.POINTER(ctypes.c_int),
                                ctypes.POINTER(ctypes.c_int),
                                ctypes.POINTER(ctypes.c_int),
                                ctypes.POINTER(ctypes.c_int),
                                ctypes.c_char_p,
                                ctypes.c_char_p,
                                ctypes.c_char_p,
                                ctypes.c_char_p,
                                ctypes.c_char_p,
                                ctypes.c_char_p,
                                ctypes.c_char_p,
                                ctypes.c_char_p,
                                ctypes.c_char_p,
                                ctypes.c_char_p,
                                ctypes.c_char_p,
                                ctypes.c_char_p,
                                ctypes.c_char_p,
                                ctypes.c_bool
                                ]
        flankMax = max(self.flanks)
        startFlank = start - flankMax
        if startFlank <=0:
            startFlank = 1

        #imgPrefix = '_'.join(map(str, [chrom, start, stop, svLen, svType]))
        #if self.imgFormat in ["png", "pdf", "svg", "eps", "jpeg", "tiff"]: # png pdf svg eps jpeg tiff
        #    c_outPathImg = os.path.join(self.outPath, imgPrefix + "." + self.imgFormat)
        #else:
        #    logger.error('The %s output format is not supported at present.'%self.imgFormat)

        c_start = ctypes.c_int(start)
        c_stop = ctypes.c_int(stop)
        c_dpi  = ctypes.c_int(self.dpi)
        c_minQuality = ctypes.c_int(self.minQuality)
        c_minSignal = ctypes.c_int(self.minSignal)
        c_chrom = chrom
        c_fStart, c_fStop = self.caculationfStartEnd(start, stop)
        svInfos = [
            "CHROM: " + chrom,
            "POS: " + str(start),
            "END: " + str(stop),
            "SVTYPE: " + svType,
            "SVLEN: " + str(svLen)
        ]
        c_svInfos = '\n'.join(svInfos)

        c_allSignals, c_samples = [],[]

        for sample, bedGzFile in zip(self.samples, self.bedGzFiles):
            imgPrefix = '_'.join(map(str, [chrom, start, stop, svLen, svType, sample]))
            if self.imgFormat in ["png", "pdf", "svg", "eps", "jpeg", "tiff"]: # png pdf svg eps jpeg tiff
                c_outPathImg = os.path.join(self.outPath, imgPrefix + "." + self.imgFormat)
            else:
                logger.error('The %s output format is not supported at present.'%self.imgFormat)

            bedGz = pysam.TabixFile(bedGzFile)
            allSignals = bedGz.fetch(reference=c_chrom, start=startFlank, end=stop+flankMax)
            c_samples.append(sample)
            c_allSignals.append('|'.join(list(allSignals)))

        if len(c_allSignals)>0:
            gtf_obj = GeneGtf(-1)
            snp_obj = ClinvarRead(-1)
            if self.supportGene:
                gtf_obj = GeneGtf(chrom, start, stop, self.gtfName)
            if self.supportClinvar:
                snp_obj = ClinvarRead(chrom, start, stop, self.clinvarName)

            lib.DealWith(
                        ':'.join(c_allSignals).encode(),
                         c_chrom.encode(),
                         svType.encode(),
                         ','.join(map(str, self.flanks)).encode(),
                         ','.join(c_fStart).encode(),
                         ','.join(c_fStop).encode(),
                         c_dpi,
                         ctypes.c_int(self.imgwith),
                         ctypes.c_int(self.imgheight),
                         c_minQuality,
                         c_minSignal,
                         c_start,
                         c_stop,
                         c_outPathImg.encode(),
                         ','.join(map(str,c_samples)).encode(),
                         ','.join(genotypes).encode(),
                         c_svInfos.encode(),
                         gtf_obj.geneIDstr.encode(),
                         gtf_obj.tranIDstr.encode(),
                         gtf_obj.strandstr.encode(),
                         gtf_obj.genePos_str.encode(),
                         gtf_obj.tranPos_str.encode(),
                         gtf_obj.exonPos_str.encode(),
                         gtf_obj.fivePos_str.encode(),
                         gtf_obj.threePos_str.encode(),
                         snp_obj.snp_str.encode(),
                         self.legendFlag,
                         )

            print("time:", time.time() - t1)
        else:
            logger.error('In the provided region, none of the samples existed alignment reads with the reference genome.')


    def run_dealSV(self, data):
        for lx in data:
            #self.dealSV(lx[0], lx[1], lx[2], lx[3], lx[4], lx[5])
            #chrom, start, stop, svLen, svType, genotypes
            print("lx:",lx)
            self.dealBedSV(lx[0], lx[1], lx[2], lx[3], lx[4], lx[5])

    def splitRecoders(self):
        nuSize = 10
        nu = 0
        datas = []
        tmp = []
        if len(self.recoders) <= nuSize:
            self.mutProces = 1
        for recoder in self.recoders:
            if nu < nuSize:
                tmp.append([recoder.chrom,recoder.start,recoder.end,recoder.svLen,recoder.svType,recoder.genotypes])
                nu+=1
            else:
                datas.append(tmp)
                tmp = [[recoder.chrom,recoder.start,recoder.end,recoder.svLen,recoder.svType,recoder.genotypes]]
                nu = 0
                continue
        if nu != 0:
            datas.append(tmp)
 
        #datas = [tmp]
        #sys.exit()
        return datas

           
    def run(self):
        print("self.mutProces:",self.mutProces)
        datas = self.splitRecoders()
        func = Func()
        multiprocessing.log_to_stderr()
        process_pool = multiprocessing.Pool(processes=self.mutProces)
        #process_pool = Pool(self.mutProces)

        for tmp in datas:
            process_pool.apply_async(self.run_dealSV, args=(tmp,),callback=func.call_back, error_callback=func.err_call_back)

        #for recoder in self.recoders:
        #    process_pool.apply_async(LogExceptions(self.dealSV), args=(recoder.chrom,
        #                                                               recoder.start,
        #                                                               recoder.end,
        #                                                               recoder.svLen,
        #                                                               recoder.svType,
        #                                                               recoder.genotypes,))
        process_pool.close()
        process_pool.join()
'''
vcfFile = "../../../compare/compare_TP_comm.vcf"
bedgz = "test/Father.bed.gz"
flanks = [500, 5000, 25000, 125000, 625000]
imgFormat = "png"
dpi = 100
outPath = "./test/"
legendFlag = True
mutProces = 1
samples = None
bamFiles = None
recoders =  readVcf(vcfFile).recoders
bedgzs = "bedGzs.bed"

plotVcf(recoders, flanks, samples , bamFiles, bedgzs, outPath, imgFormat, dpi, mutProces, legendFlag)
'''
