import sys
import os
current_path = os.path.abspath(os.path.dirname(__file__))
current_path = current_path.replace("plotSample", "")
sys.path.append("..")
from utils.LogExceptions import LogExceptions
from utils.GeneGtfSnp import GeneGtf
from utils.GeneGtfSnp import ClinvarRead
import multiprocessing
import time
import ctypes
import pysam
from loguru import logger

class plotSample:
    def __init__(self,
                 chrom,
                 start,
                 stop,
                 svLen,
                 svType,
                 minSignal,
                 minQuality,
                 genotypes,
                 flanks,
                 samples,
                 inputFiles,
                 outPath,
                 imgPrefix,
                 imgFormat,
                 dpi,
                 genome,
                 imgwith,
                 imgheight):

        self.samples = samples
        self.inputFiles = inputFiles
        self.legendFlag = True
        self.flanks = flanks
        self.outPath = outPath
        self.imgPrefix = imgPrefix
        self.imgFormat = imgFormat
        self.dpi = dpi
        self.imgwith = imgwith
        self.imgheight = imgheight
        self.supportGene = False
        self.supportClinvar = False
        self.minSignal = minSignal
        self.minQuality = minQuality
        self.genome = genome
        self.gtfName = "/home/lz/work_space/Database/VISOR_test/Chr1_two_new/compare2_500/extract/plotSV/Homo_sapiens.GRCh38.110.chr.sorted.gtf.gz"
        self.clinvarName = "/home/lz/work_space/Database/VISOR_test/Chr1_two_new/compare2_500/extract/plotSV/clinvar.vcf.gz" 
        print('self.clinvarName')
        self.run(chrom, start, stop, svLen, svType, genotypes)

    def extractCIGAR(self, chrom, startFlank, stopFlank):
        c_cigs, c_alignStarts, c_alignStops, c_nums, c_samples, c_map_quality = [], [], [], [], [], []
        sampleSize = len(self.samples)
        for i in  range(sampleSize):
            sample = self.samples[i]
            bamFile = self.inputFiles[i]
            c_samples.append(sample)

            #try:
            if self.bam_cram:
                samFile =  pysam.AlignmentFile(bamFile,"rb")
            else:
                samFile = pysam.AlignmentFile(bamFile, "rc")
            #except:

            allReads =samFile.fetch(contig=chrom, start=startFlank, stop=stopFlank)
            num = 1
            cigs, alignStarts, alignStops, nums, map_quality, query_names, query_seqs, query_quality = [],[],[],[],[], [], [], []
            out = open('a.%d.fa'%num,'w')
            for read in allReads:
                cigSeq = read.cigarstring
                min_quality = read.mapping_quality
                if min_quality < self.minQuality:# or read.is_secondary == True:
                    continue 
                #if read.flag >=256 :
                #    continue
                #print(read.query_name, read.flag, min_quality)
                #print("min_quality:", min_quality, min_quality+500)
                print("is_secondary:", read.is_secondary, read.is_reverse, read.query_name.split('|')[0])
                if min_quality >60:
                    min_quality = 60
                map_quality.append(str(min_quality+500))
                cigs.append(cigSeq)
                query_names.append(read.query_name.split('|')[0])
                query_seqs.append(read.query_sequence)
                query_quality.append(read.to_dict()['qual'])
                alignStarts.append(str(read.reference_start))
                alignStops.append(str(read.reference_end))
                nums.append(str(num))
                num += 1
                out.write(">"+read.query_name.split('|')[0]+"\n")
                out.write(read.query_sequence+"\n")
            out.close()

            c_cigs.append(':'.join(cigs))
            c_alignStarts.append(':'.join(alignStarts))
            c_alignStops.append(':'.join(alignStops))
            c_nums.append(':'.join(nums))
            c_map_quality.append(':'.join(map_quality))
            c_query_names = '|'.join(query_names)
            c_query_seqs = ','.join(query_seqs)
            c_query_quality = ' '.join(query_quality)
        return c_cigs, c_alignStarts, c_alignStops, c_nums, c_samples, c_map_quality, c_query_names, c_query_seqs, c_query_quality

    def dealOutPutFormat(self, chrom, start, stop, svLen, svType):
        if self.imgPrefix != None:
            imgPrefix = self.imgPrefix
        else:
            imgPrefix = '_'.join(map(str, [chrom, start, stop, svLen, svType]))

        if self.imgFormat in ["png", "pdf", "svg", "jpeg", "tiff"]: # png pdf svg eps jpeg tiff
            c_outPathImg = os.path.join(self.outPath, imgPrefix + "." + self.imgFormat)
        else:
            logger.error('The %s output format is currently not supported at the present time.'%self.imgFormat)

        return c_outPathImg

    def dealBamSV(self, chrom, start, stop, svLen, svType, genotype):
        t1 = time.time()
        lib = ctypes.CDLL('%s/golang/plotSample.so'%current_path)
        lib.DealWith.argtypes = [
                                ctypes.c_char_p,
                                ctypes.c_char_p,
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
        stopFlank = stop + flankMax

        c_outPathImg = self.dealOutPutFormat(chrom, start, stop, svLen, svType)

        c_start = ctypes.c_int(start)
        c_stop = ctypes.c_int(stop)
        c_svlen = ctypes.c_int(svLen)
        print("self.minSignal:",self.minSignal)
        c_minSignal = ctypes.c_int(self.minSignal)
        c_chrom = chrom
        svInfos = [
                   "CHROM: " + chrom,
                   "POS: " + str(start),
                   "END: " + str(stop),
                   "SVTYPE: " + svType,
                   "SVLEN: " + str(svLen)
                   ]
        c_svInfos = '\n'.join(svInfos)

        c_cigs, c_alignStarts, c_alignStops, c_nums, c_samples, c_map_quality, c_query_names, c_query_seqs, c_query_quality = self.extractCIGAR(chrom, startFlank, stopFlank)

        if len(c_alignStarts)>0:
            gtf_obj = GeneGtf(-1)
            snp_obj = ClinvarRead(-1)
            if self.supportGene:
                gtf_obj = GeneGtf(chrom, start, stop, self.gtfName)
            if self.supportClinvar:
                snp_obj = ClinvarRead(chrom, start, stop, self.clinvarName)
            print("gtf_obj.genePos_str:",gtf_obj.genePos_str)
            #genome = "/home/lz/work_space/Database/VISOR_test/Genome_old/Genome.fa"
            lib.DealWith(
                         ','.join(c_cigs).encode(),
                         c_query_names.encode(),
                         c_query_seqs.encode(),
                         c_query_quality.encode(),
                         c_chrom.encode(),
                         svType.encode(),
                         ','.join(map(str, self.flanks)).encode(),
                         self.genome.encode(),
                         c_start,
                         c_stop,
                         c_svlen,
                         c_minSignal,
                         ctypes.c_int(self.dpi),
                         ctypes.c_int(self.imgwith),
                         ctypes.c_int(self.imgheight),
                         c_outPathImg.encode(),
                         ','.join(map(str,c_alignStarts)).encode(),
                         ','.join(map(str,c_alignStops)).encode(),
                         ','.join(map(str,c_nums)).encode(),
                         ','.join(map(str,c_map_quality)).encode(),
                         ','.join(map(str,c_samples)).encode(),
                         ','.join(genotype).encode(),
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
        print("dealBedSV current_path:", self.minSignal, current_path)
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

        imgPrefix = '_'.join(map(str, [chrom, start, stop, svLen, svType]))
        if self.imgFormat in ["png", "pdf", "svg", "eps", "jpeg", "tiff"]: # png pdf svg eps jpeg tiff
            c_outPathImg = os.path.join(self.outPath, imgPrefix + "." + self.imgFormat)
        else:
            logger.error('The %s output format is not supported at present.'%self.imgFormat)

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

        for sample, bedGzFile in zip(self.samples, self.inputFiles):
            bedGz = pysam.TabixFile(bedGzFile)
            allSignals = bedGz.fetch(reference=c_chrom, start=startFlank, end=stop+flankMax)
            c_samples.append(sample)
            c_allSignals.append('|'.join(list(allSignals)))
            allSignals = bedGz.fetch(reference=c_chrom, start=startFlank, end=stop+flankMax)
            print('sample:', len([i for i in allSignals]), startFlank, stop, flankMax)
        #print(c_allSignals)
        print(c_samples)

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

    def run(self, chrom, start, stop, svLen, svType, genotypes):
        tmp_list = self.inputFiles[0].split(".")
        fileType = tmp_list[-1]
        fileBed  = '.'.join(tmp_list[-2:])
        print('fileType:',fileType, 'fileBed:',fileBed)
        if fileBed.lower() == "bed.gz":
            self.dealBedSV(chrom, start, stop, svLen, svType, genotypes)
        elif fileType.lower() == "bam":
            self.bam_cram = True
            print("bam===================")
            self.dealBamSV(chrom, start, stop, svLen, svType, genotypes)
        elif fileType.lower() == "cram":
            self.bam_cram = False
            self.dealBamSV(chrom, start, stop, svLen, svType, genotypes)
        else:
            logger.error('-i/--input, The input file type is incorrect. Please make sure that the input file is one of the three types of bam/cram/bed.')
