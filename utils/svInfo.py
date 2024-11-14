from loguru import logger
class svInfo:
    def __int__(self, chrom, start, stop, svLen, svType, genotypes): #chrom,start,end,svlen,svtype,gt
        try:
            self.chrom = chrom
            self.start = int(start)
            self.stop = int(stop)
            self.svLen = int(svLen)
            self.svType =svType.upper()
            self.imgPrefix = '_'.join(map(str, [chrom, start, stop, svLen, svType]))
        except ValueError:
            logger.error('Incorrect format of bed file. Note: the information in the first six columns of ' 
                         'bed file is "chrom, start, end, svlen, svtype, genotype" respectively.')


