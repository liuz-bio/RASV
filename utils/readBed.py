from utils.svInfo import svInfo
from loguru import logger
class recoder:
    def __init__(self, line):
        try:
            self.chrom = line[0]
            self.start = int(line[1])
            self.end = int(line[2])
            self.svLen = int(line[3])
            self.svType = line[4]
            self.genotypes = line[5]
            self.filter = line[6]
            self.imgPrefix = '_'.join(line[:5])
        except ValueError:
            logger.error('Incorrect format of bed file. Note: the information in the first six columns of ' 
                         'bed file is "chrom, start, end, svlen, svtype, genotype, filter" respectively.')


class readBed:
    def __init__(self, bed):
        self.recoders = self.run(bed)

    def selet(self, filter, chrom, svtype):
        if len(self.filter) == 0:
            if len(self.chroms) == 0:
                if svtype in self.svtypes:
                    return False
                else:
                    return True
            elif chrom in self.chroms and svtype in self.svtypes:
                return False
            else:
                return True
        elif len(filter) == 0:
            filter.append("PASS")
            if len( set(self.filter).intersection(filter)) >0:
                if len(self.chroms) == 0:
                    if svtype in self.svtypes:
                        return False
                    else:
                        return True
                elif chrom in self.chroms and svtype in self.svtypes:
                    return False
                else:
                    return True
            else:
                return  True
            return False
    def run(self, bed):
        recoders = []
        for line in open(bed, 'r'):
            # chrom start end svlen svtype genotype
            line = line.strip().split('\t')
            recode = recoder(line)
            if self.selet(recode.filter, recode.chrom, recode.svType):
                recoders.append(recode)
        return recoders

