import vcf
import random

class recoder:
    def __init__(self, line, record):
        try:
            #print(line)
            self.chrom = line[0]
            self.start = int(line[1])
            self.end = int(line[2])
            self.svLen = int(line[3])
            self.svType = line[4]
            self.genotypes = line[5]
            self.record = record
        except ValueError:
            logger.error('Incorrect format of bed file. Note: the information in the first six columns of '
                         'bed file is "chrom, start, end, svlen, svtype, genotype" respectively.')

class readVcf:
    def __init__(self,vcffile, chroms=[], filter=["PASS"] , svtypes=["DEL", "INS", "INV", "DUP", "BND", "TRA"]):
        self.filter = filter
        self.chroms = chroms
        self.svtypes = svtypes
        self.recoders = self.run(vcffile, chroms, filter, svtypes)

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

    def run(self, vcffile, chroms, filter, svtypes):
        recoders = []
        vcf_reader = vcf.Reader(filename=vcffile)
        for record in vcf_reader:
            chrom = record.CHROM
            svtype = record.INFO['SVTYPE'].upper() #["DEL", "INS", "INV", "DUP", "BND", "TRA"]
            filter = record.FILTER
            if self.selet(filter, chrom, svtype):
                continue

            #recodN = recod
            start = record.POS
            gt = record.samples[0].data.GT
            try:
                svlen = abs(record.INFO['SVLEN'])
                if svlen <30:
                    continue
            except:
                continue
            try:
                end = abs(record.INFO['END'])
            except KeyError:
                if svtype in ['BND','TRA','INS']:
                    end = start + 1
                elif svtype in ["DEL","DUP",'INV']:
                    end = start + svlen

            recoders.append(recoder([chrom,start,end,svlen,svtype,gt], record))
        return recoders

