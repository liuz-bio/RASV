import pysam
import ctypes

class ClinvarRead:
    def __init__(self, chrom, start=None, stop=None, fileName=None):
        if chrom == -1:
            self.snp_str = ''
        else:
            self.snp_str = self.getSnpPos(chrom, start, stop, fileName)
        #print(self.snp_str)
    def getSnpPos(self, chrom, start, stop, fileName):
        tbx = pysam.TabixFile(fileName, encoding='utf-8')
        snpPos = []
        print(chrom, start, stop)
        for row in tbx.fetch(chrom, start, stop):
            row = str(row).split('\t')
            snpPos.append(row[1])
        return ','.join(snpPos)

# geneIDstr :=  "gene1;gene2"
#tranIDstr :=  "ta1,ta2;tb1,tb2"
#strandstr :=  "+;-"
#genePos_str := "100,200;50,300"
#tranPos_str := "110,190|150,190;60,290|100,190"
#exonPos_str := "120,150:170,180|160,170:180,190;70,150:180,280|120,130:150,180"
class GeneGtf:
    def __init__(self, chrom, start=None, stop=None, gtfName=None):
        if chrom == -1:
            print('chrom:', chrom) 
            self.geneIDstr = ''
            self.genePos_str = ''
            self.strandstr = ''
            self.tranPos_str = ''
            self.fivePos_str = ''
            self.threePos_str = ''
            self.tranIDstr = ''
            self.exonPos_str = ''
        else:
            self.chrom = chrom
            self.start = start
            self.stop = stop
            self.tbx = pysam.TabixFile(gtfName)
            self.run()

    def deal_info(self, line):
        out = []
        for i in line.strip(";").split(";"):
            i = i.strip(" ").split(" ")
            out.append([i[0], ' '.join(i[1:])])
        return dict(out)

    def add_feature(self, feature, genesDict, pos, geneID, transcript_id):
        if transcript_id in genesDict[geneID]:
            if feature in genesDict[geneID][transcript_id]:
                genesDict[geneID][transcript_id][feature].append(pos)
            else:
                genesDict[geneID][transcript_id][feature] = [pos]
        else:
            genesDict[geneID][transcript_id] = {feature:[pos]}
        return genesDict

    def dictToStr(self, genesDict):
        tranIDs = []
        tranPos = []
        exonPos = []
        fivePos = []
        threePos = []
        #tranIDstr :=  "ta1,ta2;tb1,tb2"
        #tranPos_str := "110,190|150,190;60,290|100,190"
        #exonPos_str := "120,150:170,180|160,170:180,190;70,150:180,280|120,130:150,180"
        for geneID in self.geneIDs:
            transcriptDat =  genesDict[geneID]
            tranID = []
            tran_Pos = []
            exon_Pos = []
            five_utr_Pos = []
            three_utr_Pos = []
            for transcript, featureDat in transcriptDat.items():
                tranID.append(transcript)
                tran_Pos.append(featureDat['transcript_id'])
                tmp_exon_pos = []
                exon_Pos.append(':'.join(featureDat['exon']))
                try:
                    five_utr_Pos.append(':'.join(featureDat['five_prime_utr']))
                except KeyError:
                    five_utr_Pos.append("-1,-1")
                try:
                    three_utr_Pos.append(':'.join(featureDat['three_prime_utr']))
                except KeyError:
                    three_utr_Pos.append("-1,-1")

            tranIDs.append(','.join(tranID))
            tranPos.append('|'.join(tran_Pos))
            fivePos.append('|'.join(five_utr_Pos))
            threePos.append('|'.join(three_utr_Pos))
            exonPos.append('|'.join(exon_Pos))
        self.tranPos_str = ';'.join(tranPos)
        self.fivePos_str = ';'.join(fivePos)
        self.threePos_str = ';'.join(threePos)
        self.tranIDstr = ';'.join(tranIDs)
        self.exonPos_str = ';'.join(exonPos)
        #print("self.tranPos_str:",self.tranPos_str)
        #print("self.fivePos_str:",self.fivePos_str)
        #print("self.threePos_str:",self.threePos_str)
        #print("self.tranIDstr:",self.tranIDstr)
        #print("self.exonPos_str:",self.exonPos_str)

    def transcript(self, geneToPos):
        genesDict = {}
        c_transcript = []
        geneTotranscript = []
        transcriptToPos = []
        transcriptToexonPos = []
        for geneID, startStop in geneToPos.items():
            genesDict[geneID] = {}
            for row in self.tbx.fetch(self.chrom, startStop[0], startStop[1]):
                row = str(row).split('\t')
                infos = self.deal_info(row[8])
                feature = row[2].lower()
                tgeneID = infos["gene_id"].replace('"','')
                if geneID != tgeneID:
                    continue
                try:
                    transcript_id = infos["transcript_id"].replace('"','')
                except KeyError:
                    continue
                if feature == "transcript":
                    try:
                        transcript_name = infos["transcript_name"].replace('"','')
                    except:
                        transcript_name = ""

                    if transcript_id in genesDict[geneID]:
                        genesDict[geneID][transcript_id]['transcript_id'] = ','.join(row[3:5])
                    else:
                        genesDict[geneID][transcript_id] = {'transcript_id':','.join(row[3:5])}

                elif row[2] == "exon":
                    pos = ','.join(row[3:5])
                    genesDict = self.add_feature(feature, genesDict, pos, geneID, transcript_id)

                elif feature == "five_prime_utr":
                    pos = ','.join(row[3:5])
                    genesDict = self.add_feature(feature, genesDict, pos, geneID, transcript_id)

                elif feature == "three_prime_utr":
                    pos = ','.join(row[3:5])
                    genesDict = self.add_feature(feature, genesDict, pos, geneID, transcript_id)
        #print("genesDict:", genesDict)
        self.dictToStr(genesDict)

    def run(self):
        strands = []
        self.geneIDs = []
        genePos = []
        GeneToName = {}
        geneToPos = {}

        for row in self.tbx.fetch(self.chrom, self.start, self.stop):
            row = str(row).split('\t')
            #print(row)
            if row[2] == "gene":
                infos = self.deal_info(row[8])
                geneID = infos["gene_id"].replace('"','')
                try:
                    geneName = infos["gene_name"].replace('"','')
                except:
                    geneName = ""
                GeneToName[geneID] = geneName
                self.geneIDs.append(geneID)
                geneToPos[geneID] = [int(row[3]), int(row[4])]
                genePos.append(','.join(row[3:5]))
                strands.append(row[6])

        self.geneIDstr = ';'.join(self.geneIDs)
        self.genePos_str = ';'.join(genePos)
        self.strandstr = ';'.join(strands)
        #if len(self.geneIDs)<1
        #print("geneIDstr", self.geneIDstr)
        #print("genePos_str", self.genePos_str)
        #print("strandstr", self.strandstr)
        self.transcript(geneToPos)
'''
if __name__ == "__main__":
    gtfName = "Homo_sapiens.GRCh38.110.chr.sorted.gtf.gz"
    fileName = "clinvar.vcf.gz"
    chrom = "1"
    #start = 1203128
    start  = 1215053
    stop  = 1215153
    gtf_obj = GeneGtf(chrom, start, stop, gtfName)
    snp_obj = ClinvarRead(chrom, start, stop, fileName)

    lib = ctypes.CDLL('./gene.so')
    lib.DealWith.argtypes = [
                             ctypes.c_char_p,
                             ctypes.c_char_p,
                             ctypes.c_char_p,
                             ctypes.c_char_p,
                             ctypes.c_char_p,
                             ctypes.c_char_p,
                             ctypes.c_char_p,
                             ctypes.c_char_p,
                             ctypes.c_char_p
                            ]

    lib.DealWith(
                 gtf_obj.geneIDstr.encode(),
                 gtf_obj.tranIDstr.encode(),
                 gtf_obj.strandstr.encode(),
                 gtf_obj.genePos_str.encode(),
                 gtf_obj.tranPos_str.encode(),
                 gtf_obj.exonPos_str.encode(),
                 gtf_obj.fivePos_str.encode(),
                 gtf_obj.threePos_str.encode(),
                 snp_obj.snp_str.encode()
                )
'''
