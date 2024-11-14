import os
import sys
current_path = os.path.abspath(os.path.dirname(__file__))
print(current_path,'a')
current_path = current_path.replace("plotSample", "")
print(current_path,'b')
sys.path.append("..")
from plotSample.plotSample import plotSample
from plotVcf.plotVcf import plotVcf
#from callSV_new.caller import caller
from utils.readVcf import readVcf
from utils.readBed import readBed
from utils.ExtractSignal import ExtractSignals
from utils.ExtractSig import ExtractSig
import time
import pysam
import logging
from loguru import logger


def creatLogger(toolName,logPath):
    # 第一步，创建一个logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)  # Log等级总开关

    # 第二步，创建一个handler，用于写入日志文件
    NowTime = time.strftime('%Y%m%d%H%M', time.localtime(time.time()))[:-4]
    logFile =  os.path.join(logPath, NowTime + '.log')
    fh = logging.FileHandler(logFile, mode='a')
    fh.setLevel(logging.DEBUG)  # 输出到file的log等级的开关

    # 第三步，定义handler的输出格式
    formatter = logging.Formatter("%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s")
    fh.setFormatter(formatter)

    # 第四步，将logger添加到handler里面
    logger.addHandler(fh)

    # 日志
    return logger
    #logger.debug('this is a logger debug message')
    #logger.info('this is a logger info message')
    #logger.warning('this is a logger warning message')
    #logger.error('this is a logger error message')
    #logger.critical('this is a logger critical message')

def sample(parser, options, extra_args=None):
    #logger.add('runlog_{time}.sample.log',encoding='utf-8')

    #chrom, start, stop, svLen, svType, genotype, flanks, samples , bamFiles, outPath, imgFormat, legendFlag):
    os.makedirs(options.outPath, exist_ok=True)
    print("chrom, start, stop, svLen, svType, genotype, flanks, samples , bamFiles, outPath, imgFormat, legendFlag)")
    print("options.imgwith:", options.imgwith)
    print("options.imgheight:", options.imgheight)
    plotSample(
        options.chrom,
        options.start,
        options.stop,
        options.svLen,
        options.svType,
        options.minSignal,
        options.minQuality,
        options.genotype,
        options.flanks,
        options.samples,
        options.inputFiles,
        options.outPath,
        options.imgPrefix,
        options.imgFormat,
        options.dpi,
        options.genome,
        options.imgwith,
        options.imgheight,
    )

def args_sample(parent_parser):
    """Defines the allowed arguments for plot function
    """
    parser = parent_parser.add_parser(
        "sample",
        help="Plot an image of structural variant from "
             + "bam alignments, "
             + "optimized for structural variant.",
    )

    parser.add_argument(
        "-n",
        "--samples",
        help="Space-delimited list of sample name. ",
        type=str,
        nargs="+",
        required=False,
    )

    parser.add_argument(
        "-b",
        "--inputFiles",
        type=str,
        metavar='FILE',
        nargs="+",
        help="Space-delimited list of BAM/CRAM/signles_bed file names",
        required=True,
    )


    parser.add_argument(
        "-g",
        "--genome",
        help="Genome file for bam",
        type=str,
        required=False,
    )

    parser.add_argument(
        "-c",
        "--chrom",
        type=str,
        help="Chromosome name of structural variant in the genome.",
        required=True,
    )

    parser.add_argument(
        "-s",
        "--start",
        type=int,
        help="Start position of structural variant in genome.",
        required=True,
    )

    parser.add_argument(
        "-e",
        "--stop",
        type=int,
        help="Stop position of structural variant in genome.",
        required=True,
    )

    parser.add_argument(
        "-t",
        "--svType",
        type=str,
        help="Structural variant type(DEL,INS,INV,DUP,TRA,BND).",
        required=True,
    )

    parser.add_argument(
        "--minSignal",
        type=int,
        default=25,
        help="Minimum SV signal size. default:25",
        required=False,
    )

    parser.add_argument(
        "--minQuality",
        type=int,
        default=45,
        help="The lowest quality value of reads alignment.default:45",
        required=False,
    )

    parser.add_argument(
        "--genotype",
        type=str,
        default=None,
        nargs="+",
        help="Space-delimited list of Genotypes of structural variation",
        required=False,
    )
    parser.add_argument(
        "-l",
        "--svLen",
        type=int,
        help="Structural variant length.",
        required=True,
    )

    parser.add_argument(
        "-f",
        "--flanks",
        type=int,
        nargs="+",
        default=[500,5000,25000,125000,625000], #375000
        help="The gradient of the side view of the SV.",
        required=False,
    )

    parser.add_argument(
        "--imgFormat",
        type=str,
        default='png',
        help="The exported image format. Default: png. "+
             "One of png, jpeg, tiff, svg, pdf.",
        required=False,
    )
  
    parser.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="Generate image's DPI. ",
        required=False,
    )

    parser.add_argument(
        "--imgwith",
        type=int,
        default=4,
        help="Width of a single image. ",
        required=False,
    )

    parser.add_argument(
        "--imgheight",
        type=int,
        default=2,
        help="Height of a single image. ",
        required=False,
    )


    parser.add_argument(
        "--imgPrefix",
        type=str,
        default=None,
        help="Image file name; default:{chrom}_{start}_{stop}_{svLen}_{svType}.png.",
        required=False,
    )

    parser.add_argument(
        "-o",
        "--outPath",
        type=str,
        default='.',
        help="The path to the folder where the output image is stored.",
        required=False,
    )

    parser.set_defaults(func=sample)

def vcf(parser, options, extra_args=None):
    #logger.add('runlog_{time}.sample.log',encoding='utf-8')
    if options.vcfFile == None:
        if options.bedFile == None:
            logger.error("Please provide Vcf or bed file. Set parameter: '--vcfFile' or '--bedFile'")
            sys.exit()
        else:
            recoder = readBed(options.vcfFile, options.chroms, options.filter, options.svtypes)
    else:
        if options.bedFile == None:
            # vcffile, chroms, filter=["PASS"] , svtypes=["DEL", "INS", "INV", "DUP", "BND", "TRA"]):
            recoder = readVcf(options.vcfFile, options.chroms, options.filter, options.svtypes)
        else:
            logger.error("Parameters '--vcfFile' and '--bedFile' cannot be set at the same time")

    #baseData, flanks, samples , bamFiles, outPath, imgFormat, legendFlag
    #recoders, flanks, samples , bamFiles, outPath, imgFormat, dpi, mutProces
    #recoders, flanks, samples, bamFiles, bedrGzs, outPath, imgFormat, dpi, mutProces, legendFlag
    #recoder.recoders

    if len(options.samples) != len(options.inputFiles):
        logger.error("The sample name is inconsistent with the inputFiles number.")
 
    plotVcf(
            recoder.recoders, 
            options.flanks, 
            options.samples, 
            options.inputFiles, 
            options.outPath, 
            options.imgFormat, 
            options.dpi, 
            options.imgwith,
            options.imgheight,
            options.mutProces, 
            options.minSignal,
            options.minQuality,
            options.winSize,
          )

def args_vcf(parent_parser):
    """Defines the allowed arguments for plot function
    """
    """
       Defines the allowed arguments for plot function from vcf file.
       """
    parser = parent_parser.add_parser(
        "vcf",
        help="Plot multiple images of structural variant from vcf file accord to "
             + "bam alignments, "
             + "optimized for structural variant.",
    )

    parser.add_argument(
        "-f",
        "--vcfFile",
        help="The vcf file name of structural variant. ",
        type=str,
        required=False,
    )

    parser.add_argument(
        "--bedFile",
        help="The bed file name of structural variant. "+
             "bed file head is chrom, start, end, svlen, svtype, genotype.",
        type=str,
        required=False,
    )

    #parser.add_argument(
    #    "-i",
    #    "--input",
    #    type=str,
    #    default=None,
    #    help="A text file that stores two columns of information "
    #         + "sample name and signal bed file path",
    #    required=False,
    #)

    parser.add_argument(
        "-n",
        "--samples",
        help="Space-delimited list of sample name. ",
        type=str,
        nargs="+",
        required=True,
    )

    parser.add_argument(
        "-b",
        "--inputFiles",
        type=str,
        metavar='FILE',
        nargs="+",
        help="Space-delimited list of BAM/CRAM/signles_bed file names",
        required=True,
    )

    parser.add_argument(
        "--chroms",
        type=str,
        nargs="+",
        default=[],
        help="Chromosome name of structural variant in the genome.",
        required=False,
    )

    parser.add_argument(
        "--svtypes",
        type=str,
        nargs="+",
        default=["DEL", "INS", "INV", "DUP", "TRA", "BND"],
        help="Chromosome name of structural variant in the genome.",
        required=False,
    )

    parser.add_argument(
        "-g",
        "--genome",
        help="Genome file for bam",
        type=str,
        required=False,
    )

    parser.add_argument(
        "--mutProces",
        type=int,
        default=10,
        help="The number of processes involved in counting "
             + "the number of Reads alignment in bam files.",
        required=False,
    )

    parser.add_argument(
        "-r",
        "--flanks",
        type=int,
        nargs="+",
        default=[500, 5000, 25000, 125000, 375000],
        help="The gradient of the side view of the SV.",
        required=False,
    )

    parser.add_argument(
        "--minSignal",
        type=int,
        default=25,
        help="Minimum SV signal size. default:25",
        required=False,
    )

    parser.add_argument(
        "--minQuality",
        type=int,
        default=45,
        help="The lowest quality value of reads alignment.default:45",
        required=False,
    )

    parser.add_argument(
        "-w",
        "--winSize",
        type=int,
        default=10000000,
        help="The size of the chromosome interval processed "
             + "by each process. default:10000000",
        required=False,
    )

    parser.add_argument(
        "--filter",
        type=str,
        default=[],
        nargs="+",
        help="Select the SV of the specific field of the FILTER field in "
             +"the Vcf file for visualization. If not selected, visualize "
             +"all FILTER field in the Vcf file.",
        required=False,
    )

    parser.add_argument(
        "--imgFormat",
        type=str,
        default='png',
        help="The exported image format. Default: png. " +
             "One of png, jpeg, tiff, svg, pdf.",
        required=False,
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="Generate image's DPI. ",
        required=False,
    )

    parser.add_argument(
        "--imgwith",
        type=int,
        default=4,
        help="Width of a single image. ",
        required=False,
    )

    parser.add_argument(
        "--imgheight",
        type=int,
        default=2,
        help="Height of a single image. ",
        required=False,
    )

    parser.add_argument(
        "-o",
        "--outPath",
        type=str,
        default='.',
        help="The path to the folder where the output image is stored.",
        required=False,
    )
    parser.set_defaults(func=vcf)



def signals(parser, options, extra_args=None):
    if ".bam" in options.bamFiles[0]:
        bam_cram = True
    else:
        bam_cram = False
    samFile =  pysam.AlignmentFile(options.bamFiles[0], 'rb')
    chromLengths = {}
    for i in range(1000000):
        try:
            chrom = samFile.header.get_reference_name(i)
            length = samFile.header.get_reference_length(chrom)
            chromLengths[chrom] = length
        except ValueError:
            break
    os.makedirs(options.outPath, exist_ok=True)
    #samples, bamFiles, mutProces, outPath, chromLengths, minSignal, winSize=10000000
    # samples, bamFiles, mutProces, outPath, chromLengths, minSignal, bam_cram, winSize

    ext = ExtractSig(
                         options.samples, 
                         options.bamFiles, 
                         options.mutProces, 
                         options.outPath, 
                         chromLengths, 
                         options.minSignal, 
                         bam_cram,
                         options.winSize,
                        )
    ext.run()


def args_signals(parent_parser):
    """Defines the allowed arguments for plot function
    """
    parser = parent_parser.add_parser(
        "signals",
        help="Extracting SV related signals from bam files.",
    )

    parser.add_argument(
        "-n",
        "--samples",
        help="Space-delimited list of sample name.",
        type=str,
        default=None,
        nargs="+",
        required=False,
    )

    parser.add_argument(
        "-b",
        "--bamFiles",
        type=str,
        default=None,
        nargs="+",
        help="Space-delimited list of BAM/CRAM file names",
        required=False,
    )

    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default=None,
        nargs="+",
        help="A text file that stores two columns of information "
             + "sample name and bam path",
        required=False,
    )

    parser.add_argument(
        "-g",
        "--genome",
        help="Genome file for bam",
        type=str,
        required=False,
    )

    parser.add_argument(
        "-l",
        "--minSignal",
        default=25,
        type=int,
        help="The smallest SV signal in the bam file. default:25",
        required=False,
    )

    parser.add_argument(
        "-w",
        "--winSize",
        type=int,
        default=10000000,
        help="The size of the chromosome interval processed "
             + "by each process. default:10000000",
        required=False,
    )

    parser.add_argument(
        "--mutProces",
        type=int,
        default=10,
        help="The number of processes involved in counting "
             + "the number of Reads alignment in bam files. default:10",
        required=False,
    )

    parser.add_argument(
        "-o",
        "--outPath",
        type=str,
        default='.',
        help="The path to the folder where the output image is stored.",
        required=False,
    )
    parser.set_defaults(func=signals)


def call(parser, options, extra_args=None):
    pass
def args_call(parent_parser):
    """Defines the allowed arguments for plot function
    """
    parser = parent_parser.add_parser(
        "signals",
        help="Extracting SV related signals from bam files.",
    )

    parser.add_argument(
        "-n",
        "--samples",
        help="Space-delimited list of sample name.",
        type=str,
        default=None,
        nargs="+",
        required=False,
    )

    parser.add_argument(
        "-b",
        "--bedFiles",
        type=str,
        default=None,
        nargs="+",
        help="Space-delimited list of BAM/CRAM file names",
        required=False,
    )
    parser.add_argument(
        "-g",
        "--genome",
        help="Genome file for bam",
        type=str,
        required=False,
    )

    parser.add_argument(
        "-w",
        "--winSize",
        type=int,
        default=10000000,
        help="The size of the chromosome interval processed "
             + "by each process. default:10000000",
        required=False,
    )

    parser.add_argument(
        "--mutProces",
        type=int,
        default=10,
        help="The number of processes involved in counting "
             + "the number of Reads alignment in bam files. default:10",
        required=False,
    )

    parser.add_argument(
        "-o",
        "--outPath",
        type=str,
        default='.',
        help="The path to the folder where the output image is stored.",
        required=False,
    )
    parser.set_defaults(func=signals)
