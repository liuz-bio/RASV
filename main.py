#!/usr/bin/env python
import argparse
import sys
import time

from __init__ import __version__
from plotSV import args_sample
from plotSV import args_vcf
from plotSV import args_signals

def main(args=None):
    start_time = time.time()
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser(
        prog="plotSV", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-v",
        "--version",
        help="Installed version",
        action="version",
        version="%(prog)s " + str(__version__),
    )
    sub = parser.add_subparsers(title="[sub-commands]", dest="command")
    sub.required = True

    args_sample(sub)
    args_signals(sub)
    args_vcf(sub)

    args,extra_args = parser.parse_known_args(args)
    args.func(parser, args, extra_args)

if __name__ == "__main__":
    print('aaa')
    sys.exit(main())

