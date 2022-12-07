#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
#
#
######################################################################
# @author      : bjf79 (bjf79@midway2-login2.rcc.local)
# @file        : DetailedAnnotationOfIntrons
# @created     : Tuesday Dec 06, 2022 15:38:10 CST
#
# @description : 
######################################################################

import sys
import pysam
import os

# I like to script and debug with an interactive interpreter. If using
# interactive interpreter to script quick args with sys.argv, parse_args with
# hardcoded args below
if hasattr(sys, 'ps1'):
    sys.argv = ["", "Arg1" ,"Arg2"]

_, MyArg1, MyArg2 = sys.argv

def main():
    pass

if __name__ == '__main__':
    main()

tbx = pysam.TabixFile("ReferenceGenome/Annotations/gencode.v34.chromasomal.basic.annotation.gtf.gz")

for row in tbx.fetch("chr1", 1000, 20000, parser=pysam.asGTF()):
    print(row.attributes)
