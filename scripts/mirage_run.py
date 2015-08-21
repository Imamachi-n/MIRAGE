#!/usr/bin/env python

"""
MIRAGE: Comprehensive miRNA target prediction pipeline.

Created by Naoto Imamachi on 2015-04-23.
Copyright (c) 2015 Naoto Imamachi. All rights reserved.
Updated and maintained by Naoto Imamachi since Apr 2015.

Usage:
  mirage.py <analysis_type> <miRNA.fasta> <targetRNA.fasta> [options]

"""

import os, sys
import argparse
import runpy
import utils.setting_utils as utils
from parameter.common_parameters import common_parameters

def greeting(parser=None):
    print ("MIRAGE v.0.1.0-beta - Comprehensive miRNA target prediction pipeline")
    print ("-" * 20)
    if parser is not None:
        parser.print_help()

def main():
    parser = argparse.ArgumentParser(prog='mirage',description='MIRAGE - Comprehensive miRNA target prediction pipeline')
    parser.add_argument('analysis_type',action='store',help='Analysis_type: Choose estimation or prediction',choices=['estimation','prediction'])
    parser.add_argument('mirna_fasta',action='store',help='miRNA fasta file: Specify miRNA fasta file to use the analysis')
    parser.add_argument('targetrna_fasta',action='store',help='TargetRNA fasta file: Specify TargetRNA fasta file to use the analysis')
    parser.add_argument('-m','--mirna-conservation-score-file',action='store',dest='mirna_conservation',help='Conservation score file about miRNA: Specify your conservation score db file. MIRAGE preparetion toolkits enables you to make the score files about TargetRNA or miRNA bed files.')
    parser.add_argument('-t','--targetrna-conservation-score-file',action='store',dest='targetrna_conservation',help='Conservation score file about TargetRNA: Specify your conservation score db file. MIRAGE preparetion toolkits enables you to make the score files about TargetRNA or miRNA bed files.')

    args = parser.parse_args()

    #Start analysis - logging
    greeting()
    utils.now_time("MIRAGE miRNA target prediction starting...")
    analysis_type = args.analysis_type
    mirna_fasta_path = args.mirna_fasta
    targetrna_fasta_path = args.targetrna_fasta
    mirna_conservation_score = args.mirna_conservation
    targetrna_conservation_score = args.targetrna_conservation

    #Check fasta files
    if not os.path.isfile(mirna_fasta_path):
        print ("Error: miRNA fasta file does not exist...")
        sys.exit(1)
    if not os.path.isfile(targetrna_fasta_path):
        print ("Error: TargetRNA fasta file does not exist...")

    #Check conservation score db files
    #if 

    #parameters
    param = dict(
        MIRNA_FASTA_PATH = mirna_fasta_path,
        TARGETRNA_FASTA_PATH = targetrna_fasta_path,
    )
    common_parameters.update(param)
    p = utils.Bunch(common_parameters)
    print ('miRNA_Fasta_file: ' + p.MIRNA_FASTA_PATH,end="\n")
    print ('TargetRNA_Fasta_file: ' + p.TARGETRNA_FASTA_PATH,end="\n")

    '''
    mirna_dict = utils.load_fasta(mirna_fasta_path)
    #print (mirna_dict['hsa-miR-34b-5p|MIMAT0000685'],end="\n")
    #print (mirna_dict['hsa-miR-20a-5p|MIMAT0000075'],end="\n")
    targetrna_dict = utils.load_fasta(targetrna_fasta_path)
    #print (targetrna_dict['NM_000594'],end="\n")
    #print (targetrna_dict['NM_030938'],end="\n")
    
    query_mirna.update(mirna_dict)
    print (query_mirna)
    mirna = utils.Bunch(query_mirna)
    query_targetrna.update(targetrna_dict)
    targetrna = utils.Bunch(query_targetrna)
    if hasattr (mirna,'hsa-miR-34b-5p|MIMAT0000685'):
        print ("OK!!")
        print (mirna.items())
        sys.exit(0)
    else:
        print ("Error...")
        sys.exit(1)
    #test = targetrna.'NM_000594'
    #print (test,end="\n")
    #sys.exit(0)
    '''

    #runpy - choose analysis type
    if analysis_type == 'estimation':
        runpy.run_module('module.estimate',run_name="__main__",alter_sys=True)
    elif analysis_type == 'prediction':
        runpy.run_module('module.predict',run_name="__main__",alter_sys=True)
    else:
        print ('Error: Analysis type is wrong...')
        sys.exit(1)

if __name__ == '__main__':
    main()
