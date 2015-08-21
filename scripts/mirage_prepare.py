#!/usr/bin/env python

"""
MIRAGE Data Preparation: Data preparation for MIRAGE.

Created by Naoto Imamachi on 2015-04-23.

Usage:
  mirage_prepare.py <Preparation_type> [-i input_file] [-o output_file]

"""

import argparse
import os, sys
from datetime import datetime
import runpy
from parameter.common_parameters import common_parameters
import utils.setting_utils as utils

def greeting(parser=None):
    print ("MIRAGE Data Preparation Toolkit - Data preparation for MIRAGE")
    print ("-" * 20)
    if parser is not None:
        parser.print_help()

def main():
    parser = argparse.ArgumentParser(prog="mirage_prepare",description="MIRAGE preparation toolkit - Data preparation for MIRAGE")
    parser.add_argument('preparation_type',action='store',
                        choices=['bed_3UTR','mirbase_gff2bed','liftOver','phylop_score_prep','phastcons_prep','phylop_sizedown','phastcons_sizedown','phylop_score_list','phastcons_score_list','phastcons_score_R','phylop_score_R','refseq_pre','mirbase_pre','mirmark_result','cupid_result'],
                        help='Preparation Type: refseq_pre|mirbase_pre|mirmark_result|cupid_result|')
    parser.add_argument('-i','--input-file',action='store',dest='input_file',help='Input file: Specify a input file name and its path')
    parser.add_argument('-r','--reference-file',action='store',dest='reference_file',help='reference file: Specify a reference file name and its path')
    parser.add_argument('-a','--additional-file',action='store',dest='add_file',nargs=3,help='Additional_file: Specify 1-refseq_pre file dir, 2-mirbase_pre file dir, 3-error log dir)')
    parser.add_argument('-o','--ouput-file',action='store',dest='output_file',help='Output file: Specify a output file name and its path')
    args = parser.parse_args()

    #Start analysis - logging
    greeting()
    utils.now_time('MIRAGE Data Preparation starting...')

    #Parameter preparation
    prep_type = args.preparation_type
    if (args.input_file or args.output_file):
        if not (os.path.isfile(args.input_file)):
            utils.now_time("ERROR: InputFile does not exist...")
            sys.exit(1)
        if not (args.output_file):
            utils.now_time("ERROR: -o option are required...")
            sys.exit(1)
        custom_params = {}
        if prep_type == 'bed_3UTR':
            custom_params['bed_3UTR_input'] = args.input_file
            custom_params['bed_3UTR_output'] = args.output_file
            common_parameters.updata(custom_params)
            p = utils.Bunch(common_parameters) 
        elif prep_type == 'mirbase_gff2bed':
            custom_params['mirbase_gff2bed_input'] = args.input_file
            custom_params['mirbase_gff2bed_output'] = args.output_file
            common_parameters.update(custom_params)
            p = utils.Bunch(common_parameters)
        elif prep_type == 'liftOver':
            custom_params['liftover_input'] = args.input_file
            custom_params['liftover_output'] = args.output_file
            common_parameters.update(custom_params)
            p = utils.Bunch(common_parameters)
        elif prep_type == 'phylop_score_prep':
            custom_params['phylop_score_prep_input'] = args.input_file
            custom_params['phylop_score_prep_output'] = args.output_file
            common_parameters.update(custom_params)
            p = utils.Bunch(common_parameters)
        elif prep_type == 'phastcons_prep':
            custom_params['phastcons_prep_input'] = args.input_file
            custom_params['phastcons_prep_output'] = args.output_file
            common_parameters.update(custom_params)
            p = utils.Bunch(common_parameters)
        elif prep_type == 'phylop_sizedown':
            if args.reference_file:
                custom_params['phylop_sizedown_bed_input'] = args.reference_file
                custom_params['phylop_sizedown_score_input'] = args.input_file
                custom_params['phylop_sizedown_score_output'] = args.output_file
                common_parameters.update(custom_params)
                p = utils.Bunch(common_parameters)
            else:
                utils.now_time("ERROR: -r option is required...")
                sys.exit(1)
        elif prep_type == 'phastcons_sizedown':
            if args.reference_file:
                custom_params['phastcons_sizedown_bed_input'] = args.reference_file
                custom_params['phastcons_sizedown_score_input'] = args.input_file
                custom_params['phastcons_sizedown_score_output'] = args.output_file
                common_parameters.update(custom_params)
                p = utils.Bunch(common_parameters)
            else:
                utils.now_time("ERROR: -r option is required...")
                sys.exit(1)
        elif prep_type == 'phylop_score_list':
            if args.reference_file:
                custom_params['phylop_score_list_reference'] = args.reference_file
                custom_params['phylop_score_list_db_input'] = args.input_file
                custom_params['phylop_score_list_db_output'] = args.output_file
                common_parameters.update(custom_params)
                p = utils.Bunch(common_parameters)
            else:
                utils.now_time("ERROR: -r option is required...")
                sys.exit(1)
        elif prep_type == 'phastcons_score_list':
            if args.reference_file:
                custom_params['phastcons_score_list_reference'] = args.reference_file
                custom_params['phastcons_score_list_db_input'] = args.input_file
                custom_params['phastcons_score_list_db_output'] = args.output_file
                common_parameters.update(custom_params)
                p = utils.Bunch(common_parameters)
            else:
                utils.now_time("ERROR: -r option is required...")
                sys.exit(1)
        elif prep_type == 'phastcons_score_R':
            custom_params['phastcons_score_R_input'] = args.input_file
            custom_params['phastcons_score_R_output'] = args.output_file
            common_parameters.update(custom_params)
            p = utils.Bunch(common_parameters)
        elif prep_type == 'phylop_score_R':
            custom_params['phylop_score_R_input'] = args.input_file
            custom_params['phylop_score_R_output'] = args.output_file
            common_parameters.update(custom_params)
            p = utils.Bunch(common_parameters)
        elif prep_type == 'refseq_pre':
            custom_params['refseq_pre_input'] = args.input_file
            custom_params['refseq_pre_output'] = args.output_file
            common_parameters.update(custom_params)
            p = utils.Bunch(common_parameters)
        elif prep_type == 'mirbase_pre':
            custom_params['mirbase_pre_input'] = args.input_file
            custom_params['mirbase_pre_output'] = args.output_file
            p = utils.Bunch(common_parameters)
        elif prep_type == 'mirmark_result':
            if args.add_file:
                custom_params['refseq_pre_output'] = args.add_file[0]
                custom_params['mirbase_pre_output'] = args.add_file[1]
                custom_params['mirmark_pos'] = args.input_file
                custom_params['mirmark_output'] = args.output_file
                custom_params['mirmark_error'] = args.add_file[2]
                p = utils.Bunch(common_parameters)
            else:
                utils.now_time("ERROR: -a option is required...")
                sys.exit(1)
        elif prep_type == 'cupid_result':
            if args.add_file:
                custom_params['refseq_pre_output'] = args.add_file[0]
                custom_params['mirbase_pre_output'] = args.add_file[1]
                custom_params['cupid_pos'] = args.input_file
                custom_params['cupid_output'] = args.output_file
                custom_params['cupid_error'] = args.add_file[2]
                p = utils.Bunch(common_parameters)
            else:
                utils.now_time("ERROR: -a option is required...")
                sys.exit(1)

        else:
            utils.now_time("ERROR: Wrong preparation type...")
            sys.exit(1)
    elif not (args.input_file and args.output_file):
        if prep_type:
            p = utils.Bunch(common_parameters)
        else:
            utils.now_time("ERROR: Wrong preparation type...")
            sys.exit(1)
    else:
        utils.now_time("ERROR: -i and -o option are required...")
        sys.exit(1)

    #Preparation type
    if  prep_type == 'bed_3UTR':
        runpy.run_module('module.preparation.bed_3UTR',run_name="__main__",alter_sys=True)
    elif prep_type == 'mirbase_gff2bed':
        runpy.run_module('module.preparation.mirbase_gff2bed',run_name="__main__",alter_sys=True)
    elif prep_type == 'liftOver':
        runpy.run_module('module.preparation.liftOver',run_name="__main__",alter_sys=True)
    elif prep_type == 'phylop_score_prep':
        runpy.run_module('module.preparation.phylop_score_prep',run_name="__main__",alter_sys=True)
    elif prep_type == 'phastcons_prep':
        runpy.run_module('module.preparation.phastcons_score_prep',run_name="__main__",alter_sys=True)
    elif prep_type == 'phastcons_sizedown':
        runpy.run_module('module.preparation.phastcons_sizedown',run_name='__main__',alter_sys=True)
    elif prep_type == 'phylop_sizedown':
        runpy.run_module('module.preparation.phylop_sizedown',run_name="__main__",alter_sys=True)
    elif prep_type == 'phylop_score_list':
        runpy.run_module('module.preparation.phylop_score_list',run_name='__main__',alter_sys=True)
    elif prep_type == 'phastcons_score_list':
        runpy.run_module('module.preparation.phastcons_score_list',run_name='__main__',alter_sys=True)
    elif prep_type == 'phastcons_score_R':
        runpy.run_module('module.preparation.phastcons_score_R',run_name='__main__',alter_sys=True)
    elif prep_type == 'phylop_score_R':
        runpy.run_module('module.preparation.phylop_score_R',run_name='__main__',alter_sys=True)
    elif prep_type == 'refseq_pre':
        runpy.run_module('module.preparation.refseq_pre',run_name="__main__",alter_sys=True)
    elif prep_type == 'mirbase_pre':
        runpy.run_module('module.preparation.mirbase_pre',run_name="__main__",alter_sys=True)
    elif prep_type == 'mirmark_result':
        runpy.run_module('module.preparation.mirmark_result',run_name="__main__",alter_sys=True)
    elif prep_type == 'cupid_result':
        runpy.run_module('module.preparation.cupid_result',run_name="__main__",alter_sys=True)
    else:
        utils.now_time("ERROR: Wrong preparation type...")
        sys.exit(1)

if __name__ == '__main__':
    main()
