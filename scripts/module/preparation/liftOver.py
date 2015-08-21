#!usr/bin/env python

import re
from parameter.common_parameters import common_parameters
import utils.setting_utils as utils

utils.now_time("liftOver script starting...")
p = utils.Bunch(common_parameters)

def main():
    input_file = p.liftover_input
    output_file = p.liftover_output
    error_file = output_file + '.error'
    command_liftover = '../software/liftOver ' + input_file + ' ../software/hg38ToHg19.over.chain ' + output_file + ' ' + error_file
    utils.run_command(command_liftover)
    utils.now_time("liftOver script was successfully finished!!")

if __name__ == '__main__':
    main()
