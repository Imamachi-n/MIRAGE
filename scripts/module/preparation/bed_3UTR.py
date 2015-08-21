#!usr/bin/env python

import re
from parameter.common_parameters import common_parameters
import utils.setting_utils as utils

utils.now_time("bed_3UTR script starting...")
p = utils.Bunch(common_parameters)

def main():
    input_file = p.bed_3UTR_input
    output_file = p.bed_3UTR_output
    command_bed_3UTR = '../software/bed12to3UTRbed.sh ' + input_file + ' > ' + output_file
    print (command_bed_3UTR)
    utils.run_command(command_bed_3UTR)
    utils.now_time("bed_3UTR script was successfully finished!!")

if __name__ == '__main__':
    main()
