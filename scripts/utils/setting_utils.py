from datetime import datetime
import os, sys
import shlex
import subprocess
from itertools import chain

def now_time(comment):
    now = datetime.now()
    nowtime = "{0:%Y-%m-%d %H:%M:%S}".format(now)
    print ('[' + nowtime + ']' + ' ' + comment)

#Require: import os
def get_absolute_path(filename):
    '''Get absolute path of a file'''
    return os.path.abspath(os.path.expanduser(filename))

class Bunch(object):
    '''Bunches parameters into a dictionary'''
    def __init__(self, adict):
        self.__dict__.update(adict)

def load_fasta(fasta_path):
    fasta_file = open(fasta_path,'r')
    fasta_dict = {}
    checker = []
    name = ''
    for line in fasta_file:
        line = line.rstrip()
        if line.startswith('#'):
            continue
        if line.startswith('>'):
            name = line[1:].strip()
            if not name in checker:
                checker.append(name)
            else:
                print ('ERROR: The same name exists =>' + name)
                sys.exit(1)
                continue
            fasta_dict[name] = ''
        else:
            trans_table = str.maketrans("ATGCatgcUu","AUGCAUGCUU")
            line_changed = line.translate(trans_table) #Translate DNA into RNA
            fasta_dict[name] += line_changed.upper()
    fasta_file.close()
    return fasta_dict

def reverse_complement(seq):
    if 'U' in seq or 'u' in seq:
        trans_table = str.maketrans("AUGCaugc","UACGUACG")
    else:
        trans_table = str.maketrans("ATGCatgc","UACGUACG")
    rev_seq = seq[-1::-1]
    rev_comp_seq = rev_seq.translate(trans_table)
    return rev_comp_seq

def rm_duplicate_list(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]

def list_counter(args):
    counter = {}
    for x in args:
        if x in counter: #if python2, has_key() is used.
            counter[x] += 1
        else:
            counter[x] = 1
    return counter #dict

def flatten(listOfLists):
    "Flatten one level of nesting"
    return list(chain.from_iterable(listOfLists))

def run_command(command_line):
    args = shlex.split(command_line)
    p = subprocess.Popen(args,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    out = stdout.decode('utf-8')
    err = stderr.decode('utf-8')
    exitcode = int(p.returncode)
    if not exitcode == 0:
        print(err)
        sys.exit(1)
    print (out)
    print (err)
    
