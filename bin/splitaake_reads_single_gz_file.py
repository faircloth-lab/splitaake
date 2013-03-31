"""
File: dmux.py
Author: Brant Faircloth

Created by Brant Faircloth on 08 February 2012 17:02 PST (-0800)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import numpy
import argparse
import ConfigParser
from itertools import izip
from seqtools.sequence import fastq


from jellyfish import hamming_distance as hamming


import pdb

def get_args():
    parser = argparse.ArgumentParser(description='demux')
    parser.add_argument('r1', help='The input read1 fastq')
    parser.add_argument('r2', help='The input index read fastq')
    parser.add_argument('r3', help='The input read2 fastq')
    parser.add_argument('tagmap', help='A conf file mapping name:tag sequence')
    parser.add_argument('--section', help='''If multisection conf, use this
    section''')
    parser.add_argument('--no-correct', help='''Do not correct certain tags
        in the [no correct] section of the tagmap config file''',
        action='store_true', default = False)
    args = parser.parse_args()
    assert len(set([args.r1,args.r2,args.r3])) == 3, \
            """Fastq file names appear to be the same"""
    return args

#def hamming(s1, s2):
#    return sum([ch1 != ch2 for ch1, ch2 in zip(s1, s2)])

def distance(a,b):
    """Pure python version to compute the levenshtein distance between a and b.
    The Levenshtein distance includes insertions, deletions, substitutions; 
    unlike the Hamming distance, which is substitutions only.

    ORIGINALLY FROM:  http://hetland.org/coding/python/levenshtein.py

    """
    n, m = len(a), len(b)
    if n > m:
        # Make sure n <= m, to use O(min(n,m)) space
        a,b = b,a
        n,m = m,n
    current = range(n+1)
    for i in range(1,m+1):
        previous, current = current, [i]+[0]*n
        for j in range(1,n+1):
            add, delete = previous[j]+1, current[j-1]+1
            change = previous[j-1]
            if a[j-1] != b[i-1]:
                change = change + 1
            current[j] = min(add, delete, change)
    return current[n]

class Tags:
    def __init__(self, conf, section, no_correct = False):
        config = ConfigParser.ConfigParser()
        config.read(conf)
        tags = config.items(section)
        if no_correct:
            nc = config.items('no correct')
            self.__get_no_correct_set(nc)
        else:
            self.no_correct = None
        self.__name_to_tag(tags)
        self.__tag_to_name()
        self.__get_min_length()
        self.seqs = self.seq_dict.keys()
        self.count = len(self.name_dict)

    def __name_to_tag(self, tags):
        self.name_dict = dict(tags)

    def __tag_to_name(self):
        self.seq_dict = {k:v for v,k in self.name_dict.iteritems()}

    def __get_min_length(self):
        l = set([len(s) for s in self.seq_dict.keys()])
        assert len(l) == 1
        self.length = l.pop()

    def __get_no_correct_set(self, nc):
        if nc:
            self.no_correct = set([v for k,v in nc])
        else:
            self.no_correct = None

    def _create_zip_files(self, out = os.getcwd()):
        pth = os.path.join(out, 'dmux')
        if not os.path.exists(pth):
            os.makedirs(pth)
        self.files = {}
        for name, seq in self.name_dict.iteritems():
            name = "{}.fastq.gz".format(name)
            self.files[seq] = fastq.FasterFastqWriter(os.path.join(pth, name))
        unknown = ["unknown-r1.fastq.gz", "unknown-r2.fastq.gz", "unknown-r3.fastq.gz"]
        for f in unknown:
            self.files[f.split('.')[0]] = fastq.FasterFastqWriter(os.path.join(pth, f))
        lowqual = ["lowqual-r1.fastq.gz", "lowqual-r2.fastq.gz", "lowqual-r3.fastq.gz"]
        for l in lowqual:
            self.files[l.split('.')[0]] = fastq.FasterFastqWriter(os.path.join(pth, l))

    def _close_zip_files(self):
        for f in self.files.values():
            f.close()

def change_read_num(sequence, read, index):
    si = sequence[0].split(' ')
    sil = list(si[1])
    if read == 2:
        sil[0] = '2'
    new_si = "{0}{1}".format(''.join(sil), index)
    new_head = ' '.join([si[0], new_si])
    sl = list(sequence)
    sl[0] = new_head
    #pdb.set_trace()
    return tuple(sl)

def main():
    args = get_args()
    #pdb.set_trace()
    # setu tags object
    tags = Tags(args.tagmap, args.section, args.no_correct)
    #pdb.set_trace()
    # create output files
    tags._create_zip_files()
    #pdb.set_trace()
    # vectorize the hamming function to we only call once per read
    vector_hamming = numpy.vectorize(hamming)
    #vector_levenshtein = numpy.vectorize(distance)
    # read all of our files into fastq iterators
    read1 = fastq.FasterFastqReader(args.r1)
    index = fastq.FasterFastqReader(args.r2)
    read2 = fastq.FasterFastqReader(args.r3)
    #pdb.set_trace()
    read = 0
    for r1,i,r2 in izip(read1, index, read2):
        if read % 10000 == 0:
            print read
        if i[0].split(' ')[1].split(':')[1] == 'N':
            # get index sequence diff from tags
            dist = vector_hamming(tags.seqs, i[2])
            positions = numpy.where(dist == 0)[0]
            if positions is None and i[2] not in tags.no_correct:
                positions = numpy.where(dist == 1)[0]
            if (positions is not None and len(positions)) == 1:
                # assert headers match
                assert r1[0].split(' ')[0] \
                        == r2[0].split(' ')[0], "Header mismatch"
                # get tag for match
                match = tags.seqs[positions[0]]
                # write to output
                r1 = change_read_num(r1, 1, match)
                r2 = change_read_num(r2, 2, match)
                tags.files[match].write(r1) 
                tags.files[match].write(r2)
            else:
                #pdb.set_trace()
                tags.files['unknown-r1'].write(r1)
                tags.files['unknown-r2'].write(i)
                tags.files['unknown-r3'].write(r2)
        else:
            tags.files['lowqual-r1'].write(r1)
            tags.files['lowqual-r2'].write(i)
            tags.files['lowqual-r3'].write(r2)

        read += 1
    tags._close_zip_files()

if __name__ == '__main__':
    main()
