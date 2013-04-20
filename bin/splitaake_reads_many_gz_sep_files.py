"""
File: dmux.py
Author: Brant Faircloth

Created by Brant Faircloth on 08 February 2012 17:02 PST (-0800)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import glob
import numpy
import argparse
import ConfigParser
from itertools import izip
from seqtools.sequence import fastq
from collections import defaultdict

from jellyfish import levenshtein_distance as levenshtein
from jellyfish import hamming_distance as hamming


import pdb


def get_args():
    parser = argparse.ArgumentParser(description='demux')
    parser.add_argument('--reads', 
        required = True,
        help='The input fastq files'
    )
    parser.add_argument(
        '--tagmap',
        required = True,
        help='A conf file mapping name:tag sequence'
    )
    parser.add_argument('--output',
        required = True,
        help='The path to the output directory'
    )
    parser.add_argument('--section',
        default='tags',
        help="If multisection conf, use this section"
    )
    parser.add_argument('--no-correct',
        action='store_true',
        default=False,
        help="Do not correct certain tags " +
        "in the [no correct] section of the tagmap config file",
        )
    parser.add_argument('--tag-length',
        type=int,
        default=None,
        help="Length of tags if tags are < read length"
    )
    parser.add_argument('--min-qual',
        type=int,
        default=20,
        help="The minimum average quality (Q) to accept (overall min Q must be >= 10)"
    )
    args = parser.parse_args()
    return args


class Tags:
    def __init__(self, conf, section, no_correct=False):
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
        self.seq_dict = {k:v for v, k in self.name_dict.iteritems()}

    def __get_min_length(self):
        l = set([len(s) for s in self.seq_dict.keys()])
        assert len(l) == 1
        self.length = l.pop()

    def __get_no_correct_set(self, nc):
        if nc:
            self.no_correct = set([v for k, v in nc])
        else:
            self.no_correct = None

    def create_zip_files(self, pth):
        if not os.path.exists(pth):
            os.makedirs(pth)
        self.files = defaultdict(dict)
        for name, seq in self.name_dict.iteritems():
            for i in [1, 2]:
                r = "{0}_R{1}.fastq.gz".format(name, i)
                self.files[seq][i] = fastq.FasterFastqWriter(os.path.join(pth, r))
        for name in ['unknown', 'lowqual']:
            for i in [1, 2, 3]:
                r = "{0}-R{1}.fastq.gz".format(name, i)
                self.files[name][i] = fastq.FasterFastqWriter(os.path.join(pth, r))

    def close_zip_files(self):
        #pdb.set_trace()
        for f in self.files.values():
            for i in [1, 2]:
                f[i].close()


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


def get_quality(string, ascii_min=33, ascii_max=126, quality_min=0, quality_max=93):
    qa = numpy.fromstring(string, dtype='uint8')
    sanger = qa + ascii_min - quality_min
    # restrict to max
    return sanger.clip(max=quality_max)


def main():
    args = get_args()
    # setu tags object
    tags = Tags(args.tagmap, args.section, args.no_correct)
    # create output files
    tags.create_zip_files(args.output)
    # vectorize the levenshtein function so we only call once per read
    distance = numpy.vectorize(levenshtein)
    read = 0
    for f in glob.glob(os.path.join(args.r1, '*_R1_*')):
        print "Working on ", f
        read1 = fastq.FasterFastqReader(f)
        # get basename of R1 file
        read1_basename = os.path.basename(f)
        index = fastq.FasterFastqReader(os.path.join(args.r1, read1_basename.replace('R1', 'R2')))
        read2 = fastq.FasterFastqReader(os.path.join(args.r1, read1_basename.replace('R1', 'R3')))
        # read all of our files into fastq iterators
        for r1, i, r2 in izip(read1, index, read2):
            if read % 10000 == 0:
                print read
            if i[0].split(' ')[1].split(':')[1] == 'N':
                # get index sequence diff from tags
                dist = distance(tags.seqs, i[2])
                quality = get_quality(i[3])
                positions = numpy.where(dist == 0)[0]
                if positions is None and min(quality) >= 20 and i[2] not in tags.no_correct:
                    positions = numpy.where(dist == 1)[0]
                if (positions is not None) and (len(positions) == 1) and min(quality) >= 20:
                    # assert headers match
                    assert r1[0].split(' ')[0] == r2[0].split(' ')[0], "Header mismatch"
                    # get tag for match
                    match = tags.seqs[positions[0]]
                    # write to output
                    r1 = change_read_num(r1, 1, match)
                    r2 = change_read_num(r2, 2, match)
                    tags.files[match][1].write(r1)
                    tags.files[match][2].write(r2)
                else:
                    #pdb.set_trace()
                    tags.files['unknown'][1].write(r1)
                    tags.files['unknown'][2].write(i)
                    tags.files['unknown'][3].write(r2)
            else:
                tags.files['lowqual'][1].write(r1)
                tags.files['lowqual'][2].write(i)
                tags.files['lowqual'][3].write(r2)
            read += 1
        read1.close()
        index.close()
        read2.close()
    tags.close_zip_files()

if __name__ == '__main__':
    main()
