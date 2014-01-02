#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 02 January 2014 12: 32, PST (-0800)
"""


from __future__ import absolute_import
import os
import glob
import gzip
import shutil
from splitaake.split import many
from splitaake.fastq import readfq
from splitaake.core import mkdir_p


class MockArgs:
    def __init__(self):
        s = os.path.dirname(os.path.realpath(__file__))
        self.reads = os.path.join(s, "data")
        self.tagmap = os.path.join(s, "test_split_tags.config")
        self.output = os.path.join(s, "trash")
        self.section = "tags"
        self.quiet = True
        self.no_correct = False
        self.hamming = False
        self.min_qual = 10
        self.min_mean_qual = 20
        self.tag_length = None
        mkdir_p(self.output)


class TestFunction:
    def test_main(self):
        args = MockArgs()
        many.main(args, None)
        # count sequences in each file
        obs = {}
        exp = {
            "bfidt-001_R1.fastq.gz": 6,
            "bfidt-001_R2.fastq.gz": 6,
            "bfidt-006_R1.fastq.gz": 22,
            "bfidt-006_R2.fastq.gz": 22,
            "bfidt-009_R1.fastq.gz": 14,
            "bfidt-009_R2.fastq.gz": 14,
            "bfidt-010_R1.fastq.gz": 27,
            "bfidt-010_R2.fastq.gz": 27,
            "bfidt-012_R1.fastq.gz": 13,
            "bfidt-012_R2.fastq.gz": 13,
            "bfidt-013_R1.fastq.gz": 12,
            "bfidt-013_R2.fastq.gz": 12,
            "bfidt-014_R1.fastq.gz": 6,
            "bfidt-014_R2.fastq.gz": 6,
            "bfidt-015_R1.fastq.gz": 16,
            "bfidt-015_R2.fastq.gz": 16,
            "bfidt-017_R1.fastq.gz": 16,
            "bfidt-017_R2.fastq.gz": 16,
            "bfidt-018_R1.fastq.gz": 16,
            "bfidt-018_R2.fastq.gz": 16,
            "bfidt-019_R1.fastq.gz": 9,
            "bfidt-019_R2.fastq.gz": 9,
            "bfidt-020_R1.fastq.gz": 16,
            "bfidt-020_R2.fastq.gz": 16,
            "bfidt-023_R1.fastq.gz": 22,
            "bfidt-023_R2.fastq.gz": 22,
            "bfidt-026_R1.fastq.gz": 3,
            "bfidt-026_R2.fastq.gz": 3,
            "bfidt-030_R1.fastq.gz": 9,
            "bfidt-030_R2.fastq.gz": 9,
            "bfidt-032_R1.fastq.gz": 22,
            "bfidt-032_R2.fastq.gz": 22,
            "bfidt-033_R1.fastq.gz": 13,
            "bfidt-033_R2.fastq.gz": 13,
            "bfidt-034_R1.fastq.gz": 24,
            "bfidt-034_R2.fastq.gz": 24,
            "bfidt-036_R1.fastq.gz": 14,
            "bfidt-036_R2.fastq.gz": 14,
            "bfidt-038_R1.fastq.gz": 14,
            "bfidt-038_R2.fastq.gz": 14,
            "bfidt-039_R1.fastq.gz": 28,
            "bfidt-039_R2.fastq.gz": 28,
            "bfidt-041_R1.fastq.gz": 29,
            "bfidt-041_R2.fastq.gz": 29,
            "bfidt-042_R1.fastq.gz": 17,
            "bfidt-042_R2.fastq.gz": 17,
            "bfidt-043_R1.fastq.gz": 26,
            "bfidt-043_R2.fastq.gz": 26,
            "bfidt-044_R1.fastq.gz": 12,
            "bfidt-044_R2.fastq.gz": 12,
            "bfidt-045_R1.fastq.gz": 29,
            "bfidt-045_R2.fastq.gz": 29,
            "bfidt-101_R1.fastq.gz": 25,
            "bfidt-101_R2.fastq.gz": 25,
            "bfidt-105_R1.fastq.gz": 21,
            "bfidt-105_R2.fastq.gz": 21,
            "bfidt-110_R1.fastq.gz": 29,
            "bfidt-110_R2.fastq.gz": 29,
            "bfidt-113_R1.fastq.gz": 11,
            "bfidt-113_R2.fastq.gz": 11,
            "bfidt-117_R1.fastq.gz": 17,
            "bfidt-117_R2.fastq.gz": 17,
            "bfidt-125_R1.fastq.gz": 0,
            "bfidt-125_R2.fastq.gz": 0,
            "lowqual-R1.fastq.gz": 64,
            "lowqual-R2.fastq.gz": 64,
            "lowqual-R3.fastq.gz": 64,
            "unknown-R1.fastq.gz": 2398,
            "unknown-R2.fastq.gz": 2398,
            "unknown-R3.fastq.gz": 2398
        }
        for file in glob.glob(os.path.join(args.output, "*.gz")):
            fq = readfq(gzip.open(file, 'rb'))
            obs[os.path.basename(file)] = len(list(fq))
        for obsk, obsv in obs.iteritems():
            assert exp[obsk] == obsv
        shutil.rmtree(args.output)
