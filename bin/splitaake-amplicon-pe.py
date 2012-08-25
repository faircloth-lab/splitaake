#!/usr/bin/env python
# encoding: utf-8
"""
File: demuxipy-pe.py
Author: Brant Faircloth

Created by Brant Faircloth on 16 June 2012 15:06 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import pdb

#import os
import sys
import re
#import gzip
import time
#import numpy
#import string
#import cPickle
#import sqlite3
#import argparse
#import itertools
import ConfigParser

from jellyfish import levenshtein_distance as levenshtein
from itertools import islice
from multiprocessing import Process, Queue, JoinableQueue, cpu_count

from seqtools.sequence.fastq import FastqReader
from seqtools.sequence.fasta import FastaQualityReader
from seqtools.sequence.transform import DNA_reverse_complement

from splitaake.core import *
from splitaake import pe_db as db
#from demuxipy import pairwise2


class Tag:
    def __init__(self, typ):
        self.type = typ
        self.name = None
        self.seq = None
        self.start = None
        self.end = None
        self.match = None
        self.match_type = None

    def __str__(self):
        return "{0}({1})".format(self.__class__, self.__dict__['name'])

    def __repr__(self):
        return "<{0} instance at {1}>".format(self.__class__, hex(id(self)))


class Tags:
    '''Trimming, tag, and sequence data for individual reads'''
    def __init__(self, identifier):
        # super(Params, self).__init__()
        #assert isinstance(sequence, FastaSequence), \
        #    'The Record class must be instantiated with a FastaSequence object'
        # a biopython sequence object
        self.name = identifier
        self.outer = Tag('outer')
        self.f = Tag('forward')
        self.r = Tag('reverse')
        self.concat = Tag('concat')
        self.cluster = None
        self.cluster_name = None
        self.amb = False
        self.fprimer = Tag('fprimer')
        self.rprimer = Tag('rprimer')

    def __str__(self):
        return "{0}({1})".format(self.__class__, self.name)

    def __repr__(self):
        return "<{0} instance at {1}>".format(self.__class__, hex(id(self)))


def find_which_reads_have_tags(params, r1, r2, tags, p=False):
    fmatch, rmatch, forward, reverse = [None] * 4
    pair = [r1, r2]
    for index, read in enumerate(pair):
        fmatch = find_left_tag(
                read.sequence,
                params.sequence_tags.inners['None']['forward_regex'],
                params.sequence_tags.inners['None']['forward_string'],
                params.sequence_tags.inner_gap,
                params.sequence_tags.inner_len,
                params.inner_fuzzy,
                params.inner_errors
            )
        if fmatch:
            forward = read
            tags.f.name = read.identifier.split(' ')[1]
            tags.f.seq, tags.f.match_type, tags.f.start, tags.f.end, tags.f.match = fmatch
            if p:
                print tags.f.name
                print "forward = ", read.sequence[:20], fmatch
            #pdb.set_trace()
            trim_tags_from_reads(read, tags.f.start, tags.f.end)
            # remove the read from further consideration
            pair.pop(index)
            break
        else:
            forward = None
    for index, read in enumerate(pair):
        rmatch = find_left_tag(
                read.sequence,
                params.sequence_tags.inners['None']['reverse_regex'],
                params.sequence_tags.inners['None']['reverse_string'],
                params.sequence_tags.inner_gap,
                params.sequence_tags.inner_len,
                params.inner_fuzzy,
                params.inner_errors
            )
        if rmatch:
            reverse = read
            tags.r.name = read.identifier.split(' ')[1]
            tags.r.seq, tags.r.match_type, tags.r.start, tags.r.end, tags.r.match = rmatch
            if p:
                print tags.r.name
                print "reverse = ", read.sequence[:20], rmatch
            trim_tags_from_reads(read, tags.r.start, tags.r.end)
            break
        else:
            reverse = None
    if p:
        pdb.set_trace()
    return forward, reverse, tags


def get_cluster(params, tags):
    tags.cluster_name = "{},{}".format(tags.f.seq, tags.r.seq)
    tags.cluster = params.sequence_tags.cluster_map['None'][tags.cluster_name]
    return tags


def trim_tags_from_reads(read, start, end):
    return read.slice(end - start, len(read.sequence), False)


def trim_primer_sequences(params, forward, reverse, tags):
    fprimer = find_left_tag(
            forward.sequence,
            params.primers.f['regex'],
            params.primers.f['string'],
            params.primers.buff,
            params.primers.r['len'],
            params.primers.fuzzy,
            params.primers.errors
        )
    if fprimer:
        tags.fprimer.seq, tags.fprimer.match_type, tags.fprimer.start, tags.fprimer.end, tags.fprimer.match = fprimer
        trim_tags_from_reads(forward, tags.fprimer.start, tags.fprimer.end)
    rprimer = find_left_tag(
            reverse.sequence,
            params.primers.r['regex'],
            params.primers.r['string'],
            params.primers.buff,
            params.primers.r['len'],
            params.primers.fuzzy,
            params.primers.errors
        )
    if rprimer:
        tags.rprimer.seq, tags.rprimer.match_type, tags.rprimer.start, tags.rprimer.end, tags.rprimer.match = rprimer
        trim_tags_from_reads(reverse, tags.rprimer.start, tags.rprimer.end)
    return forward, reverse, tags


def singleproc(sequences, results, params, interval=1000, big_interval=10000):
    for pair in sequences:
        assert len(pair) == 2
        r1, r2 = pair
        identifiers = [i.identifier.split(' ')[0] for i in [r1, r2]]
        # make sure fastq headers are the same
        assert identifiers[0] == identifiers[1]
        # create a tag object to hold ID information
        tags = Tags(identifiers[0])
        if params.qual_trim:
            r1.read, r2.read = [r.trim(params.min_qual, False) for r in [r1, r2]]
        if params.drop:
            if 'N' in r1.sequence or 'N' in r2.sequence:
                tags.amb = True
                forward, reverse = None, None
        if not tags.amb and params.search == 'InnerCombinatorial':
            forward, reverse, tags = find_which_reads_have_tags(params, r1, r2, tags)
            if tags.f.name and tags.r.name:
                #pdb.set_trace()
                tags = get_cluster(params, tags)
                if params.primers.f and params.primers.r:
                    forward, reverse, tags = trim_primer_sequences(params, forward, reverse, tags)
        #write_results_out(tags, forward, reverse, fout, rout, count)
        results.put([forward, reverse, tags])
    #del sequences
    #return results


def multiproc(jobs, results, params):
    """locate linker sequences in a read, returning a record object"""
    while True:
        job = jobs.get()
        if job is None:
            break
        singleproc(job, results, params)


def split_every(n, iterable):
    i = iter(iterable)
    piece = list(islice(i, n))
    while piece:
        yield piece
        piece = list(islice(i, n))


def get_work(params):
    if params.r1 and params.r2:
        num_reads = get_sequence_count(params.r1, 'fastq')
        reads1 = FastqReader(params.r1)
        reads2 = FastqReader(params.r2)
        if params.num_procs > 1:
            merged = merge_fastq(reads1, reads2)
            work = split_every(10000, merged)
        else:
            work = merge_fastq(reads1, reads2)
        #pdb.set_trace()
    else:
        raise IOError("Cannot find r1 and r2 parameters in configuration file")
    return num_reads, work


def write_results_out(tags, forward, reverse, fout, rout, count):
    if tags.fprimer.seq and tags.rprimer.seq:
        fdist = levenshtein(tags.f.seq, tags.f.match)
        fout.write(">{0}_{1} {2} orig_bc={3} new_bc={4} bc_diff={5}\n{6}\n".format(
            tags.cluster,
            count,
            forward.read.identifier.split(' ')[0],
            tags.f.match,
            tags.f.seq,
            fdist,
            forward.read.sequence
            ))
        rdist = levenshtein(tags.r.seq, tags.r.match)
        rout.write(">{0}_{1} {2} orig_bc={3} new_bc={4} bc_diff={5}\n{6}\n".format(
            tags.cluster,
            count,
            reverse.read.identifier.split(' ')[0],
            tags.r.match,
            tags.r.seq,
            rdist,
            reverse.read.sequence
            ))


def main():
    """Main loop"""
    start_time = time.time()
    motd()
    args = get_args()
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    # build our configuration object w/ input params
    conf = ConfigParser.ConfigParser()
    conf.read(args.config)
    params = Parameters(conf)
    #pdb.set_trace()
    # create the db and tables, returning connection
    # and cursor
    conn, cur = db.create_db_and_new_tables(params.db)
    # get num reads and split up work
    print "Splitting reads into work units..."
    num_reads, work = get_work(params)
    print "There are {:,d} reads".format(num_reads)
    # give some indication of progress for longer runs
    if num_reads > 999:
        sys.stdout.write('Running...\n')
    #pdb.set_trace()
    # MULTICORE
    fout = open('forward-reads.fasta', 'w', 1)
    rout = open('reverse-reads.fasta', 'w', 1)
    if params.multiprocessing and params.num_procs > 1:
        jobs = Queue()
        results = JoinableQueue()
        # We're stacking groups of jobs on the work
        # Queue, conceivably to save the overhead of
        # placing them on there one-by-one.
        print "Adding jobs to work queue..."
        for unit in work:
            jobs.put(unit)
        print "There are {} jobs...".format(num_reads / 10000)
        # setup the processes for the jobs
        print "Starting {} workers...".format(params.num_procs)
        # start the worker processes
        [Process(target = multiproc, args=(jobs, results, params)).start()
            for i in xrange(params.num_procs)]
        # we're putting single results on the results Queue so
        # that the db can (in theory) consume them at
        # a rather consistent rate rather than in spurts
        #for unit in xrange(num_reads):
        #count = 0
        for unit in xrange(num_reads):
            forward, reverse, tags = results.get()
            rowid = db.insert_record_to_db(cur, tags)
            write_results_out(tags, forward, reverse, fout, rout, rowid)
            results.task_done()
            progress(rowid, 10000, 100000)
            #count += 1
        # make sure we put None at end of Queue
        # in an amount equiv. to num_procs
        for unit in xrange(params.num_procs):
            jobs.put(None)
        # join the results, so that they can finish
        results.join()
        # close up our queues
        jobs.close()
        results.close()

    # SINGLECORE
    else:
        # fake a multiprocessing queue, so stacking and accessing results
        # is identical.
        fake_queue = ListQueue()
        results = singleproc(work, fake_queue, params)
        for result in results:
            forward, reverse, tags = result
            rowid = db.insert_record_to_db(cur, tags)
            write_results_out(tags, forward, reverse, fout, rout, rowid)
            progress(rowid, 10000, 100000)
    fout.close()
    rout.close()
    conn.commit()
    cur.close()
    conn.close()
    end_time = time.time()
    pretty_end_time = time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print "\nEnded: {} (run time {} minutes)".format(pretty_end_time,
            round((end_time - start_time)/60, 3))

if __name__ == '__main__':
    main()