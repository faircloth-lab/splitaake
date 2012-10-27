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

from itertools import islice
from jellyfish import levenshtein_distance as levenshtein
from multiprocessing import Process, Queue, JoinableQueue, cpu_count

from seqtools.sequence.fastq import FastqReader
from seqtools.sequence.fasta import FastaQualityReader
from seqtools.sequence.transform import DNA_reverse_complement

from splitaake.core import *
from splitaake import pe_db2 as db


class Demux:
    '''Trimming, tag, and sequence data for individual reads'''
    def __init__(self, identifier):
        # super(Params, self).__init__()
        #assert isinstance(sequence, FastaSequence), \
        #    'The Record class must be instantiated with a FastaSequence object'
        # a biopython sequence object
        self.name = identifier
        self.r1 = None
        self.r2 = None
        self.individual = None
        self.drop = False
        self.r1site = SeqSearchResult('site')
        self.r2site = SeqSearchResult('site')
        self.r1tag = SeqSearchResult('tag')
        self.r2tag = SeqSearchResult('tag')

    def __str__(self):
        return "{0}({1})".format(self.__class__, self.name)

    def __repr__(self):
        return "<{0} instance at {1}>".format(self.__class__, hex(id(self)))


class SeqSearchResult:
    def __init__(self, typ):
        self.type = typ
        self.name = None
        self.seq = None
        self.tag = None
        self.start = None
        self.end = None
        self.match = False
        self.match_type = None

    def __str__(self):
        return "{0}({1})".format(self.__class__, self.__dict__['name'])

    def __repr__(self):
        return "<{0} instance at {1}>".format(self.__class__, hex(id(self)))


def find_which_reads_have_tags(params, r1, r2, dmux, p=False):
    pair = [r1, r2]
    for index, read in enumerate(pair):
        dmux.r1tag = find_tag(read.sequence, params.tags, 'r1', 'regex', dmux.r1tag)
        if dmux.r1tag.match:
            dmux.r1 = read
            #dmux.r1tag = forward
            trim_tags_from_reads(dmux.r1, dmux.r1tag.start, dmux.r1tag.end)
            if p:
                print dmux.r1.identifier
                print "regex forward: tag={0} seq={1} type={2} sequence={3}".format(dmux.r1tag.tag, dmux.r1tag.seq, dmux.r1tag.match_type, dmux.r1.sequence[:20])
            #pdb.set_trace()
            # remove the read from further consideration
            pair.pop(index)
            break
    for index, read in enumerate(pair):
        dmux.r2tag = find_tag(read.sequence, params.tags, 'r2', 'regex', dmux.r2tag)
        if dmux.r2tag.match:
            dmux.r2 = read
            #dmux.r2tag = reverse
            trim_tags_from_reads(dmux.r2, dmux.r2tag.start, dmux.r2tag.end)
            if p:
                print dmux.r2.identifier
                print "regex reverse: tag={0} seq={1} type={2} sequence={3}".format(dmux.r2tag.tag, dmux.r2tag.seq, dmux.r2tag.match_type, dmux.r2.sequence[:20])
            pair.pop(index)
            break
    if not dmux.r1tag.match:
        for index, read in enumerate(pair):
            dmux.r1tag = find_tag(read.sequence, params.tags, 'r1', 'fuzzy', dmux.r1tag)
            if dmux.r1tag.match:
                dmux.r1 = read
                #dmux.r1tag = forward
                trim_tags_from_reads(dmux.r1, dmux.r1tag.start, dmux.r1tag.end)
                if p:
                    print dmux.r1.identifier
                    print "fuzzy forward: tag={0} seq={1} type={2} sequence={3}".format(dmux.r1tag.tag, dmux.r1tag.seq, dmux.r1tag.match_type, dmux.r1.sequence[:20])
                #pdb.set_trace()
                # remove the read from further consideration
                pair.pop(index)
                break
    if not dmux.r2tag.match:
        for index, read in enumerate(pair):
            dmux.r2tag = find_tag(read.sequence, params.tags, 'r2', 'fuzzy', dmux.r2tag)
            if dmux.r2tag.match:
                dmux.r2 = read
                #dmux.r2tag = reverse
                trim_tags_from_reads(dmux.r2, dmux.r2tag.start, dmux.r2tag.end)
                if p:
                    print dmux.r2.identifier
                    print "fuzzy reverse: tag={0} seq={1} type={2} sequence={3}".format(dmux.r2tag.tag, dmux.r2tag.seq, dmux.r2tag.match_type, dmux.r2.sequence[:20])
                break
    #if p:
    #    pdb.set_trace()
    return dmux


def get_cluster(params, tags):
    tags.cluster_name = "{},{}".format(tags.f.seq, tags.r.seq)
    tags.cluster = params.sequence_tags.cluster_map['None'][tags.cluster_name]
    return tags


def trim_tags_from_reads(read, start, end):
    return read.slice(end - start, len(read.sequence), False)


def find_which_reads_have_sites(params, dmux):
    find_site(dmux.r1.sequence, params.site, 'r1', dmux.r1site)
    if dmux.r1site.match:
        trim_tags_from_reads(dmux.r1, dmux.r1site.start, dmux.r1site.end)
    find_site(dmux.r2.sequence, params.site, 'r2', dmux.r2site)
    if dmux.r2site.match:
        trim_tags_from_reads(dmux.r2, dmux.r2site.start, dmux.r2site.end)
    return dmux


def split_and_check_reads(pair):
    assert len(pair) == 2
    r1, r2 = pair
    identifiers = [i.identifier.split(' ')[0] for i in [r1, r2]]
    # make sure fastq headers are the same
    assert identifiers[0] == identifiers[1]
    return r1, r2, identifiers[0]


def singleproc(sequences, results, params, interval=1000, big_interval=10000):
    for pair in sequences:
        r1, r2, header = split_and_check_reads(pair)
        #print "\n"
        #print r1.sequence, r1.quality
        #print r2.sequence, r2.quality
        # create a Demux metadata object to hold ID information
        dmux = Demux(header)
        # quality trim the ends of reads before proceeding
        if params.quality.trim:
            r1.read, r2.read = [r.trim(params.quality.min, False) for r in [r1, r2]]
        # drop anything with an N in either read after quality trimming
        if params.quality.drop_n:
            if 'N' in r1.sequence or 'N' in r2.sequence:
                dmux.drop = True
        # keep all reads when params.quality.drop_n == FALSE
        else:
            dmux.drop = False
        if not dmux.drop:
            dmux = find_which_reads_have_tags(params, r1, r2, dmux)
            # if we've assigned reads (which means we've assigned tags)
            if dmux.r1tag.match and dmux.r2tag.match:
                dmux.individual = params.tags.combo_lookup(dmux.r1tag.name, dmux.r2tag.name)
                #print dmux.individual
                #pdb.set_trace()
                dmux = find_which_reads_have_sites(params, dmux)
                #if dmux.r1site and dmux.r2site:
                #    good += 1
                #try:
                #    print "r1 primer= ", dmux.r1site.tag, dmux.r1site.match_type, dmux.r1.sequence[:20]
                #except:
                #    pass
                #try:
                #    print "r2 primer= ", dmux.r2site.tag, dmux.r2site.match_type, dmux.r2.sequence[:20]
                #except:
                #    pass
                #pdb.set_trace()
        results.put(dmux)
    return results


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
    try:
        num_reads = get_sequence_count(params.reads.r1, 'fastq')
        reads1 = FastqReader(params.reads.r1)
        reads2 = FastqReader(params.reads.r2)
        if params.parallelism.cores > 1:
            merged = merge_fastq(reads1, reads2)
            work = split_every(10000, merged)
        else:
            work = merge_fastq(reads1, reads2)
        #pdb.set_trace()
    except:
        raise IOError("Cannot find r1 and r2 parameters in configuration file")
    return num_reads, work


def write_results_out(r1out, r2out, dmux, rowid):
    if dmux.r1tag.match and dmux.r1site.match and dmux.r2tag.match and dmux.r2site.match:
        #pdb.set_trace()
        #fdist = levenshtein(
        #        dmux.r1tag.seq[dmux.r1tag.start:dmux.r1tag.end],
        #        dmux.r1tag.tag[dmux.r1tag.start:dmux.r1tag.end]
        #    )
        r1out.write(">{0}_{1} {2} expected_index={3} observed_index={4} match_type={5}\n{6}\n".format(
                dmux.individual,
                rowid,
                dmux.r1.identifier.split(' ')[0],
                dmux.r1tag.tag,
                dmux.r1tag.seq,
                dmux.r1tag.match_type,
                dmux.r1.sequence
            ))
        #rdist = levenshtein(tags.r.seq, tags.r.match)
        r2out.write(">{0}_{1} {2} expected_index={3} observed_index={4} match_type={5}\n{6}\n".format(
                dmux.individual,
                rowid,
                dmux.r2.identifier.split(' ')[0],
                dmux.r2tag.tag,
                dmux.r2tag.seq,
                dmux.r2tag.match_type,
                dmux.r2.sequence
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
    # create the db and tables, returning connection
    # and cursor
    if params.db.create:
        conn, cur = db.create_db_and_new_tables(params.db.name)
    # get num reads and split up work
    print "Splitting reads into work units..."
    num_reads, work = get_work(params)
    print "There are {:,d} reads".format(num_reads)
    # give some indication of progress for longer runs
    if num_reads > 999:
        sys.stdout.write('Running...\n')
    #pdb.set_trace()
    r1out = open('r1-reads.fasta', 'w', 1)
    r2out = open('r2-reads.fasta', 'w', 1)
    # MULTICORE
    if params.parallelism.cores > 1:
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
        print "Starting {} workers...".format(params.parallelism.cores)
        # start the worker processes
        [Process(target=multiproc, args=(jobs, results, params)).start()
            for i in xrange(params.parallelism.cores)]
        # we're putting single results on the results Queue so
        # that the db can (in theory) consume them at
        # a rather consistent rate rather than in spurts
        #for unit in xrange(num_reads):
        count = 0
        for unit in xrange(num_reads):
            dmux = results.get()
            rowid = db.insert_record_to_db(cur, dmux)
            write_results_out(r1out, r2out, dmux, rowid)
            results.task_done()
            progress(rowid, 10000, 100000)
            count += 1
        # make sure we put None at end of Queue
        # in an amount equiv. to num_procs
        for unit in xrange(params.parallelism.cores):
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
        count = 0
        for dmux in results:
            rowid = db.insert_record_to_db(cur, dmux)
            write_results_out(r1out, r2out, dmux, rowid)
            progress(rowid, 10000, 100000)
            count += 1
    r1out.close()
    r2out.close()
    conn.commit()
    cur.close()
    conn.close()
    end_time = time.time()
    pretty_end_time = time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print "\nEnded: {} (run time {} minutes)".format(pretty_end_time,
            round((end_time - start_time)/60, 3))

if __name__ == '__main__':
    main()