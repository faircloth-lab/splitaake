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

import re
import sys
import time
import ConfigParser
from multiprocessing import Process, Queue, JoinableQueue, cpu_count
from seqtools.sequence.fastq import FastqReader
from seqtools.sequence.fasta import FastaQualityReader
from seqtools.sequence.transform import DNA_reverse_complement
from splitaake.core import *
from splitaake import db


def get_args():
    """get arguments (config file location)"""
    parser = argparse.ArgumentParser(description = "splitaake.py:  sequence " + \
        "demultiplexing for hierarchically-tagged samples")
    parser.add_argument('config', 
            help="The input configuration file",
            action=FullPaths
        )
    parser.add_argument('--job-size',
            help="The number of sequence to add to each job",
            type=int,
            default=10000
        )
    return parser.parse_args()


def singleproc(sequences, results, params, interval=1000, big_interval=10000):
    for pair in sequences:
        r1, r2, header = split_and_check_reads(pair)
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
        if not dmux.drop:
            #pdb.set_trace()
            dmux = find_which_reads_have_tags(params, r1, r2, dmux)
            # if we've assigned reads (which means we've assigned tags)
            if dmux.r1tag.match and dmux.r2tag.match:
                try:
                    dmux.individual = params.tags.combo_lookup(dmux.r1tag.name, dmux.r2tag.name)
                except KeyError:
                    dmux.individual = None
                if dmux.individual is not None:
                    dmux = find_which_reads_have_sites(params, dmux)
                    if params.overruns.trim == True:
                        dmux = find_which_reads_have_overruns(params, dmux)
        results.put(dmux)
    return results


def multiproc(jobs, results, params):
    """locate linker sequences in a read, returning a record object"""
    while True:
        job = jobs.get()
        if job is None:
            break
        singleproc(job, results, params)
	del job


def get_work(params, job_size=10000):
    try:
        num_reads = get_sequence_count(params.reads.r1, 'fastq')
        reads1 = FastqReader(params.reads.r1)
        reads2 = FastqReader(params.reads.r2)
        if params.parallelism.cores > 1:
            merged = merge_fastq(reads1, reads2)
            work = split_every(job_size, merged)
        else:
            work = merge_fastq(reads1, reads2)
    except:
        raise IOError("Cannot find r1 and r2 parameters in configuration file")
    return num_reads, work


def write_results_out(cur, rowid, params, dmux):
    if dmux.r1tag.match and dmux.r1site.match and dmux.r2tag.match and dmux.r2site.match:
        if len(dmux.r1) >= params.quality.drop_len and len(dmux.r2) >= params.quality.drop_len:
            params.storage.output[dmux.individual]['R1'].write(dmux.r1)
            params.storage.output[dmux.individual]['R2'].write(dmux.r2)
            cur.execute("UPDATE tags SET written = 1 WHERE rowid = ?", (rowid,))


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
    # alert user we're dropping reads
    print "[WARN] Dropping all demultiplexed reads ≤ {0} bp long".format(params.quality.drop_len)
    # get num reads and split up work
    print "Splitting reads into work units..."
    num_reads, work = get_work(params, args.job_size)
    print "There are {:,d} reads".format(num_reads)
    # give some indication of progress for longer runs
    if num_reads > 999:
        sys.stdout.write('Running...\n')
    #pdb.set_trace()
    #r1out = open('r1-reads.fasta', 'w', 1)
    #r2out = open('r2-reads.fasta', 'w', 1)
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
        print "There are {} jobs...".format(num_reads / args.job_size)
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
            write_results_out(cur, rowid, params, dmux)
            results.task_done()
            progress(rowid, 10000, 100000)
            # delete obj
            del dmux
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
            write_results_out(cur, rowid, params, dmux)
            progress(rowid, 10000, 100000)
            count += 1
    params.storage.close()
    conn.commit()
    cur.close()
    conn.close()
    end_time = time.time()
    pretty_end_time = time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print "\nEnded: {} (run time {} minutes)".format(pretty_end_time,
            round((end_time - start_time)/60, 3))

if __name__ == '__main__':
    main()
