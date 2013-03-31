#!/usr/bin/env python
# encoding: utf-8
"""
File: get_reads_from_db.py
Author: Brant Faircloth

Created by Brant Faircloth on 01 November 2012 11:11 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import sqlite3
import argparse

import pdb

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Get reads out of the splitaake database""")
    parser.add_argument(
            "db",
            help="""Help text"""
        )
    parser.add_argument(
            "r1",
            help="""Help text"""
        )
    parser.add_argument(
            "r2",
            help="""Help text"""
        )
    parser.add_argument(
            "--outname",
            type=str,
            help="""The name to prepend to the outfiles""",
        )
    parser.add_argument(
            "--query",
            type=str,
            help="""The query string to use""",
        )
    return parser.parse_args()


def get_read_identifiers(cur, query):
    cur.execute(query)
    results = set([name[0] for name in cur.fetchall()])
    return results


def get_reads(ids, r1, r2, outname):
    #r1
    r1in = open(r1, 'rU')
    r2in = open(r2, 'rU')
    r1out = open('{}splitaake-query-output-r1.fastq'.format(outname), 'w')
    r2out = open('{}splitaake-query-output-r2.fastq'.format(outname), 'w')
    while True:
        r1split = [r1in.readline() for line in range(4)]
        r2split = [r2in.readline() for line in range(4)]
        if not r1split[1]:
            break
        r1name = r1split[0].split(' ')[0]
        r2name = r2split[0].split(' ')[0]
        #pdb.set_trace()
        assert r1name == r2name, "Files are not in the same order"
        if r1name in ids:
            r1out.write("{0}".format(''.join(r1split)))
            r2out.write("{0}".format(''.join(r2split)))
        else:
            pass
    for i in [r1in, r2in, r1out, r2out]:
        i.close()


def main():
    args = get_args()
    conn = sqlite3.connect(args.db)
    cur = conn.cursor()
    ids = get_read_identifiers(cur, args.query)
    get_reads(ids, args.r1, args.r2, args.outname)
    conn.close()

if __name__ == '__main__':
    main()