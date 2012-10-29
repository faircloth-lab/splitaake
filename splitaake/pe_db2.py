"""
File: db.py
Author: Brant Faircloth

Created by Brant Faircloth on 02 October 2011 15:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: Database creation and entry functions for demuxi.py 

"""

import os
import sys
import sqlite3
#try:
#    import cPickle as pickle
#except:
#    print "Using pickle instead of cPickle"
#    import pickle

#import pdb


def create_db_and_new_tables(db_name):
    conn = sqlite3.connect(db_name)
    cur = conn.cursor()
    cur.execute("PRAGMA foreign_keys = ON")
    try:
        cur.execute('''CREATE TABLE tags (
                id integer PRIMARY KEY AUTOINCREMENT,
                name text,
                f text,
                fseq text,
                fmatch text,
                ftype text,
                r text,
                rseq text,
                rmatch text,
                rtype text,
                fpseq text,
                fpmatch text,
                fptype text,
                rpseq text,
                rpmatch text,
                rptype text,
                foverrun text,
                roverrun text,
                individual text
            )'''
        )
        cur.execute("CREATE INDEX idx_sequence_cluster on tags(individual)")
    except sqlite3.OperationalError, e:
        #pdb.set_trace()
        if "already exists" in e[0]:
            answer = raw_input("\nDatabase already exists.  Overwrite [Y/n]? ")
            #pdb.set_trace()
            if answer == "Y" or answer == "YES":
                os.remove(db_name)
                conn, cur = create_db_and_new_tables(db_name)
            else:
                sys.exit()
        else:
            raise sqlite3.OperationalError, e
    return conn, cur


def insert_record_to_db(cur, dmux):
    cur.execute('''INSERT INTO tags (
            name,
            f,
            fseq,
            fmatch,
            ftype,
            r,
            rseq,
            rmatch,
            rtype,
            fpseq,
            fpmatch,
            fptype,
            rpseq,
            rpmatch,
            rptype,
            foverrun,
            roverrun,
            individual
        )
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''',
        (
            dmux.name,
            dmux.r1tag.name,
            dmux.r1tag.tag,
            dmux.r1tag.seq,
            dmux.r1tag.match_type,
            dmux.r2tag.name,
            dmux.r2tag.tag,
            dmux.r2tag.seq,
            dmux.r2tag.match_type,
            dmux.r1site.tag,
            dmux.r1site.seq,
            dmux.r1site.match_type,
            dmux.r2site.tag,
            dmux.r2site.seq,
            dmux.r2site.match_type,
            dmux.r1overrun.match,
            dmux.r2overrun.match,
            dmux.individual
        )
    )
    return cur.lastrowid
