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
try:
    import cPickle as pickle
except:
    print "Using pickle instead of cPickle"
    import pickle

def create_db_and_new_tables(db_name):
    conn = sqlite3.connect(db_name)
    cur = conn.cursor()
    cur.execute("PRAGMA foreign_keys = ON")
    try:
        cur.execute('''CREATE TABLE tags (
                id integer PRIMARY KEY AUTOINCREMENT,
                name text,
                outer text,
                outer_seq text,
                outer_match text,
                outer_method text,
                inner text,
                inner_seq text,
                inner_match text,
                inner_method text,
                cluster text,
                concat_seq text,
                concat_match text,
                concat_method text
            )'''
        )
        cur.execute('''CREATE TABLE sequence (
                id INTEGER,
                untrimmed_len integer,
                trimmed_len integer,
                n_count integer,
                seq_trimmed text,
                record blob,
                FOREIGN KEY(id) REFERENCES tags(id) DEFERRABLE INITIALLY
                DEFERRED
            )'''
        )
        cur.execute("CREATE INDEX idx_sequence_cluster on tags(cluster)")
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

def insert_record_to_db(cur, tagged):
    name = tagged.read.identifier.split(' ')[0].lstrip('>')
    cur.execute('''INSERT INTO tags (
            name,
            outer,
            outer_seq,
            outer_match,
            outer_method,
            inner,
            inner_seq,
            inner_match,
            inner_method,
            cluster,
            concat_seq,
            concat_match,
            concat_method
        ) 
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)''', 
        (
            name,
            tagged.outer_name,
            tagged.outer_seq,
            tagged.outer_match,
            tagged.outer_type,
            tagged.inner_name,
            tagged.inner_seq,
            tagged.inner_match,
            tagged.inner_type,
            tagged.cluster,
            tagged.concat_seq,
            tagged.concat_match,
            tagged.concat_type
        )
    )
    key = cur.lastrowid
    # pick the actual sequence
    sequence_pickle = pickle.dumps(tagged.read,1)
    cur.execute('''INSERT INTO sequence (
            id,
            trimmed_len,
            n_count,
            seq_trimmed,
            record
        ) 
        VALUES (?,?,?,?,?)''',
        (
            key,
            len(tagged.read.sequence),
            tagged.read.sequence.lower().count('n'),
            tagged.read.sequence,
            sequence_pickle
        )
    )
