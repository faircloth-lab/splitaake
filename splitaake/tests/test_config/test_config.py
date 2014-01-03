"""
File: test_demuxipy.py
Author: Brant Faircloth

Created by Brant Faircloth on 09 October 2011 13:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: Tests for demuxipy/lib.py

"""

import os
import ConfigParser
import pytest

from seqtools.sequence.transform import DNA_reverse_complement

from splitaake import config

import pdb


class TestListQueue:
    def test_put(self, lq):
        lq.put(4)
        assert lq == [4]
        lq.put(6)
        assert lq == [4, 6]

    def test_get(self, lq):
        assert lq.get() == 6
        assert len(lq) == 1
        assert lq.get() == 4
        assert len(lq) == 0


class TestParallelism:
    def test_parallel(self, conffile):
        p = config.Parallelism(conffile)
        assert p.multiprocessing is False
        assert p.cores == 1

    def test_improper_config_1(self):
        with pytest.raises(config.ConfigurationError) as excinfo:
            config.Parallelism(True, None)
            assert excinfo.value == ("[Multiprocessing] You have specified "
                                     "multiprocessing but no number of cores "
                                     "[or 'auto']")

    def test_improper_config_2(self):
        with pytest.raises(config.ConfigurationError) as excinfo:
            config.Parallelism(True, 48)
            assert excinfo.value == ("[Multiprocessing] You have specified "
                                     "more compute cores than you have")


class TestDb:
    def test_db(self, conffile, cwd):
        p = config.Db(conffile)
        assert p.create is True
        assert p.name == os.path.join(
            cwd,
            'test.sqlite'
        )

    def test_improper_config_1(self):
        with pytest.raises(config.ConfigurationError) as excinfo:
            config.Db(True, "")
            assert excinfo.value == "[Database] NAME must be specified"


class TestReads:
    def test_get_reads(self, conffile, cwd):
        p = config.Reads(conffile)
        assert os.path.basename(p.r1) == "test1-r1.fastq"
        assert os.path.basename(p.r2) == "test1-r2.fastq"
        assert p.r1 == os.path.join(cwd, "data/test1-r1.fastq")
        assert p.r2 == os.path.join(cwd, "data/test1-r2.fastq")

    def test_read_failure(self):
        with pytest.raises(IOError) as excinfo:
            config.Reads('this/fails-R1.fastq', 'this/fails-R2.fastq')
            assert excinfo.value == "[Sequence] File {0} does not exist."


class TestQuality:
    def test_get_quality(self, conffile):
        p = config.Quality(conffile)
        assert p.trim is False
        assert p.min == 15
        assert p.dropn is True
        assert p.drop_len == 50

    def test_quality_failure(self):
        with pytest.raises(config.ConfigurationError) as excinfo:
            config.Quality(True, True, 100, 1)
            assert excinfo.value == ("[Quality] Min quality must be "
                                     "0 <= min <= 93.")


class TestTags:
    def test_name_d(self, tags):
        assert tags.name_d == {
            'p2,n1': 'combo5',
            'p2,n2': 'combo6',
            'p4,n2': 'combo10',
            'p3,n4': 'combo8',
            'p3,n3': 'combo7',
            'p1,n4': 'combo4',
            'p4,n1': 'combo9',
            'p1,n1': 'combo1',
            'p1,n2': 'combo2',
            'p1,n3': 'combo3'
        }

    def test_seq_d(self, tags):
        assert tags.seq_d == {
            'p2': 'CGATgg',
            'p3': 'GTCAca',
            'p1': 'ACTG',
            'p4': 'TAGCaca',
            'n1': 'CGATc',
            'n2': 'GTCAgc',
            'n3': 'TAGCagc',
            'n4': 'ACTGtagc'
        }

    def test_rev_seq_d(self, tags):
        assert tags.rev_seq_d == {
            'CGATgg': 'p2',
            'GTCAca': 'p3',
            'ACTG': 'p1',
            'TAGCaca': 'p4',
            'CGATc': 'n1',
            'GTCAgc': 'n2',
            'TAGCagc': 'n3',
            'ACTGtagc': 'n4'
        }

    def test_max_tag_length_1(self, tags):
        assert tags.max_tag_length == 8

    def test_get_sequences(self, test_tags):
        assert test_tags.seq_d == {
            'p1': 'TTTT',
            'n1': 'GGGG',
            'n2': 'AAAA'
        }
        assert test_tags.max_tag_length == 4

    def test_get_reverse_sequences(self, test_tags):
        test_tags._get_reverse_sequences()
        assert test_tags.rev_seq_d == {
            'TTTT': 'p1',
            'GGGG': 'n1',
            'AAAA': 'n2'
        }

    def test_get_r1_tag_name(self, test_tags):
        assert test_tags._get_r1_tag_names() == set(['n1', 'p1'])

    def test_get_r1_tag_d(self, test_tags):
        faux = set(['n1', 'p1'])
        assert test_tags._get_read_tag_dicts(faux) == {
            'n1': 'GGGG',
            'p1': 'TTTT'
        }

    def test_get_r2_tag_name(self, test_tags):
        assert test_tags._get_r2_tag_names() == set(['n1', 'n2', 'p1'])

    def test_get_r2_tag_d(self, test_tags):
        faux = set(['n1', 'n2', 'p1'])
        assert test_tags._get_read_tag_dicts(faux) == {
            'n1': 'GGGG',
            'n2': 'AAAA',
            'p1': 'TTTT'
        }

    def test_name_lookup(self, test_tags):
        assert test_tags.name_lookup('TTTT') == 'p1'
        assert test_tags.name_lookup('GGGG') == 'n1'
        assert test_tags.name_lookup('AAAA') == 'n2'

    def test_combo_lookup(self, test_tags):
        assert test_tags.combo_lookup('p1', 'n1') == 'combo1'
        assert test_tags.combo_lookup('p1', 'n2') == 'combo2'
        assert test_tags.combo_lookup('n1', 'p1') == 'combo3'
'''
@pytest.fixture(scope="module")
def p(conf_file):
    return config.Parameters(conf_file)

class TestParametersMethods:
    def test_first(self, p):
        if p:
            return
'''
