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


class MockMultiprocessing:
    def __init__(self):
        pass

    def get(self, name, value):
        if name == 'Multiprocessing':
            return str(100)

    def getboolean(self, name, value):
        return True


class TestConfParallelism:
    def test_parallel(self, parallelism):
        assert parallelism.multiprocessing is False
        assert parallelism.cores == 1

    def test_bad_parallel_values(self):
        with pytest.raises(ValueError):
            mp = MockMultiprocessing()
            config.ConfParallelism(mp)


class TestConfDb:
    def test_db(self, db):
        assert db.create is True
        assert os.path.basename(db.name) == "test.sqlite"


class TestConfReads:
    def test_get_reads(self, reads):
        assert os.path.basename(reads.r1) == "test1-r1.fastq"
        assert os.path.basename(reads.r2) == "test1-r2.fastq"

    def test_r1_failure(self, reads):
        with pytest.raises(IOError):
            reads._check_and_get_path(0, 'test100-r1.fastq')

    def test_r2_failure(self, reads):
        with pytest.raises(IOError):
            reads._check_and_get_path(1, 'test100-r2.fastq')


class TestConfQuality:
    def test_get_quality(self, quality):
        assert quality.trim is False
        assert quality.min == 15
        assert quality.dropn is True
        assert quality.drop_len == 50



'''
@pytest.fixture(scope="module")
def p(conf_file):
    return config.Parameters(conf_file)

class TestParametersMethods:
    def test_first(self, p):
        if p:
            return
'''
