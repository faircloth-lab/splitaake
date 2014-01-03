"""
File: test_demuxipy.py
Author: Brant Faircloth

Created by Brant Faircloth on 09 October 2011 13:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: Tests for demuxipy/lib.py

"""

import os
import unittest
import ConfigParser
from demuxipy import *
from seqtools.sequence.transform import DNA_reverse_complement

import pdb

class TestParamsInstance(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

class TestListQueue(unittest.TestCase):
    def setUp(self):
        pass

    def test_put(self):
        lq = ListQueue()
        lq.put(4) 
        assert lq == [4]
        lq.put(6)
        assert lq == [4,6]

    def test_get(self):
        lq = ListQueue()
        lq.put(4)
        lq.put(6)
        assert lq.get() == 6
        assert len(lq) == 1
        assert lq.get() == 4
        assert len(lq) == 0

class TestParametersMethods(unittest.TestCase):
    def setUp(self):
        conf = ConfigParser.ConfigParser()
        conf.read('./test-data/demuxi-test.conf')
        self.p = Parameters(conf)

    def test_get_all_inner(self):
        expected = {
                'simplex2': 'GCTGCTGGCGAATC',
                'simplex3': 'CGTGCTGCGGAACT',
                'simplex1': 'CGTCGTGCGGAATC'
            }
        observed = self.p._get_all_inner()
        assert observed == expected

    def test_get_all_outer(self):
        expected = {
                'mid15': 'ATACGACGTA',
                'mid16': 'TCACGTACTA',
                'mid17': 'CGTCTAGTAC'
            }
        observed = self.p._get_all_outer()
        assert observed == expected

    def test_wrong_outer_type(self):
        self.p.outer_type = 'Bob'
        self.assertRaises(AssertionError, self.p._check_values)
        self.p.outer_type = 'Single'
        self.p._check_values()
        self.p.outer_type = 'Both'
        self.p._check_values()

    def test_wrong_inner_type(self):
        self.p.inner_type = 'Bob'
        self.assertRaises(AssertionError, self.p._check_values)
        self.p.inner_type = 'Single'
        self.p._check_values()
        self.p.inner_type = 'Both'
        self.p._check_values()

    def test_wrong_inner_orientation(self):
        self.p.inner_orientation = 'Bob'
        self.assertRaises(AssertionError, self.p._check_values)
        self.p.inner_orientation = 'Forward'
        self.p._check_values()
        self.p.inner_orientation = 'Reverse'
        self.p._check_values()

    def test_wrong_outer_orientation(self):
        self.p.outer_orientation = 'Bob'
        self.assertRaises(AssertionError, self.p._check_values)
        self.p.outer_orientation = 'Forward'
        self.p._check_values()
        self.p.outer_orientation = 'Reverse'
        self.p._check_values()

class TestSequenceTagsSimpleMethods(unittest.TestCase):
    def setUp(self):
        conf = ConfigParser.ConfigParser()
        conf.read('./test-data/demuxi-test.conf')
        self.p = Parameters(conf)
        self.inner = {
                'simplex2': 'GCTGCTGGCGAATC',
                'simplex3': 'CGTGCTGCGGAACT',
                'simplex1': 'CGTCGTGCGGAATC'
            }
        self.outer = {
                'mid15': 'ATACGACGTA',
                'mid16': 'TCACGTACTA',
                'mid17': 'CGTCTAGTAC'
            }

    def test_get_length_with_inner(self):
        observed = self.p.sequence_tags._get_length(self.inner)
        expected = len(self.inner['simplex1'])

    def test_get_length_with_outer(self):
        observed = self.p.sequence_tags._get_length(self.outer)
        expected = len(self.outer['mid15'])

    def test_build_forward_regex(self):
        t = self.outer.values()
        gap = 5
        observed = self.p.sequence_tags._build_regex(t, gap)
        for k,v in enumerate(observed):
            assert v.pattern == '^[acgtnACGTN]{{0,{0}}}{1}'.format(gap, t[k])

    def test_build_reverse_regex(self):
        t = self.outer.values()
        gap = 5
        observed = self.p.sequence_tags._build_regex(t, gap, True)
        for k,v in enumerate(observed):
            assert v.pattern == '{1}[acgtnACGTN]{{0,{0}}}$'.format(gap, t[k])

    def test_reverse_dict(self):
        expected = {
                'ATACGACGTA':'mid15',
                'TCACGTACTA':'mid16',
                'CGTCTAGTAC':'mid17'
            }
        observed = self.p.sequence_tags._reverse_dict(self.outer)
        assert expected == observed

    def test_parse_group(self):
        rows = self.p.conf.items('TestCombos')
        expected = ('simplex1','simplex2','duck')
        for k,v in enumerate(rows):
            m,l,org = self.p.sequence_tags._parse_group(v)
            assert (m,l,org) == expected

    def test_make_combinatorial_cluster_map(self):
        name = self.p.sequence_tags._make_combinatorial_cluster_map('ATACGACGTA','TCACGTACTA')
        assert name == 'ATACGACGTA,TCACGTACTA'

class TestSequenceTagsRegexAndStringMethods(unittest.TestCase):
    def setUp(self):
        conf = ConfigParser.ConfigParser()
        conf.read('./test-data/demuxi-test.conf')
        self.p = Parameters(conf)
        self.o_gap = self.p.outer_buffer
        self.i_gap = self.p.inner_buffer
        self.inner = {
                'simplex2': 'GCTGCTGGCGAATC',
                'simplex3': 'CGTGCTGCGGAACT',
                'simplex1': 'CGTCGTGCGGAATC'
            }
        self.outer = {
                'mid15': 'ATACGACGTA',
                'mid16': 'TCACGTACTA',
                'mid17': 'CGTCTAGTAC'
            }
        self.tags = set(['TCACGTACTA', 'ATACGACGTA'])


    def test_generate_outer_regex_single(self):
        self.p.sequence_tags.outers = {
                'forward_string':self.tags,
                'reverse_string':self.tags
                }
        forward = ['^[acgtnACGTN]{{0,{0}}}{1}'.format(self.o_gap, t) for t in
                self.tags]
        #pdb.set_trace()
        self.p.sequence_tags._generate_outer_regex('single')
        assert set(self.p.sequence_tags.outers.keys()) == set([
                'forward_string',
                'reverse_string',
                'forward_regex'
            ])
        for k,v in enumerate(self.p.sequence_tags.outers['forward_regex']):
            assert forward[k] == v.pattern


    def test_generate_outer_regex_both(self):
        self.p.sequence_tags.outers = {
                'forward_string':self.tags,
                'reverse_string':self.tags
                }
        forward = ['^[acgtnACGTN]{{0,{0}}}{1}'.format(self.o_gap, t) for t in
                self.tags]
        reverse = ['{1}[acgtnACGTN]{{0,{0}}}$'.format(self.o_gap, t) for t in
                self.tags]
        #pdb.set_trace()
        self.p.sequence_tags._generate_outer_regex('both')
        assert set(self.p.sequence_tags.outers.keys()) == set([
                'forward_string',
                'reverse_string',
                'forward_regex',
                'reverse_regex'
            ])
        for k,v in enumerate(self.p.sequence_tags.outers['forward_regex']):
            assert forward[k] == v.pattern
        for k,v in enumerate(self.p.sequence_tags.outers['reverse_regex']):
            assert reverse[k] == v.pattern


    def test_generate_inner_regex_single(self):
        self.p.sequence_tags.inners = {
                'mid15':{
                    'forward_string':self.tags,
                    'reverse_string':self.tags
                    }
                }
        forward = ['^[acgtnACGTN]{{0,{0}}}{1}'.format(self.o_gap, t) for t in
                self.tags]
        #pdb.set_trace()
        self.p.sequence_tags._generate_inner_regex('single')
        assert set(self.p.sequence_tags.inners['mid15'].keys()) == set([
                'forward_string',
                'reverse_string',
                'forward_regex'
            ])
        for k,v in enumerate(self.p.sequence_tags.inners['mid15']['forward_regex']):
            assert forward[k] == v.pattern

    def test_generate_inner_regex_both(self):
        self.p.sequence_tags.inners = {
                'mid15':{
                    'forward_string':self.tags,
                    'reverse_string':self.tags
                    }
                }
        forward = ['^[acgtnACGTN]{{0,{0}}}{1}'.format(self.o_gap, t) for t in
                self.tags]
        reverse = ['{1}[acgtnACGTN]{{0,{0}}}$'.format(self.o_gap, t) for t in
                self.tags]
        #pdb.set_trace()
        self.p.sequence_tags._generate_inner_regex('both')
        assert set(self.p.sequence_tags.inners['mid15'].keys()) == set([
                'forward_string',
                'reverse_string',
                'forward_regex',
                'reverse_regex'
            ])
        for k,v in enumerate(self.p.sequence_tags.inners['mid15']['forward_regex']):
            assert forward[k] == v.pattern
        for k,v in enumerate(self.p.sequence_tags.inners['mid15']['reverse_regex']):
            assert reverse[k] == v.pattern

    def test_generate_outer_reverse_string(self):
        old = self.p.sequence_tags.outers
        revtags = set([DNA_reverse_complement(t) for t in self.tags])
        self.p.sequence_tags.outers = defaultdict(set)
        for t in self.tags:
            self.p.sequence_tags._generate_outer_reverse_strings(t, 'both', 'reverse')
        assert len(self.p.sequence_tags.outers.keys()) == 1
        assert self.p.sequence_tags.outers.keys()[0] == 'reverse_string'
        assert self.p.sequence_tags.outers['reverse_string'] == revtags
        self.p.sequence_tags.outers = old

    def test_generate_outer_reverse_string_forward(self):
        old = self.p.sequence_tags.outers
        revtags = self.tags
        self.p.sequence_tags.outers = defaultdict(set)
        for t in self.tags:
            self.p.sequence_tags._generate_outer_reverse_strings(t, 'both',
            'forward')
        assert len(self.p.sequence_tags.outers.keys()) == 1
        assert self.p.sequence_tags.outers.keys()[0] == 'reverse_string'
        assert self.p.sequence_tags.outers['reverse_string'] == revtags
        self.p.sequence_tags.outers = old

    def test_generate_inner_reverse_string(self):
        old = self.p.sequence_tags.inners
        revtags = set([DNA_reverse_complement(t) for t in self.tags])
        m = 'mid15'
        self.p.sequence_tags.inners = defaultdict(lambda : defaultdict(set))
        for t in self.tags:
            self.p.sequence_tags._generate_inner_reverse_strings(m, t, 'both', 'reverse')
        #pdb.set_trace()
        assert len(self.p.sequence_tags.inners[m].keys()) == 1
        assert self.p.sequence_tags.inners[m].keys()[0] == 'reverse_string'
        assert self.p.sequence_tags.inners[m]['reverse_string'] == revtags
        self.p.sequence_tags.inners = old

    def test_generate_inner_reverse_string_forward(self):
        old = self.p.sequence_tags.inners
        revtags = self.tags
        m = 'mid15'
        self.p.sequence_tags.inners = defaultdict(lambda : defaultdict(set))
        for t in self.tags:
            self.p.sequence_tags._generate_inner_reverse_strings(m, t, 'both',
                    'forward')
        #pdb.set_trace()
        assert len(self.p.sequence_tags.inners[m].keys()) == 1
        assert self.p.sequence_tags.inners[m].keys()[0] == 'reverse_string'
        assert self.p.sequence_tags.inners[m]['reverse_string'] == revtags
        self.p.sequence_tags.inners = old

class TestSequenceTagsGroupMethods(unittest.TestCase):
    def setUp(self):
        self.conf = ConfigParser.ConfigParser()
        self.conf.read('./test-data/demuxi-test.conf')
        self.outers = set(['ATACGACGTA', 'TCACGTACTA'])
        self.inners = set(['CGTCGTGCGGAATC', 'GCTGCTGGCGAATC'])
        self.rev_outers = set([DNA_reverse_complement(i) for i in self.outers])
        self.rev_inners = set([DNA_reverse_complement(i) for i in self.inners])
        self.mp = {
                'TCACGTACTA': 
                    {
                        'CGTCGTGCGGAATC': 'pony'
                    },
                'ATACGACGTA':
                    {
                        'CGTCGTGCGGAATC': 'cat',
                        'GCTGCTGGCGAATC': 'dog'
                    }
                }
        self.tests = (('single','reverse'), ('single','forward'),
                ('both','reverse'), ('both','forward'))

    def refresh(self, obj):
        return obj._get_sequence_tags(obj._get_all_outer(), 
                obj._get_all_inner())

    def forward_regex(self, tags, g):
        return set(['^[acgtnACGTN]{{0,{0}}}{1}'.format(g, t) for t in
                tags])

    def reverse_regex(self, tags, g):
        return set(['{1}[acgtnACGTN]{{0,{0}}}$'.format(g, t) for t in
                tags])

    def regex(self, observed, expected, buff, reverse = False):
        observed_patterns = set([r.pattern for r in \
                    observed])
        if not reverse:
            expected_patterns = self.forward_regex(expected, buff)
        else:
            expected_patterns = self.reverse_regex(expected, buff)
        return observed_patterns == expected_patterns

    def test_inner_outer_outers(self):
        p = Parameters(self.conf)
        # reset group
        p.search = 'OuterInnerGroups'
        clust = {
                'TCACGTACTA': 
                    {
                        'CGTCGTGCGGAATC': 'pony'
                    },
                'ATACGACGTA':
                    {
                        'CGTCGTGCGGAATC': 'cat',
                        'GCTGCTGGCGAATC': 'dog'
                    }
                }
        self.outer_groups(p, clust)

    def test_outers(self):
        p = Parameters(self.conf)
        p.search = 'OuterGroups'
        clust = {
                'ATACGACGTA':'shrimp',
                'TCACGTACTA':'lobsters'
                }
        self.outer_groups(p, clust)

    def test_outer_map(self):
        p = Parameters(self.conf)
        p.search = 'OuterGroups'
        clust = {
                'ATACGACGTA':
                    {
                        'None':'shrimp',
                    },
                'TCACGTACTA':
                    {
                        'None':'lobsters'
                    }
                }
        self.check_map(p, clust)

    def test_inner_outer_inners(self):
        p = Parameters(self.conf)
        # reset group
        p.search = 'OuterInnerGroups'
        clust = {
                'TCACGTACTA': 
                    {
                        'CGTCGTGCGGAATC': 'pony'
                    },
                'ATACGACGTA':
                    {
                        'CGTCGTGCGGAATC': 'cat',
                        'GCTGCTGGCGAATC': 'dog'
                    }
                }
        self.inner_groups(p, clust)

    def test_inner_outer_map(self):
        p = Parameters(self.conf)
        # reset group
        p.search = 'OuterInnerGroups'
        clust = {
                'TCACGTACTA': 
                    {
                        'CGTCGTGCGGAATC': 'pony'
                    },
                'ATACGACGTA':
                    {
                        'CGTCGTGCGGAATC': 'cat',
                        'GCTGCTGGCGAATC': 'dog'
                    }
                }
        self.check_map(p, clust)

    def test_inners(self):
        p = Parameters(self.conf)
        p.search = 'InnerGroups'
        clust = {
                'None':
                    {
                        'CGTCGTGCGGAATC':'walrus',
                        'GCTGCTGGCGAATC':'whale',
                        'CGTGCTGCGGAACT':'porpoise'
                    }
                }
        self.inner_groups(p, clust)
    
    def test_inner_map(self):
        p = Parameters(self.conf)
        p.search = 'InnerGroups'
        clust = {
                'None':
                    {
                        'CGTCGTGCGGAATC':'walrus',
                        'GCTGCTGGCGAATC':'whale',
                        'CGTGCTGCGGAACT':'porpoise'
                    }
                }
        self.check_map(p, clust)

    def check_map(self, p, mp):
        #if p.search == 'HierarchicalCombinatorial':
        #   pdb.set_trace()
        p.sequence_tags = None
        p.sequence_tags = self.refresh(p)
        for k,v in p.sequence_tags.cluster_map.iteritems():
            if type(v) is defaultdict:
                for i,j in v.iteritems():
                    assert mp[k][i] == j

    def outer_groups(self, p, mp):
        for t in self.tests:
            p.outer_type, p.outer_orientation = t[0], t[1]
            p.outer_buffer = 5
            p.sequence_tags = None
            p.sequence_tags = self.refresh(p)
            #pdb.set_trace()
            if t[0] == 'single':
                assert set(p.sequence_tags.outers.keys()) == \
                    set([
                        'forward_string',
                        'forward_regex'
                    ])
                assert p.sequence_tags.outers['forward_string'] == \
                    self.outers
                assert self.regex(p.sequence_tags.outers['forward_regex'],
                        self.outers, p.outer_buffer)
            else:
                assert set(p.sequence_tags.outers.keys()) == \
                    set([
                        'forward_string',
                        'forward_regex',
                        'reverse_string',
                        'reverse_regex'
                        ])
                if t[1] == 'reverse':
                    assert p.sequence_tags.outers['forward_string'] == \
                        self.outers
                    assert self.regex(p.sequence_tags.outers['forward_regex'],
                        self.outers, p.outer_buffer)
                    assert p.sequence_tags.outers['reverse_string'] == \
                        self.rev_outers
                    assert self.regex(p.sequence_tags.outers['reverse_regex'],
                        self.rev_outers, p.outer_buffer, True)

                elif t[1] == 'forward':
                    assert p.sequence_tags.outers['forward_string'] == \
                        self.outers
                    assert self.regex(p.sequence_tags.outers['forward_regex'],
                        self.outers, p.outer_buffer)
                    assert p.sequence_tags.outers['reverse_string'] == \
                        self.outers
                    assert self.regex(p.sequence_tags.outers['reverse_regex'],
                        self.outers, p.outer_buffer, True)

    def inner_groups(self, p, mp):
        for t in self.tests:
            p.inner_type, p.inner_orientation = t[0], t[1]
            p.inner_buffer = 5
            p.sequence_tags = None
            p.sequence_tags = self.refresh(p)
            for outer in p.sequence_tags.inners:
                outer = str(outer)
                if t[0] == 'single':
                    assert set(p.sequence_tags.inners[outer].keys()) == \
                    set([
                        'forward_string',
                        'forward_regex'
                    ])
                    assert p.sequence_tags.inners[outer]['forward_string'] == \
                            set(mp[outer].keys())
                    assert self.regex(p.sequence_tags.inners[outer]['forward_regex'],
                        mp[outer].keys(), p.inner_buffer)
                else:
                    if t[1] == 'reverse':
                        assert p.sequence_tags.inners[outer]['forward_string'] == \
                            set(mp[outer].keys())
                        assert self.regex(p.sequence_tags.inners[outer]['forward_regex'],
                            mp[outer].keys(), p.inner_buffer)
                        revs = set([DNA_reverse_complement(i) for i in \
                                    mp[outer].keys()])
                        assert p.sequence_tags.inners[outer]['reverse_string'] == \
                            revs
                        assert self.regex(p.sequence_tags.inners[outer]['reverse_regex'],
                            revs, p.inner_buffer, True)

                    elif t[1] == 'forward':
                        assert p.sequence_tags.inners[outer]['forward_string'] == \
                            set(mp[outer].keys())
                        assert self.regex(p.sequence_tags.inners[outer]['forward_regex'],
                            mp[outer].keys(), p.inner_buffer)
                        revs = set(mp[outer].keys())
                        assert p.sequence_tags.inners[outer]['reverse_string'] == \
                            revs
                        assert self.regex(p.sequence_tags.inners[outer]['reverse_regex'],
                            revs, p.inner_buffer, True)

    def test_outer_combo(self):
        p = Parameters(self.conf)
        p.search = 'OuterCombinatorial'
        clust = {
                'None':
                    {
                        'ATACGACGTA,TCACGTACTA': 'cat',
                        'TCACGTACTA,ATACGACGTA': 'dog',
                    }
                }
        self.outer_combo_groups(p)

    def test_outer_combo_map(self):
        p = Parameters(self.conf)
        p.search = 'OuterCombinatorial'
        clust = {
                'None':
                    {
                        'ATACGACGTA,TCACGTACTA': 'cat',
                        'TCACGTACTA,ATACGACGTA': 'dog',
                    }
                }
        self.check_map(p, clust)

    def test_inner_combo(self):
        p = Parameters(self.conf)
        p.search = 'InnerCombinatorial'
        clust = {
                'None':
                    {
                        'CGTCGTGCGGAATC,GCTGCTGGCGAATC': 'duck',
                        'GCTGCTGGCGAATC,CGTCGTGCGGAATC': 'goose',
                    }
                }
        self.inner_combo_groups(p, clust)

    def test_inner_combo_map(self):
        p = Parameters(self.conf)
        p.search = 'InnerCombinatorial'
        clust = {
                'None':
                    {
                        'CGTCGTGCGGAATC,GCTGCTGGCGAATC': 'duck',
                        'GCTGCTGGCGAATC,CGTCGTGCGGAATC': 'goose',
                    }
                }
        self.check_map(p, clust)

    def test_hierarchical_combo_outers(self):
        p = Parameters(self.conf)
        p.search = 'HierarchicalCombinatorial'
        self.outer_combo_groups(p)

    def test_hierarchical_combo_inners(self):
        p = Parameters(self.conf)
        p.search = 'HierarchicalCombinatorial'
        clust = {
                'TCACGTACTA':
                    {
                        'CGTCGTGCGGAATC,GCTGCTGGCGAATC': 'bushbaby',
                        'GCTGCTGGCGAATC,CGTCGTGCGGAATC': 'hedgehog'
                    },
                'ATACGACGTA':
                    {
                        'GCTGCTGGCGAATC,CGTCGTGCGGAATC': 'opossum',
                        'CGTCGTGCGGAATC,GCTGCTGGCGAATC': 'newt'
                    }
                }
        self.inner_combo_groups(p, clust)


    def test_hierarchical_combo_map(self):
        p = Parameters(self.conf)
        p.search = 'HierarchicalCombinatorial'
        clust = {
                'TCACGTACTA':
                    {
                        'CGTCGTGCGGAATC,GCTGCTGGCGAATC': 'bushbaby',
                        'GCTGCTGGCGAATC,CGTCGTGCGGAATC': 'hedgehog'
                    },
                'ATACGACGTA':
                    {
                        'GCTGCTGGCGAATC,CGTCGTGCGGAATC': 'opossum',
                        'CGTCGTGCGGAATC,GCTGCTGGCGAATC': 'newt'
                    }
                }
        self.check_map(p, clust)


    def outer_combo_groups(self, p):
        for t in ['reverse', 'forward']:
            # t[0] == 'single' has no meaning in combinatorial context
            p.outer_type = 'both'
            p.outer_orientation = t
            p.outer_buffer = 5
            p.sequence_tags = None
            p.sequence_tags = self.refresh(p)
            #pdb.set_trace()
            assert set(p.sequence_tags.outers.keys()) == \
                    set([
                        'forward_string',
                        'forward_regex',
                        'reverse_string',
                        'reverse_regex'
                        ])
            if t == 'reverse':
                assert p.sequence_tags.outers['forward_string'] == \
                    self.outers
                assert self.regex(p.sequence_tags.outers['forward_regex'],
                    self.outers, p.outer_buffer)
                assert p.sequence_tags.outers['reverse_string'] == \
                    self.rev_outers
                assert self.regex(p.sequence_tags.outers['reverse_regex'],
                    self.rev_outers, p.outer_buffer, True)
            else:
                assert p.sequence_tags.outers['forward_string'] == \
                    self.outers
                assert self.regex(p.sequence_tags.outers['forward_regex'],
                    self.outers, p.outer_buffer)
                assert p.sequence_tags.outers['reverse_string'] == \
                    self.outers
                assert self.regex(p.sequence_tags.outers['reverse_regex'],
                        self.outers, p.outer_buffer, True)

    def inner_combo_groups(self, p, mp):
        for t in ['reverse', 'forward']:
            # t[0] == 'single' has no meaning in combinatorial context
            p.inner_type = 'both'
            p.inner_orientation = t
            p.inner_buffer = 5
            p.sequence_tags = None
            p.sequence_tags = self.refresh(p)
            #pdb.set_trace()
            for outer in p.sequence_tags.inners:
                assert set(p.sequence_tags.inners[outer].keys()) == \
                        set([
                            'forward_string',
                            'forward_regex',
                            'reverse_string',
                            'reverse_regex'
                            ])
                #if p.search == 'HierarchicalCombinatorial':
                #    pdb.set_trace()
                if t == 'reverse':
                    assert p.sequence_tags.inners[outer]['forward_string'] == \
                        self.inners
                    assert self.regex(p.sequence_tags.inners[outer]['forward_regex'],
                        self.inners, p.inner_buffer)
                    assert p.sequence_tags.inners[outer]['reverse_string'] == \
                        self.rev_inners
                    assert self.regex(p.sequence_tags.inners[outer]['reverse_regex'],
                        self.rev_inners, p.inner_buffer, True)
                else:
                    assert p.sequence_tags.inners[outer]['forward_string'] == \
                        self.inners
                    assert self.regex(p.sequence_tags.inners[outer]['forward_regex'],
                        self.inners, p.inner_buffer)
                    #pdb.set_trace()
                    revs = set(','.join(mp[outer].keys()).split(','))
                    assert p.sequence_tags.inners[outer]['reverse_string'] == \
                        revs
                    assert self.regex(p.sequence_tags.inners[outer]['reverse_regex'],
                        revs, p.inner_buffer, True)

class TestSequenceTagsAllPossibleTags(unittest.TestCase):
    def setUp(self):
        conf = ConfigParser.ConfigParser()
        conf.read('./test-data/demuxi-test.conf')
        self.p = Parameters(conf)

    def refresh(self, obj):
        return obj._get_sequence_tags(obj._get_all_outer(), 
                obj._get_all_inner())

    def test_1(self):
        p.sequence_tags = None
        p.sequence_tags = self.refresh(p)


'''
class TestSequenceTagsCombinatorialMethods(unittest.TestCase):
    def setUp(self):
        conf = ConfigParser.ConfigParser()
        conf.read('./test-data/demuxi-test.conf')
        self.p = Parameters(conf)

    def _combinatorial_tester(self, o, t, e):
        self.p.inner_orientation = o
        self.p.inner_type = t
        test = self.p._get_sequence_tags(self.p._get_all_outer(), 
                self.p._get_all_inner())
        assert test.cluster_map.keys() == ['None']
        assert len(test.cluster_map['None']) == len(e)
        for k,v in e.iteritems():
            assert k in test.cluster_map['None'].keys()
            assert test.cluster_map['None'][k] == v

    def test_combinatorial_cluster_single_reverse(self):
        """
        duck = CGTCGTGCGGAATC - read - GATTCGCCAGCAGC

        """
        self.p.search = 'InnerCombinatorial'
        expected = {
                'GCTGCTGGCGAATC,AGTTCCGCAGCACG': 'goose',
                'CGTCGTGCGGAATC,AGTTCCGCAGCACG': 'brant',
                'CGTCGTGCGGAATC,GATTCGCCAGCAGC': 'duck',
            }
        self._combinatorial_tester('Reverse','Single', expected)

    def test_combinatorial_cluster_single(self):
        """
        duck = CGTCGTGCGGAATC - read - GCTGCTGGCGAATC

        """
        self.p.search = 'InnerCombinatorial'
        expected = {
                'GCTGCTGGCGAATC,CGTGCTGCGGAACT': 'goose',
                'CGTCGTGCGGAATC,CGTGCTGCGGAACT': 'brant',
                'CGTCGTGCGGAATC,GCTGCTGGCGAATC': 'duck',
            }
        self._combinatorial_tester('Forward','Single', expected)

    def test_combinatorial_cluster_both_reverse(self):
        """
        duck = CGTCGTGCGGAATC - read - GATTCGCCAGCAGC
        duck = GCAGCACGCCTTAG - read - CTAAGCGGTCGTCG

        """
        self.p.search = 'InnerCombinatorial'
        expected = {
                'GCTGCTGGCGAATC,AGTTCCGCAGCACG': 'goose',
                'CGACGACCGCTTAG,TCAAGGCGTCGTGC': 'goose',
                'GCAGCACGCCTTAG,TCAAGGCGTCGTGC': 'brant',
                'CGTCGTGCGGAATC,AGTTCCGCAGCACG': 'brant',
                'CGTCGTGCGGAATC,GATTCGCCAGCAGC': 'duck',
                'GCAGCACGCCTTAG,CTAAGCGGTCGTCG': 'duck'
            }
        self._combinatorial_tester('Reverse','Both', expected)

    def test_combinatorial_cluster_both_forward(self):
        """
        duck = CGTCGTGCGGAATC - read - GCTGCTGGCGAATC
        duck = GCAGCACGCCTTAG - read - CGACGACCGCTTAG
        """
        self.p.search = 'InnerCombinatorial'
        expected = {
                'GCTGCTGGCGAATC,CGTGCTGCGGAACT': 'goose',
                'CGACGACCGCTTAG,GCACGACGCCTTGA': 'goose',
                'CGTCGTGCGGAATC,GCTGCTGGCGAATC': 'duck',
                'GCAGCACGCCTTAG,CGACGACCGCTTAG': 'duck',
                'CGTCGTGCGGAATC,CGTGCTGCGGAACT': 'brant',
                'GCAGCACGCCTTAG,GCACGACGCCTTGA': 'brant'
            }
        self._combinatorial_tester('Forward','Both', expected)

    def tearDown(self):
        del(self.p)

'''

if __name__ == '__main__':
    unittest.main()


