"""
File: lib.py
Author: Brant Faircloth

Created by Brant Faircloth on 01 October 2011 11:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description:  common files for demuxi.py

"""
import os
import re
import sys
import argparse
import ConfigParser
from collections import defaultdict
from multiprocessing import cpu_count
from seqtools.sequence.fasta import FastaSequence
from seqtools.sequence.transform import DNA_reverse_complement, DNA_complement
from seqtools.sequence.transform import reverse as DNA_reverse

import pdb


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


class ListQueue(list):
    def __init__(self):
        list.__init__(self)

    def put(self, item):
        """append an item to the list"""
        self.append(item)

    def get(self):
        """return an item from the list"""
        return self.pop()


class Parameters:
    '''linkers.py run parameters'''
    def __init__(self, conf):
        self.conf = conf
        try:
            self.fasta = os.path.abspath(os.path.expanduser(
                    self.conf.get('Sequence', 'fasta').strip("'")))
            self.quality = os.path.abspath(os.path.expanduser( \
                    self.conf.get('Sequence', 'quality').strip("'")))
            self.fastq, self.r1, self.r2 = False, False, False, False
        #pdb.set_trace()
        except ConfigParser.NoOptionError:
            self.r1 = self.conf.get('Sequence', 'r1').strip("'")
            self.r2 = self.conf.get('Sequence', 'r2').strip("'")
            self.fasta, self.quality, self.fastq = False, False, False
        except ConfigParser.NoOptionError:
            self.fastq = self.conf.get('Sequence', 'fastq').strip("'")
            self.fasta, self.quality, self.r1, self.r2 = False, False, False, False
        except ConfigParser.NoOptionError:
            print "Cannot find valid sequence files in [Sequence] section of {}".format(self.conf)
        self.db = self.conf.get('Database', 'DATABASE')
        self.qual_trim = self.conf.getboolean('Quality', 'QualTrim')
        self.min_qual = self.conf.getint('Quality', 'MinQualScore')
        self.drop = self.conf.getboolean('Quality', 'DropN')
        self.concat_check = self.conf.getboolean('Concatemers', 'ConcatemerChecking')
        self.concat_fuzzy = self.conf.getboolean('Concatemers', 'ConcatemerFuzzyMatching')
        self.concat_allowed_errors = self.conf.getboolean('Concatemers', 'ConcatemerAllowedErrors')
        self.search = self.conf.get('Search', 'SearchFor')
        #if self.search.lower() in ['innergroups', 'outerinnergroups', 'hierarchicalcombinatorial']:
        #    assert self.conf.has_section('InnerTags')
        #elif self.search == 'OuterGroups':
        #    assert self.conf.has_section('OuterTags')
        #elif self.search == 'OuterInnerGroups':
        #    assert self.conf.has_section('OuterTags')
        #    assert self.conf.has_section('InnerTags')
        if self.conf.has_section('OuterTags'):
            self.outer = self.conf.getboolean('OuterTags', 'Search')
            self.outer_type = self.conf.get('OuterTags', 'TrimType')
            self.outer_buffer = self.conf.getint('OuterTags', 'Buffer')
            self.outer_orientation = self.conf.get('OuterTags', 'ThreePrimeOrientation')
            self.outer_trim = self.conf.getint('OuterTags', 'Trim')
            self.outer_fuzzy = self.conf.getboolean('OuterTags', 'FuzzyMatching')
            self.outer_errors = self.conf.getint('OuterTags', 'AllowedErrors')
        if self.conf.has_section('InnerTags'):
            self.inner = self.conf.getboolean('InnerTags', 'Search')
            self.inner_type = self.conf.get('InnerTags', 'TrimType')
            self.inner_buffer = self.conf.getint('InnerTags', 'Buffer')
            self.inner_orientation = self.conf.get('InnerTags', 'ThreePrimeOrientation')
            self.inner_trim = self.conf.getint('InnerTags', 'Trim')
            self.inner_fuzzy = self.conf.getboolean('InnerTags', 'FuzzyMatching')
            self.inner_errors = self.conf.getint('InnerTags', 'AllowedErrors')
        #all_outer             = self._get_all_outer()
        #all_inner             = self._get_all_inner()
        self._check_values()
        self.sequence_tags = self._get_sequence_tags(self._get_all_outer(),
                self._get_all_inner())
        if self.conf.has_section('Primers'):
            self.primers = Primers(
                    self.conf.get('Primers', 'Forward'),
                    self.conf.get('Primers', 'Reverse'),
                    self.conf.getint('Primers', 'Buffer'),
                    self.conf.getboolean('Primers', 'FuzzyMatching'),
                    self.conf.getint('Primers', 'AllowedErrors')
                )
        self.multiprocessing = conf.getboolean('Multiprocessing', 'Multiprocessing')
        # compute # cores for computation; leave 1 for db and 1 for sys
        #pdb.set_trace()
        if self.multiprocessing == True:
            if conf.get('Multiprocessing', 'processors').lower() == 'auto' and cpu_count > 2:
                self.num_procs = cpu_count() - 1
            elif conf.get('Multiprocessing', 'processors').lower() != 'auto' and \
                    cpu_count >= conf.getint('Multiprocessing', 'processors'):
                self.num_procs = conf.getint('Multiprocessing', 'processors')
        else:
            self.num_procs = 1

    def __str__(self):
        return "{0}({1})".format(self.__class__, self.__dict__)

    def __repr__(self):
        return "<{0} instance at {1}>".format(self.__class__, hex(id(self)))

    def _get_all_outer(self):
        # if only linkers, you don't need MIDs
        if self.search.lower() in ['outergroups', 'outerinnergroups',
                'outercombinatorial', 'hierarchicalcombinatorial']:
            return dict(self.conf.items('OuterTagSequences'))
        else:
            return None

    def _get_all_inner(self):
        # if only linkers, you don't need MIDs
        if self.search.lower() in ['innergroups', 'outerinnergroups',
                'innercombinatorial', 'hierarchicalcombinatorial']:
            return dict(self.conf.items('InnerTagSequences'))
        else:
            return None

    def _get_sequence_tags(self, all_outer, all_inner):
        return SequenceTags(
                all_outer,
                all_inner,
                self.search,
                self.conf.items(self.search),
                self.outer_buffer,
                self.inner_buffer,
                self.concat_check,
                self.outer_type,
                self.outer_orientation,
                self.inner_type,
                self.inner_orientation
            )

    def _check_values(self):
        assert self.outer_type.lower() in ['single', 'both'], \
                "Outer type must be one of ['Single','Both']"
        assert self.outer_orientation.lower() in ['reverse', 'forward'], \
                "Innert Orientation must be one of ['Forward','Reverse']"
        assert self.inner_type.lower() in ['single', 'both'], \
                "Outer type must be one of ['Single','Both']"
        assert self.inner_orientation.lower() in ['reverse', 'forward'], \
                "Inner orientation must be one of ['Forward','Reverse']"
        assert self.search.lower() in \
                [
                    'innergroups',
                    'outergroups',
                    'outerinnergroups',
                    'outercombinatorial',
                    'innercombinatorial',
                    'hierarchicalcombinatorial'
                ], \
                "SearchFor must be one of ['InnerGroups','OuterGroups'," +\
                "'OuterInnerGroups']"


class Primers:
    def __init__(self, f, r, buff, fuzzy, errors):
        self.f, self.r = {}, {}
        self.buff = buff
        self.fuzzy = fuzzy
        self.errors = errors
        self.f['string'] = [f]
        self.f['len'] = len(f)
        self.r['string'] = [r]
        self.r['len'] = len(r)
        self.f['regex'] = [re.compile('^[acgtnACGTN]{{0,{}}}{}'.format(self.buff, f))]
        self.r['regex'] = [re.compile('^[acgtnACGTN]{{0,{}}}{}'.format(self.buff, r))]


class SequenceTags:
    """ """
    def __init__(self, all_outers, all_inners, search, group, outer_gap, inner_gap,
            concat, o_type, o_orientation, i_type, i_orientation):
        self.outers = None
        self.inners = None
        self.cluster_map = None
        self.all_tags = None
        self.outer_gap = outer_gap
        self.inner_gap = inner_gap
        if all_outers:
            # map mid sequences to names
            self.reverse_outer_lookup = self._reverse_dict(all_outers)
            self.outer_len = self._get_length(all_outers)
        if all_inners:
            # map linker sequences to names
            self.reverse_inner_lookup = self._reverse_dict(all_inners)
            self.inner_len = self._get_length(all_inners)
        # pare down the list of linkers and MIDS to those we've used
        self._generate_clusters_and_get_cluster_tags(all_outers, all_inners, search,
                group, o_type, o_orientation, i_type, i_orientation)
        # do we check for concatemers?
        if concat:
            self._all_possible_tags(search)

    def __str__(self):
        return "{0}({1})".format(self.__class__, self.__dict__)

    def __repr__(self):
        return "<{0} instance at {1}>".format(self.__class__, hex(id(self)))

    def _get_length(self, tags):
        #linkers         = dict(self.conf.items('Linker'))
        lset = set([len(l) for l in tags.values()])
        assert len(lset) == 1, "Your {} sequences are difference lengths".format(name)
        return lset.pop()

    def _build_regex(self, tags, gap, rev=False):
        if not rev:
            return [re.compile('^[acgtnACGTN]{{0,{}}}{}'.format(gap, seq)) for seq in
                    tags]
        else:
            return [re.compile('{}[acgtnACGTN]{{0,{}}}$'.format(seq, gap))
                    for seq in tags]

    def _reverse_dict(self, d):
        return {v:k for k, v in d.iteritems()}

    def _parse_group(self, row, length=False):
        if length == 'short':
            m = row[0].replace(' ', '')
            return m, row[1]
        elif length == 'long':
            o1, i1, i2, o2 = row[0].replace(' ', '').split(',')
            return o1, i1, i2, o2, row[1]
        else:
            m, l = row[0].replace(' ', '').split(',')
            return m, l, row[1]

    def _generate_outer_regex(self, outer_type):
        # compile regular expressions:
        self.outers['forward_regex'] = \
                self._build_regex(self.outers['forward_string'],
                self.outer_gap)
        if outer_type.lower() == 'both':
            self.outers['reverse_regex'] = \
                    self._build_regex(self.outers['reverse_string'],
                    self.outer_gap, rev=True)

    def _generate_inner_regex(self, inner_type):
        for m in self.inners:
            self.inners[m]['forward_regex'] = \
                self._build_regex(self.inners[m]['forward_string'],
                self.inner_gap)
            if inner_type.lower() == 'both':
                self.inners[m]['reverse_regex'] = \
                    self._build_regex(self.inners[m]['reverse_string'],
                    self.inner_gap, rev=True)

    def _generate_inner_combi_regex(self, inner_type):
        for m in self.inners:
            self.inners[m]['forward_regex'] = \
                self._build_regex(self.inners[m]['forward_string'],
                self.inner_gap)
            if inner_type.lower() == 'both':
                self.inners[m]['reverse_regex'] = \
                self._build_regex(self.inners[m]['reverse_string'],
                self.inner_gap)

    def _generate_outer_reverse_strings(self, m, outer_type, outer_orientation):
        if outer_type.lower() == 'both':
            if outer_orientation.lower() == 'reverse':
                self.outers['reverse_string'].add(DNA_reverse_complement(m))
            else:
                self.outers['reverse_string'].add(m)

    def _generate_inner_reverse_strings(self, m, l, inner_type,
            inner_orientation):
        #if inner_type == 'single':
        #    pdb.set_trace()
        if inner_type.lower() == 'both':
            if inner_orientation.lower() == 'reverse':
                self.inners[str(m)]['reverse_string'].add(DNA_reverse_complement(l))
            else:
                # reverse orientation really == forward orientation
                self.inners[str(m)]['reverse_string'].add(l)

    def _generate_outer_inner_groups(self, all_outers, all_inners, group,
            o_type, o_orientation, i_type, i_orientation):
        self.outers = defaultdict(set)
        self.inners = defaultdict(lambda: defaultdict(set))
        for row in group:
            m, l, org = self._parse_group(row)
            self.outers['forward_string'].add(all_outers[m])
            self._generate_outer_reverse_strings(all_outers[m],
                    o_type, o_orientation)
            self.inners[all_outers[m]]['forward_string'].add(all_inners[l])
            self._generate_inner_reverse_strings(all_outers[m],
                    all_inners[l], i_type, i_orientation)
            self.cluster_map[all_outers[m]][all_inners[l]] = org
        # compile regular expressions for outers and inners
        self._generate_outer_regex(o_type)
        self._generate_inner_regex(i_type)

    def _generate_outer_groups(self, all_outers, group, o_type, o_orientation):
        self.outers = defaultdict(set)
        for row in group:
            m, org = self._parse_group(row, length='short')
            self.outers['forward_string'].add(all_outers[m])
            self._generate_outer_reverse_strings(all_outers[m], o_type,
                    o_orientation)
            self.cluster_map[all_outers[m]]['None'] = org
        # compile regular expressions for outers
        self._generate_outer_regex(o_type)

    def _generate_inner_groups(self, all_inners, group, i_type, i_orientation):
        self.inners = defaultdict(lambda: defaultdict(set))
        for row in group:
            l, org = self._parse_group(row, length='short')
            self.inners['None']['forward_string'].add(all_inners[l])
            self._generate_inner_reverse_strings(None, all_inners[l],
                    i_type, i_orientation)
            self.cluster_map['None'][all_inners[l]] = org
        # compile regular expressions for inners
        self._generate_inner_regex(i_type)

    def _make_combinatorial_cluster_map(self, ll, rl):
        return "{},{}".format(ll, rl)

    def _generate_inner_combinatorial_groups(self, all_inners, group,
            i_type, i_orientation):
        self.inners = defaultdict(lambda: defaultdict(set))
        for row in group:
            m, l, org = self._parse_group(row)
            self.inners['None']['forward_string'].add(all_inners[m])
            #self._generate_inner_reverse_strings(None, all_inners[l],
            #        i_type, i_orientation)
            self.inners['None']['reverse_string'].add(all_inners[l])
            name = self._make_combinatorial_cluster_map(all_inners[m], all_inners[l])
            self.cluster_map['None'][name] = org
        # compile regular expressions for combinatorial inners
        #self._generate_inner_regex(i_type)
        self._generate_inner_combi_regex(i_type)

    def _generate_outer_combinatorial_groups(self, all_outers, group, o_type,
            o_orientation):
        self.outers = defaultdict(set)
        for row in group:
            m, l, org = self._parse_group(row)
            #pdb.set_trace()
            self.outers['forward_string'].add(all_outers[m])
            self._generate_outer_reverse_strings(all_outers[l], o_type,
                    o_orientation)
            name = self._make_combinatorial_cluster_map(all_outers[m], all_outers[l])
            self.cluster_map['None'][name] = org
        # compile regular expressions for combinatorial outers
        self._generate_outer_regex(o_type)

    def _generate_hierarchical_combinatorial_groups(self, all_outers, all_inners,
            group, o_type, o_orientation, i_type, i_orientation):
        self.outers = defaultdict(set)
        self.inners = defaultdict(lambda: defaultdict(set))
        for row in group:
            m, ll, rl, n, org = self._parse_group(row, length='long')
            #pdb.set_trace()
            assert m == n, "Outer tags in hierarchical combo must agree"
            self.outers['forward_string'].add(all_outers[m])
            self._generate_outer_reverse_strings(all_outers[m], o_type,
                    o_orientation)
            self.inners[all_outers[m]]['forward_string'].add(all_inners[ll])
            self._generate_inner_reverse_strings(all_outers[m], all_inners[rl],
                    i_type, i_orientation)
            name = self._make_combinatorial_cluster_map(all_inners[ll],
                    all_inners[rl])
            self.cluster_map[all_outers[m]][name] = org
        # compile regular expressions for outers and inners
        self._generate_outer_regex(o_type)
        self._generate_inner_regex(i_type)

    def _generate_clusters_and_get_cluster_tags(self, all_outers, all_inners, search,
            group, o_type, o_orientation, i_type, i_orientation):

        self.cluster_map = defaultdict(lambda: defaultdict(str))

        if search == 'OuterGroups':
            self._generate_outer_groups(all_outers, group, o_type,
                    o_orientation)

        elif search == 'InnerGroups':
            self._generate_inner_groups(all_inners, group, i_type,
                    i_orientation)

        elif search == 'OuterInnerGroups':
            self._generate_outer_inner_groups(all_outers, all_inners, group,
                    o_type, o_orientation, i_type, i_orientation)

        elif search == 'InnerCombinatorial':
            self._generate_inner_combinatorial_groups(all_inners, group,
                    i_type, i_orientation)

        elif search == 'OuterCombinatorial':
            self._generate_outer_combinatorial_groups(all_outers, group,
                    o_type, o_orientation)

        elif search == 'HierarchicalCombinatorial':
            self._generate_hierarchical_combinatorial_groups(all_outers,
                     all_inners, group, o_type, o_orientation, i_type,
                     i_orientation)

    def _all_possible_tags(self, search):
        '''Create regular expressions for the forward and reverse complements
        of all of the tags sequences used in a run'''
        # at = all tags; rat = reverse complement all tags
        self.all_tags = defaultdict(lambda: defaultdict(set))
        if self.inners:
            for m in self.inners:
                m = str(m)
                if 'reverse_string' in self.inners[m].keys():
                    self.all_tags[m]['string'] = \
                            self.inners[m]['forward_string'].union(self.inners[m]['reverse_string'])
                else:
                    self.all_tags[m]['string'] = \
                             self.all_tags[m]['string'].union(self.inners[m]['forward_string'])
        else:
            m = 'None'
            if 'reverse_string' in self.outers.keys():
                self.all_tags[m]['string'] = \
                        self.outers['forward_string'].union(self.outers['reverse_string'])
            else:
                self.all_tags[m]['string'] = \
                             self.all_tags[m]['string'].union(self.outers['forward_string'])
        # compile
        self.all_tags[m]['regex'] = [re.compile(t) for t in
                self.all_tags[m]['string']]


class Tagged:
    '''Trimming, tag, and sequence data for individual reads'''
    def __init__(self, sequence):
        # super(Params, self).__init__()
        #assert isinstance(sequence, FastaSequence), \
        #    'The Record class must be instantiated with a FastaSequence object'
        # a biopython sequence object
        self.read = sequence
        self.outer = None
        self.outer_name = None
        self.outer_seq = None
        self.outer_match = None
        self.outer_type = None
        self.inner_name = None
        self.inner_seq = None
        self.inner_match = None
        self.inner_type = None
        self.cluster = None
        self.concat_seq = None
        self.concat_type = None
        self.concat_match = None

    #def __repr__(self):
    #    return '''<linkers.record for %s>''' % self.identifier


def reverse(items, null=False):
    '''build a reverse dictionary from a list of tuples'''
    l = []
    if null:
        items += ((None, None),)
    for i in items:
        t = (i[1], i[0])
        l.append(t)
    return dict(l)
