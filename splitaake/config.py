"""
File: lib.py
Author: Brant Faircloth

Created by Brant Faircloth on 01 October 2011 11:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description:  config parameter and tag combination classes

"""

import os
import re
import sys
import argparse
import ConfigParser
from collections import defaultdict
from multiprocessing import cpu_count
from seqtools.sequence.fasta import FastaSequence
from seqtools.sequence.fastq import FastqWriter
from seqtools.sequence.transform import DNA_reverse_complement, DNA_complement
from seqtools.sequence.transform import reverse as DNA_reverse

import pdb


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


class ConfigurationError(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)


class ListQueue(list):
    def __init__(self):
        list.__init__(self)

    def put(self, item):
        """append an item to the list"""
        self.append(item)

    def get(self):
        """return an item from the list"""
        return self.pop()


class Parallelism:
    def __init__(self, *args):
        self.multiprocessing = None
        self.cores = None
        if isinstance(args[0], ConfigParser.ConfigParser):
            self.get_conf_values(args[0])
            self.check_values()
        else:
            self.get_manual_values(args)
            self.check_values()

    def get_conf_values(self, conf):
        self.multiprocessing = conf.getboolean(
            'Multiprocessing',
            'multiprocessing'
        )

        if self.multiprocessing:
            self.cores = conf.get(
                'Multiprocessing',
                'cores'
            )
        else:
            self.cores = 1

    def get_manual_values(self, args):
        self.multiprocessing = args[0]
        self.cores = args[1]

    def check_values(self):
        if isinstance(self.cores, str):
            if self.cores == 'auto':
                self.cores = cpu_count() - 1
            elif self.cores.isdigit():
                self.cores = int(self.cores)
        elif isinstance(self.cores, int):
            self.cores = self.cores
        if self.multiprocessing and self.cores is None:
            raise ConfigurationError("[Multiprocessing] You have specified "
                                     "multiprocessing but no number of cores "
                                     "[or 'auto']")
        try:
            assert self.cores <= cpu_count()
        except:
            raise ConfigurationError("[Multiprocessing] You have specified "
                                     "more compute cores than you have")


class Db:
    def __init__(self, *args):
        self.create = None
        self.name = None
        if isinstance(args[0], ConfigParser.ConfigParser):
            self.get_conf_values(args[0])
            self.check_values()
        else:
            self.get_manual_values(args)
            self.check_values()

    def get_conf_values(self, conf):
        self.create = conf.getboolean(
            'Database',
            'database'
        )
        if self.create:
            self.name = conf.get('Database', 'name').strip("'")

    def get_manual_values(self, args):
        self.create = args[0]
        self.name = args[1]

    def check_values(self):
        if self.name is None or self.name == '':
            raise ConfigurationError("[Database] NAME must be specified")
        else:
            self.name = os.path.abspath(os.path.expanduser(self.name))


class Reads:
    def __init__(self, *args):
        self.r1 = None
        self.r2 = None
        if isinstance(args[0], ConfigParser.ConfigParser):
            self.get_conf_r1_and_r2(args[0])
        else:
            self.get_manual_r1_and_r2(args)
            self.check_values()

    def get_conf_r1_and_r2(self, conf):
        sequence_data = [
            conf.get('Sequence', 'r1'),
            conf.get('Sequence', 'r2')
        ]
        for pos, file in enumerate(sequence_data):
            self.check_and_get_path(pos, file)

    def get_manual_r1_and_r2(self, args):
        sequence_data = [
            args[0],
            args[1]
        ]
        for pos, file in enumerate(sequence_data):
            self.check_and_get_path(pos, file)

    def check_and_get_path(self, pos, file):
        fullpath = os.path.abspath(os.path.expanduser(file.strip("'")))
        try:
            assert os.path.isfile(fullpath)
        except:
            raise IOError("[Sequence] File {0} does not exist.".format(file))
        if pos == 0:
            self.r1 = fullpath
        else:
            self.r2 = fullpath


class Quality:
    def __init__(self, *args):
        self.trim = None
        self.min = None
        self.drop_n = None
        self.drop_len = None
        if isinstance(args[0], ConfigParser.ConfigParser):
            self.get_conf_quality_values(args[0])
            self.check_values()
        else:
            self.get_manual_quality_values(args)
            self.check_values()

    def get_conf_quality_values(self, conf):
        self.trim = conf.getboolean('QualitySetup', 'trim')
        self.dropn = conf.getboolean('QualitySetup', 'dropn')
        self.min = conf.getint('QualitySetup', 'min')
        try:
            self.drop_len = conf.getint('QualitySetup', 'droplen')
        except:
            self.drop_len = 1

    def get_manual_quality_values(self, args):
        self.trim = args[0]
        self.dropn = args[1]
        self.min = args[2]
        self.drop_len = args[3]

    def check_values(self):
        try:
            assert 0 <= self.min <= 93
        except:
            raise ConfigurationError("[Quality] Min quality must be "
                                     "0 <= min <= 93.")


class Tags:
    def __init__(self, *args):
        self.type = None
        self.five_p_buffer = None
        self.three_p_buffer = None
        self.three_p_orient = None
        self.fuzzy = None
        self.errors = None
        if isinstance(args[0], ConfigParser.ConfigParser):
            self.get_tag_values(args[0])
            combos = args[0].items('Combinations')
            tags = args[0].items('TagSequences')
            self.check_values()
            self.build(combos, tags)
        else:
            self.get_manual_tag_values(args)
            self.check_values()

    def build(self, combos, tags):
        self._get_combinations(combos)
        self._get_sequences(tags)
        self._get_reverse_sequences()
        self._get_r1_tags()
        self._get_r2_tags()
        #self._get_five_p_overrun(conf)
        #self._get_three_p_overrun(conf)

    def get_tag_values(self, conf):
        self.fuzzy = conf.getboolean('TagSetup', 'FuzzyMatching')
        if self.fuzzy:
            self.errors = conf.getint('TagSetup', 'AllowedErrors')
        else:
            self.errors = 0
        self.three_p_orient = conf.get(
            'TagSetup',
            'ThreePrimeOrientation').lower()

    def get_manual_tag_values(self, args):
        self.fuzzy = args[0]
        self.errors = args[1]
        self.three_p_orient = args[2]

    def check_values(self):
        try:
            assert self.three_p_orient in ['forward', 'reverse']
        except:
            raise ConfigurationError("[TagSetup] ThreePrimeOrientation must "
                                     "be 'Forward' or 'Reverse'")

    def _get_combinations(self, combos):
        self.name_d = dict(combos)

    def _get_sequences(self, tags):
        self.seq_d = dict(tags)
        self.max_tag_length = max([len(v) for k, v in self.seq_d.iteritems()])

    def _get_reverse_sequences(self):
        self.rev_seq_d = {v: k for k, v in self.seq_d.iteritems()}

    def _get_r1_tag_names(self):
        return set([k.split(',')[0] for k in self.name_d.keys()])

    def _get_r2_tag_names(self):
        return set([k.split(',')[1] for k in self.name_d.keys()])

    def _get_read_tag_dicts(self, tag_names):
        return {k: v for k, v in self.seq_d.iteritems() if k in tag_names}

    def _get_r2_tag_d(self, r2_tag_names):
        return {k: v for k, v in self.seq_d.iteritems() if k in r2_tag_names}

    def _get_r1_tags(self):
        # make sure we reduce set to only uniques
        r1_tag_names = self._get_r1_tag_names()
        r1_tag_d = self._get_read_tag_dicts(r1_tag_names)
        #pdb.set_trace()
        self.r1 = [TagMeta(k, v, 'TagSetup') for k, v in r1_tag_d.iteritems()]

    def _get_r2_tags(self):
        # make sure we reduce set to only uniques
        r2_tag_names = self._get_r2_tag_names()
        r2_tag_d = self._get_read_tag_dicts(r2_tag_names)
        self.r2 = [TagMeta(k, v, 'TagSetup') for k, v in r2_tag_d.iteritems()]

    def name_lookup(self, tag):
        return self.rev_seq_d[tag]

    def combo_lookup(self, r1tag, r2tag):
        return self.name_d["{0},{1}".format(r1tag, r2tag)]


class TagMeta:
    def __init__(self, name, string, section):
        self._check_bases(string)
        self.name = name
        self.string = string
        self._get_tag_parts()
        self.string_len = len(self.string)
        self.tag_len = len(self.tag)
        wild_five = ''.join(['{}?'.format(i) for i in self.five_p])
        wild_three = ''.join(['{}?'.format(i) for i in self.three_p])
        self.regex = re.compile(
            '^{0}({1}){2}'.format(wild_five, self.tag, wild_three),
            re.IGNORECASE
        )
        if self.five_p != '':
            temp_format = '^{0}+'.format('*'.join(list(self.five_p)))
            self.five_p_start = re.compile(temp_format, re.IGNORECASE)
        else:
            self.five_p_start = None
        #pdb.set_trace()

    def _check_bases(self, string):
        nucleotides = set(list('ACGTacgt'))
        for base in string:
            try:
                assert base in nucleotides
            except:
                raise ValueError("[{0}] Foward/Reverse bases must be in the "
                                 "alphabet [ACGTacgt]".format(section))

    def _get_tag_parts(self):
        regex = re.compile('^([acgt]*)([ACGT]+)([acgt]*)$')
        result = regex.search(self.string)
        try:
            self.five_p, self.tag, self.three_p = result.groups()
        except:
            pdb.set_trace()


class ConfOverruns(dict):
    def __init__(self, conf, tags, site):
        dict.__init__(self)
        self.trim = conf.getboolean('TagSetup', 'TrimOverruns')
        for combo in tags.name_d.keys():
            r1name, r2name = combo.split(',')
            r1t = tags.seq_d[r1name]
            r2t = tags.seq_d[r2name]
            both = {'r1':[OverrunsMeta(r2t, i.string) for i in site.r2], 'r2':[OverrunsMeta(r1t, i.string) for i in site.r1]}
            dict.__setitem__(self, combo, both)
        #pdb.set_trace()

    def __getitem__(self, key):
        val = dict.__getitem__(self, key)
        return val


class OverrunsMeta:
    def __init__(self, tag, site):
        #nucleotides = set(list('ACGTacgt'))
        #for base in string:
        #    assert base in nucleotides, ValueError("[{0}] Foward/Reverse bases must be in the alphabet [ACGTacgt]".format(section))
        #self.name = name
        self.string = DNA_reverse_complement("{0}{1}".format(tag, site)).upper()
        self.string_len = len(self.string)
        self.regex = re.compile("({0})".format(self.string), re.IGNORECASE)


class ConfSite:
    """get site parameters"""
    def __init__(self, conf):
        self.r1 = None
        self.r2 = None
        self.fuzzy = None
        self.errors = None
        self.iupac = {
            'R': ('A', 'G'),
            'Y': ('C', 'T'),
            'S': ('G', 'C'),
            'W': ('A', 'T'),
            'K': ('G', 'T'),
            'M': ('A', 'C'),
            'B': ('C', 'G', 'T'),
            'D': ('A', 'G', 'T'),
            'H': ('A', 'C', 'T'),
            'V': ('A', 'C', 'G'),
        }
        self._get_tag_values(conf)
        self._get_sequences(conf)
        self._get_forward(conf)
        self._get_reverse(conf)

    def _get_tag_values(self, conf):
        self.fuzzy = conf.getboolean('SiteSetup', 'FuzzyMatching')
        if self.fuzzy:
            self.errors = conf.getint('SiteSetup', 'AllowedErrors')
        else:
            self.errors = 0

    def _get_sequences(self, conf):
        self.max_site_length = max([len(site) for site in [conf.get('SiteSetup', 'forward'), conf.get('SiteSetup', 'reverse')]])

    def _get_forward(self, conf):
        forward = conf.get('SiteSetup', 'forward')
        forwards = self._convert_degenerate_primers(forward)
        #pdb.set_trace()
        self.r1 = [TagMeta('r1site', s, 'SiteSetup') for s in forwards]

    def _get_reverse(self, conf):
        reverse = conf.get('SiteSetup', 'reverse')
        reverses = self._convert_degenerate_primers(reverse)
        self.r2 = [TagMeta('r2site', s, 'SiteSetup') for s in reverses]

    def _convert_degenerate_primers(self, sites, lst=False):
        """Convert a given degenerate primer sequence to standard IUPAC nuclteotides"""
        if not lst:
            sites = [sites]
        new_sites = []
        still_degen = False
        for site in sites:
            for degen, bases in self.iupac.iteritems():
                if degen in site:
                    still_degen = True
                    for base in bases:
                        new_sites.append(site.replace(degen, base))
            if not still_degen:
                new_sites.append(site)
        #pdb.set_trace()
        if still_degen:
            new_sites = self._convert_degenerate_primers(new_sites, lst=True)
        new_sites = list(set(new_sites))
        #pdb.set_trace()
        return new_sites


class ConfStorage:
    def __init__(self, conf, tags):
        self.output = defaultdict(lambda: defaultdict(str))
        tld = 'splitaake-output'
        os.mkdir(tld)
        for indiv in tags.name_d.values():
            fpath = os.path.join(tld, indiv)
            os.mkdir(fpath)
            for read in ['R1', 'R2']:
                name = "{0}-{1}.fastq".format(indiv, read)
                fname = os.path.join(fpath, name)
                self.output[indiv][read] = FastqWriter(fname)

    def close(self):
        for name, files in self.output.iteritems():
            for file in files.values():
                file.close()


class Parameters:
    '''linkers.py run parameters'''
    def __init__(self, conf):
        self.db = ConfDb(conf)
        self.parallelism = ConfParallelism(conf)
        self.reads = ConfReads(conf)
        self.quality = ConfQuality(conf)
        self.tags = ConfTag(conf)
        self.site = ConfSite(conf)
        self.overruns = ConfOverruns(conf, self.tags, self.site)
        self.storage = ConfStorage(conf, self.tags)


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
