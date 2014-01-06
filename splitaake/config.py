"""
File: lib.py
Author: Brant Faircloth

Created by Brant Faircloth on 01 October 2011 11:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description:  config parameter and tag combination classes

"""

import os
import re
import gzip
import argparse
import ConfigParser
from collections import defaultdict
from multiprocessing import cpu_count

from seqtools.sequence.fasta import FastaSequence
from seqtools.sequence.fastq import FastqWriter

from splitaake.core import reverse_complement

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
    def __init__(self, *args, **kwargs):
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
            if "auto" in kwargs and kwargs["auto"] is True:
                self.build(args[3], args[4])

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
        self.trim_overruns = conf.getboolean(
            'TagSetup',
            'TrimOverruns'
        )

    def get_manual_tag_values(self, args):
        self.fuzzy = args[0]
        self.errors = args[1]
        self.three_p_orient = args[2]
        self.trim_overruns = args[3]

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
    def __init__(self, name, string, section, auto=True):
        self._check_bases(string, section)
        self.name = name
        self.string = string
        self.string_len = len(self.string)
        # extract the components of a given tag - preceding
        # lowercase spacer middle uppercase tag or trailing
        # lowercase space
        if auto:
            self._get_tag_parts(section)
            self._get_tag_wildcards()
            #pdb.set_trace()

    def _check_bases(self, string, section):
        nucleotides = set(list('ACGTacgt'))
        for base in string:
            try:
                assert base in nucleotides
            except:
                raise ValueError("[{0}] Foward/Reverse bases must be in the "
                                 "alphabet [ACGTacgt]".format(section))

    def _get_tag_parts(self, section):
        result = re.search('^([acgt]*)([ACGT]+)([acgt]*)$', self.string)
        try:
            self.five_p, self.tag, self.three_p = result.groups()
        except:
            raise ConfigurationError("[{}] The tags are not structured "
                                     "correctly".format(section))
        self.tag_len = len(self.tag)

    def _get_tag_wildcards(self):
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


class Sites:
    """get site parameters"""
    def __init__(self, *args, **kwargs):
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
            'N': ('A', 'C', 'G', 'T'),
        }
        if isinstance(args[0], ConfigParser.ConfigParser):
            self.get_tag_values(args[0])
            self.check_values()
            self.build()
        else:
            self.get_manual_tag_values(args)
            self.check_values()
            if "auto" in kwargs and kwargs["auto"] is True:
                self.build()

    def get_tag_values(self, conf):
        self.fuzzy = conf.getboolean('SiteSetup', 'FuzzyMatching')
        if self.fuzzy:
            self.errors = conf.getint('SiteSetup', 'AllowedErrors')
        else:
            self.errors = 0
        self.orig_forward_site = conf.get('SiteSetup', 'forward')
        self.orig_reverse_site = conf.get('SiteSetup', 'reverse')

    def get_manual_tag_values(self, args):
        self.fuzzy = args[0]
        self.errors = args[1]
        self.orig_forward_site = args[2]
        self.orig_reverse_site = args[3]

    def check_values(self):
        try:
            assert self.errors <= 5
        except:
            raise ConfigurationError("[Sites] You almost certainly want to "
                                     "keep the number of allowable errors < 5")

    def build(self):
        self._get_sequences()
        self._get_forward()
        self._get_reverse()

    def _get_sequences(self):
        sites = [self.orig_forward_site, self.orig_reverse_site]
        self.max_site_length = max([len(site) for site in sites])

    def _get_forward(self):
        forwards = self._convert_degenerate_primers(self.orig_forward_site)
        self.r1 = [TagMeta('r1site', s, 'SiteSetup') for s in forwards]

    def _get_reverse(self):
        reverses = self._convert_degenerate_primers(self.orig_reverse_site)
        self.r2 = [TagMeta('r2site', s, 'SiteSetup') for s in reverses]

    def _convert_degenerate_primers(self, sites, lst=False):
        """
        Convert a given degenerate primer sequence to standard
        IUPAC nuclteotides
        """
        if not lst:
            sites = [sites]
        new_sites = []
        still_degen = False
        for site in sites:
            for degen, bases in self.iupac.iteritems():
                if degen in site:
                    still_degen = True
                    for base in bases:
                        # replace only first degen base we find; let recursion
                        # deal with those that follow
                        new_sites.append(site.replace(degen, base, 1))
            if not still_degen:
                new_sites.append(site)
        if still_degen:
            new_sites = self._convert_degenerate_primers(new_sites, lst=True)
        new_sites = list(set(new_sites))
        return new_sites


class Overruns(dict):
    def __init__(self, tags, sites):
        dict.__init__(self)
        self.trim = tags.trim_overruns
        for combo in tags.name_d.keys():
            r1name, r2name = combo.split(',')
            r1t = tags.seq_d[r1name]
            r2t = tags.seq_d[r2name]
            both = {
                'r1': [OverrunsMeta(r2t, i.string) for i in sites.r2],
                'r2': [OverrunsMeta(r1t, i.string) for i in sites.r1]
            }
            dict.__setitem__(self, combo, both)

    def __getitem__(self, key):
        val = dict.__getitem__(self, key)
        return val


class OverrunsMeta:
    def __init__(self, tag, site):
        self.string = reverse_complement("{0}{1}".format(tag, site)).upper()
        self.string_len = len(self.string)
        self.regex = re.compile("({0})".format(self.string), re.IGNORECASE)


class Storage:
    def __init__(self, tags):
        self.output = defaultdict(lambda: defaultdict(str))
        tld = 'splitaake-output'
        os.mkdir(tld)
        for indiv in tags.name_d.values():
            fpath = os.path.join(tld, indiv)
            os.mkdir(fpath)
            for read in ['R1', 'R2']:
                name = "{0}-{1}.fastq.gz".format(indiv, read)
                self.output[indiv][read] = self._open_zip_file(
                    os.path.join(fpath, name)
                )

    def _open_zip_file(self, pth):
        return gzip.open(pth, 'wb')

    def close(self):
        for name, files in self.output.iteritems():
            for file in files.values():
                file.close()


class Parameters:
    """All run parameters"""
    def __init__(self, conf):
        self.db = Db(conf)
        self.parallelism = Parallelism(conf)
        self.reads = Reads(conf)
        self.quality = Quality(conf)
        self.tags = Tags(conf)
        self.sites = Sites(conf)
        self.overruns = Overruns(self.tags, self.sites)
        self.storage = Storage(self.tags)
