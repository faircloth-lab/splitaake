#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 02 January 2014 09:31 PST (-0800)
"""

from __future__ import absolute_import


from splitaake.split import many


descr = "Split a directory containing many files of gzip reads."


def configure_parser(sub_parsers):
    sp = sub_parsers.add_parser(
        "many",
        description=descr,
        help=descr
    )

    sp.add_argument(
        "--reads",
        required=True,
        help="The input fastq files"
    )
    sp.add_argument(
        "--tagmap",
        required=True,
        help="A conf file mapping name:tag sequence"
    )
    sp.add_argument(
        "--output",
        required=True,
        help="The path to the output directory"
    )
    sp.add_argument(
        "--section",
        default="tags",
        help="If multisection conf, use this section"
    )
    sp.add_argument(
        "--no-correct",
        action="store_true",
        default=False,
        help="Do not correct certain tags " +
        "in the [no correct] section of the tagmap config file",
        )
    sp.add_argument(
        "--tag-length",
        type=int,
        default=None,
        help="Length of tags if tags are < read length"
    )
    sp.add_argument(
        "--min-qual",
        type=int,
        default=10,
        help="The minimum single-base quality (Q) to accept"
    )
    sp.add_argument(
        "--min-mean-qual",
        type=int,
        default=20,
        help="The minimum average quality (Q) to accept"
    )
    sp.add_argument(
        "--hamming",
        action="store_true",
        default=False,
        help="Use the Hamming (substitution-only) distance",
    )
    sp.add_argument(
        "--quiet",
        action="store_true",
        default=False,
        help="Run in quiet mode",
    )

    sp.set_defaults(func=split_many_reads)


def split_many_reads(args, parser):
    many.main(args, parser)
