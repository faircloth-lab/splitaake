#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 02 January 2014 09:29 PST (-0800)
"""

from __future__ import absolute_import
import sys

from splitaake.cli import sub_split_many


descr = "Split many files of gzip reads."


def configure_parser(sub_parsers):
    if len(sys.argv) == 2:
        sys.argv.append("-h")
    p = sub_parsers.add_parser(
        'split',
        description=descr,
        help=descr
    )

    sub_parsers = p.add_subparsers(
        metavar="command",
        dest="cmd",
    )

    sub_split_many.configure_parser(sub_parsers)
