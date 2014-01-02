#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 02 January 2014 09:15 PST (-0800)
"""

from __future__ import absolute_import


descr = "Get help info on a splitaake command."


def configure_parser(sub_parsers):
    p = sub_parsers.add_parser(
        'help',
        description=descr,
        help=descr
    )
    p.add_argument(
        "command",
        metavar="COMMAND",
        action="store",
        nargs='?',
        help="print help information for COMMAND "
             "(same as: phyluce COMMAND -h)",
    )
    p.set_defaults(func=execute)


def execute(args, parser):
    if not args.command:
        parser.print_help()
        return

    import sys
    import subprocess

    subprocess.call([sys.executable, sys.argv[0], args.command, '-h'])
