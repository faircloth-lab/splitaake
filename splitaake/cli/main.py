#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 02 January 2014 09:13 PST (-0800)
"""

from __future__ import absolute_import
import sys
import argparse

from splitaake.cli import main_help
#from phyluce.cli import main_amplicon
#from phyluce.cli import main_rad
#from phyluce.cli import main_split
#from phyluce.cli import main_db


def main():
    # print same output as help if no arguments
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    # setup main program args
    p = argparse.ArgumentParser(
        description="splitaake is a software package for splitting or "
                    "demultiplexing massively parallel sequencing data "
                    "that have been tagged combinatorially."
    )
    p.add_argument(
        "-V", "--version",
        action="version",
        version="splitaake {}".format("1.0.0")
    )

    sub_parsers = p.add_subparsers(
        metavar="command",
        dest="cmd",
    )

    main_help.configure_parser(sub_parsers)

    try:
        import argcomplete
        argcomplete.autocomplete(p)
    except ImportError:
        pass
    except AttributeError:
        # On Python 3.3, argcomplete can be an empty namespace package when
        # argcomplete is not installed. Not sure why, but this fixes it.
        pass

    args = p.parse_args()

    try:
        args.func(args, p)
    except RuntimeError as e:
        sys.exit("Error: %s" % e)
    except Exception as e:
        raise


if __name__ == "__main__":
    main()
