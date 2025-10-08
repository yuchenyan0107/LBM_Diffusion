#!/usr/bin/env python
# encoding: utf-8
# Copyright (c) 2017-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import argparse
import glob
from module import Module

__author__ = "Peter Vitt"
__copyright__ = "Copyright 2017, University of Siegen"
__version__ = "0.0.1"

def _check_implicit(module, args):
    pass

def _check_private(module, args):
    pass

def _check_publics(module, args):
    publics = module.getPublicSymbols()
    prefix = args['prefix']
    violations = [p for p in publics if not p.startswith(prefix)]
    for v in violations:
        print('{0}: The public symbol {1} is not prefixed correctly'.format(
            module.Name,
            v))

def _check_uses(module, args):
    uses = module.getUses()
    lines = module.getCombinedLines()
    for u in uses:
        if uses[u] is None:
            continue
        for ustmt in uses[u]:
            nUses = len([True for l in lines if ustmt.lower() in l.lower()])
            if 'operator' in ustmt.lower():
                # We imported an operator here, thus we need to look for the
                # operator itself instead of the operator(.xxx.) string.
                op = '.{0}.'.format(ustmt.split('.')[1].lower())
                nUses = len([True for l in lines if op in l.lower()])
            if nUses == 1:
                print("{0}: Import of {1} from {2} is not used within the module".format(
                    module.Name,
                    ustmt,
                    u))
            elif nUses == 0:
                print("{0}: There appears to be no import for {1} from {2} anymore. Something is really fishiy here.".format(
                    module.Name,
                    ustmt,
                    u))

def _check_type_prefix(module, args):
    types = module.getTypes()
    prefix = args['prefix']
    if prefix:
        violations = [t for t in types if not t.startswith(prefix)]
        for v in violations:
            print('{0}: The type {1} is not prefixed correctly'.format(
                module.Name,
                v))


checks = {
    'implicit': _check_implicit,
    'private': _check_private,
    'publics': _check_publics,
    'types': _check_type_prefix,
    'uses': _check_uses
}

def _get_checks(args):
    """
    This routine creates a dictionary of checks that were requested by the user.

    When dish is called, the user can specify which checks should be performed.
    For each check selected, we add an entry to the dictionary that contains the
    name of the check as key and a function that accepts the module as argument
    to call the check.
    """
    if 'all' in args['checks']:
        selected_checks = checks
    else:
        selected_checks = {}
        for c in args['checks']:
            selected_checks[c] = checks[c]
    return selected_checks

def main():
    parser = argparse.ArgumentParser(
            description='Check fortran files for coding style compliance')
    parser.add_argument("files", nargs='+', help='files to parse')
    parser.add_argument("--prefix", 
            default="",
            help="The prefix used for types, modules and routines")

    available_checks = ['all']
    [available_checks.append(c) for c in list(checks.keys())]
    parser.add_argument("-c", "--checks",
            action='append',
            choices=available_checks,
            help="Checks to perform")
    args = vars(parser.parse_args())
    if args['checks'] is None:
        args['checks'] = ['all']
    selected_checks = _get_checks(args)

    for g in args['files']:
        for f in glob.glob(g):
            m = Module(f)
            for check, func in selected_checks.items():
                func(m, args)


if __name__ == "__main__":
    main()
