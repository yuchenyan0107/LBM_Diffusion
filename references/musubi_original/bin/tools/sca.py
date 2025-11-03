#!/usr/bin/env python3

import os
import argparse
import logging

from fparser.common.readfortran import *
from fparser.two.parser import ParserFactory
from fparser.two.utils import walk
from fparser.two.Fortran2003 import *

parser = argparse.ArgumentParser(description='Perform a static code analysis on fortran files.')
parser.add_argument('files', help='The absolute or relative path of the files to parse', nargs='+')
args = parser.parse_args()

f2008_parser = ParserFactory().create(std="f2008")

def iequal(a, b):
    try:
        return a.upper() == b.upper()
    except AttributeError:
        return a == b

def PrivateModules_Rule(tree):
    # check whether this is a module
    res = walk(tree, Module)
    if (len(res)>0):
        # Check for exactly one private statement
        res = walk(tree, Access_Stmt)
        private = [n for n in res if iequal(n.children[0], "PRIVATE")]
        if len(private) < 1:
            raise AssertionError("No private access statement found")
        elif len(private) > 1:
            raise AssertionError("Too many private access statements found")

def Prefix_Rule(tree):
    # check whether this is a module
    res = walk(tree, Module)
    if (len(res)>0):
        res = walk(tree, Access_Stmt)
        publics = [n for n in res if iequal(n.children[0], "PUBLIC")]
        for p in publics:
            access_id_list = p.children[1] # the part after the ::
            for aid in access_id_list.children:
                if not str(aid)[:4] == "atl_":
                    raise AssertionError("Public entity {} not prefixed with atl_".format(aid))

rules = (PrivateModules_Rule, Prefix_Rule)

for f in args.files:
    reader = FortranFileReader(os.path.abspath(f))
    tree = f2008_parser(reader)

    for rule in rules:
        try:
            rule(tree)
        except AssertionError as ae:
            print(ae)
