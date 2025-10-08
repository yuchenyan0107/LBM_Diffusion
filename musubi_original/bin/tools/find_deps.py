#! /usr/bin/env python
import sys
import glob
import fileinput
import re
import pygraphviz as pgv
from optparse import OptionParser

USE_REGEX = """(?:^|;)\s*USE(?:\s+|(?:(?:\s*,\s*(?:NON_)?INTRINSIC)?\s*::))\s*(\w+)"""
MOD_REGEX = """(?:^|;)\s*MODULE(?!\s*PROCEDURE)(?:\s+|(?:(?:\s*,\s*(?:NON_)?INTRINSIC)?\s*::))\s*(\w+)"""
PRG_REGEX = """(?:^|;)\s*PROGRAM\s*(\w+)"""

re_use = re.compile(USE_REGEX, re.I)
re_mod = re.compile(MOD_REGEX, re.I)
re_prg = re.compile(PRG_REGEX, re.I)

cm = 'NONE'
GV = pgv.AGraph(directed=True)

nodes = dict()

usage = "usage: %prog [options] files"
parser = OptionParser(usage=usage)
parser.add_option("-e", 
                  "--exclude", 
                  metavar="EXCLUDE",
                  help="Excludes modules starting with EXCLUDE")
parser.add_option("-i", 
                  "--include", 
                  metavar="INCLUDE",
                  help="Include only modules starting with INCLUDE")

(options, args) = parser.parse_args()

for line in fileinput.input(args):
	m = re_prg.search(line)
	if m:
		cm = m.group(1)
		print 'PROGRAM ' + cm
		GV.add_node(cm)

	m = re_mod.search(line)
	if m:
		cm = m.group(1)
		print 'MODULE ' + cm
		GV.add_node(cm)

	m = re_use.search(line)
	if m:
                if options.include and not m.group(1).startswith(options.include):
                        continue
                if options.exclude and m.group(1).startswith(options.exclude):
                        continue
		GV.add_edge(cm, m.group(1))
                if cm in nodes:
                    nodes[cm].append(m.group(1))
                else:
                    nodes[cm] = [m.group(1)]
		print '  uses ' + m.group(1)

GV.layout('dot')
GV.draw('moddeps.png')

def find_circle(start, nodename, nodelist):
    # First we have to check whether the module is one we found during the file
    # parse. If not, then we cannot detect any circle. If we know the module,
    # proceed with the search.
    if nodename in nodelist:
        # We know the node, so we can check all the used modules
        for child in nodelist[nodename]:
            # Is the node already in the start list?
            if child in start:
                # If so, we have a circular dependency
                print 'Circular reference found: [%s]' % ', '.join(map(str,start))
                return True

            # We haven't passed the current node yet, so we have to check it's 
            # decendants
            new_start = list(start)
            new_start.append(child)
            if find_circle(new_start, child, nodelist):
                return True

    return False

for node in nodes.keys():
    find_circle([], node, nodes)

