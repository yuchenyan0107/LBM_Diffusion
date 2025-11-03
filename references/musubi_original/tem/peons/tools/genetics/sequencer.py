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
__author__ = "Peter Vitt"
__copyright__ = "Copyright 2017, University of Siegen"
__version__ = "0.0.1"

def _get_public(line):
    """
    Determines all public symbols a line contains.
    """
    if "public" in line:
        return line.replace("public","").replace("::","").strip().split(",")

def _get_use(line):
    """
    Determines whether this line is a use statement and returns the module name
    as well as all explicitly stated imports, if any.
    """
    if line.startswith('use '):
        iuse = line.index('use')
        try:
            icolon = line.index(',')
        except ValueError:
            icolon = len(line)
            pass
        try:
            ionly = line.index('only:')
        except ValueError:
            ionly = None
            pass
        module = line[iuse+3:icolon].strip()
        if ionly is not None:
            uses = line[ionly+5:].strip().split(',')
            # filter operatores (they contain ()) and remove whitespace
            uses = [u.strip() for u in uses if not ("(" in u or ")" in u)]
        else:
            uses = None
#        for i in range(len(uses)):
#            uses[i] = uses[i].strip()
        return True, module, uses
    return False, None, None


def _strip_comments(lines):
    """
    Removes all lines that contain only comments.
    """

    code_lines = []
    for line in lines:
        # Check whether the whole line is a comment
        if not line.strip().startswith('!'):
            # remove comments that do not fill the whole line
            p = line.find('!')
            if p > 0:
                line = line[:p]
            code_lines.append(line)
    return code_lines


def make_long_lines(lines):
    """
    Takes a bunch of lines and collapses all continued lines into a single one.
    """

    long_lines = []
    for line in _strip_comments(lines):
        if len(line) > 0:
            if line[0] == '&':
                lli = len(long_lines)-1
                # continuation line, add to last line
                long_lines[lli] += line.replace('&', '').strip()
                long_lines[lli] = long_lines[lli].replace('&','')
            else:
                long_lines.append(line)
    return long_lines


def get_publics(lines):
    """
    Extract public symbols
    """

    long_lines = make_long_lines(_strip_comments(lines))

    publics = []
    for l in lines:
        if "public" in l:
            [publics.append(p.strip()) for p in _get_public(l)]

    return publics

def get_uses(lines):
    """
    Extract use statements.

    This routine extracts all use statements and stores them in a dictionary.
    The module name is the key, and all explicitly stated imports are stored in
    a list as the key's value.

    """

    long_lines = make_long_lines(_strip_comments(lines))
    uses = {}

    for l in long_lines:
        found, module, ustmts = _get_use(l)
        if found:
            uses[module] = 0
            uses[module] = ustmts
    return uses


def get_types(lines):
    """
    Extract all type definitions.

    This routine extracts all type definitions and returns their names in a
    list.

    """
    import re

    long_lines = make_long_lines(_strip_comments(lines))

    types = []
    for l in long_lines:
        if l.startswith("type"):
            p = l.find("type")
            t = l[p:].strip()
            # check for single word
            if re.match(r'\A[\w-]+\Z', t):
                types.append(t)
    return types



