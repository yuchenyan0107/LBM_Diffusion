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
import os
import sequencer

__author__ = "Peter Vitt"
__copyright__ = "Copyright 2017, University of Siegen"
__version__ = "0.0.1"

class Module(object):

    def __init__(self, filepath):
        self._file = filepath
        self.Name = os.path.splitext(os.path.basename(filepath))[0]

    def _get_lines(self):
        """Read all lines from _file"""
        # read all lines from file
        with open(self._file, 'r') as f:
            lines = f.readlines()
        # strip whitespace from line endings
        return [x.strip() for x in lines]

    def getPublicSymbols(self):
        lines = sequencer.make_long_lines(self._get_lines())

        all_private = any([True for l in lines if l == 'private'])
        if all_private:
            return sequencer.get_publics(lines)
        else:
            # Module is not private by default, thus we have to append all
            # module's top level items
            print("Not yet implemented: Module "
                    + self.Name
                    + " is not private by default, thus we"
                    " have to append all module's top level items!")
            return []


    def getUses(self):
        lines = self._get_lines()
        return sequencer.get_uses(lines)

    def getCombinedLines(self):
        return sequencer.make_long_lines(self._get_lines())

    def getTypes(self):
        lines = self._get_lines()
        return sequencer.get_types(lines)
