#! /usr/bin/env python
# encoding: utf-8

import os
import re
import sys
from waflib import Configure, Context, Logs, Utils
from waflib.Task import Task
from waflib.TaskGen import feature, extension

def configure(conf):
    conf.add_os_flags('COCOFLAGS')
    conf.find_program('coco', var='COCO', mandatory=False)
    if not conf.env.COCO:
        mycoco = os.path.abspath(Context.run_dir)+'/coco'

        if not conf.find_program(mycoco, mandatory=False):
            Logs.warn('Need to compile coco itself...')
            if not conf.env.FC:
                conf.fatal('No Fortran Compiler defined yet!')
            buildcoco = ' '.join(conf.env.FC)+' -o '+mycoco+' '+Context.waf_dir+'/waflib/extras/coco.f90'
            conf.exec_command(buildcoco)
            Logs.warn('Done with compilation of coco.')

        conf.env.COCO = [ mycoco ]

    conf.env.COCO[0] = os.path.abspath(conf.env.COCO[0])
    Logs.warn('Found CoCo: '+ ' '.join(conf.env.COCO))

@feature('coco')
def dummy(self):
    pass

class coco(Task):
    run_str = '${COCO} ${COCOFLAGS} ${TGT} ${SRC}'
    color = 'BLUE'
    ext_in = ['.fpp']
    ext_out = ['.f90']

    def scan(self):
        if Logs.verbose:
            Logs.debug('deps: dependencies for %r' % self.inputs)
        # print("Looking for dependencies of:")
        # print(self.inputs[0])

        INC_REGEX = """^\?\?\s*INCLUDE\s+(?:\w+_)?[<"'](.+?)(?=["'>])"""
        re_inc = re.compile(INC_REGEX, re.I)
        node = self.inputs[0]

        dnodes = []
        incs = []
        names = []
        seen = []
        y = 0
        deplist = iterativecheck(node,re_inc,incs,names,seen)
        if Logs.verbose:
            Logs.debug('deps: %r' % deplist)

        bld_dir = self.generator.bld.bldnode
        for x in deplist:
            inc_node = bld_dir.find_node(incs[y])
            if not inc_node:
                inc_node = bld_dir.make_node(incs[y])
                srctxt=x.read()
                inc_node.parent.mkdir()
                if Logs.verbose:
                    Logs.debug('build COCO INCLUDED file '
                               + inc_node.abspath()
                               + ' from ' + x.abspath())
                inc_node.write(srctxt)
            inc_node.sig = inc_node.cache_sig = Utils.h_file(inc_node.abspath())
            dnodes.append(inc_node)
            y = y+1
        return(deplist, names)

@extension('.fpp')
def process_coco(self, node):
    "Preprocess the .fpp files with coco."
    f90_node = node.change_ext('.f90')
    self.create_task('coco', node, [f90_node])
    # If this node is to be compiled afterwards
    # (TaskGenerator has the fc feature),
    # add it to the list of sources to process.
    if "fc" in self.features:
        self.source.append(f90_node)

# Look for dependency files recursively
def iterativecheck(node,re_inc,incs,names,seen):
    if node:
        txt = node.read()
        for line in txt.splitlines():
            m = re_inc.search(line)
            if m:
                incs.append(m.group(1))
        for x in incs:
            if x in seen:
                continue
            seen.append(x)
            np = node.parent
            found = np.find_resource(x)
            # search for depend coco file recursively
            while not found:
                found=np.parent.find_resource(x)
                np=np.parent
            # print(found)
            names.append(found)
            iterativecheck(found,re_inc,incs,names,seen)
    return names
