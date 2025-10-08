#! /usr/bin/env python
# encoding: utf-8

import os
import re
import sys
from waflib import Configure, Context, Logs, Utils
from waflib.Task import Task
from waflib.TaskGen import feature, extension, before_method

def options(opt):
    opt.add_option('--coco_reports', action='store_true', default=False,
                   help='Activate output of CoCo reports')
    opt.add_option('--coco_set', action='store', help='File with coco settings')

def find_in_parents(node, name):
    '''Find the given name in the given node directory or one of its parents'''
    np = node
    found = np.find_resource(name)
    while (not found) and np.parent:
        np = np.parent
        found = np.find_resource(name)
    return found

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
            conf.cmd_and_log(buildcoco)
            Logs.warn('Done with compilation of coco.')

        conf.env.COCO = [ mycoco ]

    conf.env.COCO[0] = os.path.abspath(conf.env.COCO[0])
    Logs.warn('Found CoCo: '+ ' '.join(conf.env.COCO))

    if conf.options.coco_set:
        # Using find_in_parents here to check for the definition in the top most
        # directory if there are nested wscripts.
        cocoset_src = find_in_parents(conf.path, conf.options.coco_set).path_from(conf.path)
        if cocoset_src:
            conf.env["COCOSET"] = cocoset_src
        else:
            Logs.warn(f'Coco settings file "{conf.options.coco_set}" not found!')

    if not conf.options.coco_reports:
        # Make coco silent, if not explicitly asked for reports:
        if conf.env.COCOFLAGS:
            conf.env.COCOFLAGS.insert(0, '-s')
            conf.env.COCOFLAGS.append('-ad')
        else:
            conf.env.COCOFLAGS = ['-s', '-ad']

def build(bld):
    '''Need to copy the coco.set file if requested, before any coco task starts.'''
    cocoset = None
    if bld.env.COCOSET:
        cocoset = bld.env.COCOSET
    if bld.options.coco_set:
        cocoset = find_in_parents(bld.path, bld.options.coco_set).path_from(bld.path)
    if cocoset:
        topnode = bld.root.find_node(bld.top_dir)
        bld(rule = "cp ${SRC} ${TGT}",
            name = 'bld_cocoset',
            source = bld.path.find_resource(cocoset),
            target = topnode.find_or_declare('coco.set'),
            before = 'coco')

@feature('coco')
@before_method('process_source', 'process_rule')
def add_cocoset_dep(self):
    '''Ensure that the task to copy the coco.set file is posted!

       This is needed for partial builds, where the coco.set target
       may not be included.
       It's still needed for the coco tasks, though and has to be
       posted if any coco task is to be built.
    '''
    if self.env.COCOSET:
        ccc_set = self.bld.get_tgen_by_name('bld_cocoset')
        ccc_set.post()

class coco(Task):
    run_str = '${COCO} ${COCOFLAGS} ${TGT} ${SRC}'
    color = 'BLUE'
    ext_in = ['.fpp']
    ext_out = ['.f90']

    def scan(self):
        if Logs.verbose:
            Logs.debug('deps: dependencies for %r' % self.inputs)

        INC_REGEX = r"""^\?\?\s*INCLUDE\s+(?:\w+_)?[<"'](.+?)(?=["'>])"""
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
        if self.env.COCOSET:
            deplist.append(bld_dir.find_node('coco.set'))
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
            found = find_in_parents(node.parent, x)
            names.append(found)
            iterativecheck(found,re_inc,incs,names,seen)
    return names
