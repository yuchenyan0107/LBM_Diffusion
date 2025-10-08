def options(opt):
    opt.add_option('--fordonline', action='store_true', default=False,
                   help='Use links to online documentation to refer to external projects in FORD')

def configure(conf):
    conf.env.fordonline = conf.options.fordonline

def build(bld):
    bld.env.fordonline = bld.env.fordonline or bld.options.fordonline

def gendoc(task):
    """ A Task rule to generate Documentation with FORD

        The rule expects the location of the mainpage with the
        FORD project description as task generator argument 'mainpage'.
        The location of the index.html that is to be generated, should
        be the first target argument.
        All unique parent nodes of the files given as sources will be
        used as src_dir for the FORD command.
        Additional paths can be provided via the `src_paths` option,
        which has to be a list of paths to look for sources in.
        Pathes to external projects in the local file system can be
        provided to the generator via the 'extern' argument.
        URLs of external projects that are available online can be
        provided via the the 'extern_urls' arguments.
    """
    import os
    tgt_path = os.path.dirname(task.outputs[0].abspath())
    src_paths = set()
    use_urls = task.env.fordonline
    if not hasattr(task.generator, 'extern'):
        task.generator.extern = []
    if not hasattr(task.generator, 'extern_urls'):
        task.generator.extern_urls = []
    if hasattr(task.generator, 'src_paths'):
        for sp in task.generator.src_paths:
            src_paths.add(sp)
    for srcfile in task.inputs:
        if srcfile not in task.generator.extern:
            src_paths.add(srcfile.parent.abspath())
    cmd = ['ford', '--externalize', '-o', tgt_path]
    if hasattr(task.env, 'revision_string'):
        cmd.append('-r')
        cmd.append(task.env.revision_string)
    for spath in src_paths:
        cmd.append('-d')
        cmd.append(spath)
    if task.env.fordonline:
        for elink in task.generator.extern_urls:
            cmd.append('-L')
            cmd.append(elink)
    else:
        for elink in task.generator.extern:
            cmd.append('-L')
            cmd.append(elink.parent.abspath())
    cmd.append(task.generator.mainpage)
    return task.exec_command( cmd, shell=False, cwd=task.generator.bld.top_dir )
