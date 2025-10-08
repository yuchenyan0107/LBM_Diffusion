#! /usr/bin/env python

"""
Print the commands as executed into a file.

When compiling with a single thread (-j1), this yields
a simple sequence of commands to build the project.
"""

from waflib import Context

orig_exec = Context.Context.exec_command

def exec_command(self, cmd, **kw):
	kw['shell'] = isinstance(cmd, str)

	if isinstance(cmd, str):
		txt = cmd
	else:
		txt = ' '.join(repr(x) if ' ' in x else x for x in cmd)

	with open('command_sequence.txt', 'a') as seqfile:
		seqfile.write('{}\n'.format(txt))
	return orig_exec(self,cmd,**kw)

Context.Context.exec_command = exec_command
