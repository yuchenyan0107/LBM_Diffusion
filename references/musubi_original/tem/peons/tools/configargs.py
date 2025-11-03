# Copyright (c) 2016 Harald Klimach <harald.klimach@uni-siegen.de>
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
'''
This module provides a way to use argument parsing along with a configuration
parser.
It allows you to use a basic configuration file to define settings, instead of
pure definition on the command line.
This allows the user to provide typical settings in a file, but override them
at runtime by the command line.
Such a possibility is especially useful when dealing with many configuration
settings, where only a small subset might be of interest to vary.

The argument parsing is implemented by the Configurator class, which is derived
from the argparse.ArgumentParser.
Most aspects of the ArgumentParser remain, but the args have to be provided in
the initialization, and the parse_args only takes the namespace argument now.
In the initialization it is now possible to provide a default file name for the
configuration file to use.
The user may override the location of the actual configuration file by the
--config or -c option.
Note, that the reading of the config file will not fail on the definition of
unknown arguments. Those will just be ignored, and the application would need
to check for them afterwards if desired.

See the configparse documentation for some details on the file format.
In sample.conf is a small example how the config file may look like.
'''

__version__ = "1.0.0"

import argparse
import itertools
try:
  import configparser
except ImportError:
  import ConfigParser as configparser

class Configurator(argparse.ArgumentParser):
  ''' A class to represent a configuration that can be provided in a file or
      on the command line.
      It utilizes the argparse module to obtain configuration settings but also
      configparser to prescribe settings in a file.
      Configurator is derived from argparse.ArgumentParser, and most routines
      behave are the same, but a major difference is that the arguments need to
      be passed to the initialization already.
      The config file has to have a [config] section, only this section will be
      considered, other sections are ignored.
      If there is no section header at all in the file, the complete file will
      be considered as being the config section.

      Note that the '-c'/'--config' argument is preoccupied and can't be used
      as an option anymore. This argument describes the location of the config
      file.

      Example:

      Setup the Configurator with a list of arguments, if you do not provide
      a list of arguments, sys.argv will be used. Here we use the argument list
      to define the location of the configuration file, and to set the two
      options '-foo' and '--set'.
      Definitions in the command line arguments take precedence over both, the
      given defaults in add_argument and the settings in the configuration file.
      >>> conf = Configurator(args=['-c', 'sample.conf', '--foo', 'boo', '--set', 'go'])

      Parameter set by argument list.
      >>> conf.add_argument('--foo', default='bar')
      _StoreAction(option_strings=['--foo'], dest='foo', nargs=None, const=None, default='bar', type=None, choices=None, help=None, metavar=None)

      Not provided parameter, here the given default will be used.
      >>> conf.add_argument('--param', default='nix')
      _StoreAction(option_strings=['--param'], dest='param', nargs=None, const=None, default='nix', type=None, choices=None, help=None, metavar=None)

      Parameter that is set in the config file (sample.conf).
      >>> conf.add_argument('--tree', default='oak')
      _StoreAction(option_strings=['--tree'], dest='tree', nargs=None, const=None, default='oak', type=None, choices=None, help=None, metavar=None)

      Parameter that is set in the config file (sample.conf), but overriden by
      the argument list.
      >>> conf.add_argument('--set', default='empty')
      _StoreAction(option_strings=['--set'], dest='set', nargs=None, const=None, default='empty', type=None, choices=None, help=None, metavar=None)

      Actually parse the arguments and obtain their values.
      >>> conf.parse_args()
      Namespace(config='script.conf', foo='boo', param='nix', set='go', tree='elm')
  '''

  def __init__(self, default_file='script.conf', args=None, **kwargs):
    ''' (String, List, options of argparse.ArgumentParser) -> Configurator
        Initialize the configurator object with the default_file as name to
        choose as default for the configuration file. If it is not provided,
        'script.conf' will be used as default.
        args are the arguments that are to be parsed, if not provided, this
        will default to sys.args.
        All remaining arguments by keywords are passed on to argparse.ArgParser.
    '''

    # The following approach is taken from
    # http://stackoverflow.com/a/5826167/577108
    # Parse any config specification
    # We make this parser with add_help=False so that
    # it doesn't parse -h and print help.
    conf_parser = argparse.ArgumentParser(
        description=__doc__, # printed with -h/--help
        # Don't mess with format of description
        formatter_class=argparse.RawDescriptionHelpFormatter,
        # Turn off help, so we print all options in response to -h
        add_help=False
        )
    conf_parser.add_argument("-c", "--config",
                        help='''
Configuration file to define any of the available arguments.
Defaults to {0}.
                        '''.format(default_file),
                        metavar="CONFIGFILE", default=default_file)

    # Get the (possible) config file
    foundargs, self.remaining_argv = conf_parser.parse_known_args(args)

    self.config_file = foundargs.config

    # Flag to indicate whether a config file was loaded.
    self.loaded_config_file = False

    # Add the conf_parser to the kwargs to inherit help from the config_file
    # option.
    if 'parents' in kwargs:
      kwargs['parents'] += conf_parser
    else:
      kwargs['parents'] = [conf_parser]

    super(Configurator, self).__init__(**kwargs)


  def parse_args(self, namespace=None):
    ''' Convert argument strings to objects and assign them as attributes of
        the namespace.
        This operates on the remaining arguments after processing the config
        file argument already in the initializiation.
        The defaults are updated with the values from the config file, before
        passing the arguments on to parse_args of the argparse.ArgumentParser.
        Thus, values on the command line override values from the configuration
        file which in turn override the defaults given in the add_argument
        calls.
    '''

    # Update defaults according to the settings in the configuration file.
    config = configparser.SafeConfigParser()
    try:
      # Try to ordinarily read the configuration from the provided config file.
      read_configs = config.read([self.config_file])
    except configparser.MissingSectionHeaderError:
      # If the section header is missing from the file, insert it as first line
      try:
        # With Python 3 we can use the read_file attribute of the ConfigParser
        config.read_file(itertools.chain(['[config]'], open(self.config_file, 'r')))
        read_configs = self.config_file
      except AttributeError:
        # For Python 2 we have to use readfb, and modify the file handle
        # see: http://stackoverflow.com/a/2819788/577108
        config.readfp(_FakeSecHead(open(self.config_file, 'r')))
        read_configs = self.config_file

    if read_configs:
      self.loaded_config_file = True
      file_defaults = dict(config.items("config"))
      super(Configurator, self).set_defaults(**file_defaults)

    return super(Configurator, self).parse_args(self.remaining_argv, namespace)


class _FakeSecHead(object): 
  """ A class to provide the config file with an prepended section header, if 
      none is provided in the file itself. 
      Provided for Python 2 at: http://stackoverflow.com/a/2819788/577108 
  """ 
  def __init__(self, fp): 
    self.fp = fp 
    self.sechead = '[config]\n' 
 
  def readline(self): 
    if self.sechead: 
      try: 
          return self.sechead 
      finally: 
          self.sechead = None 
    else: 
      return self.fp.readline() 


if __name__ == '__main__':
  import doctest
  doctest.testmod()
