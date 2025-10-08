import subprocess
import argparse
import glob
import csv

import os
import sys
# Make the working directory available to import modules from.
sys.path.append(os.getcwd())

def authors_of_file(f):
    ''' Parse the history of the given file f in the repository
        and return the contributing authors as a dictionary with lists
        comprised of the last log line they appeared in (first commit) and a
        set of contributing years.

        Users might employ an authormap.py file in the working directory
        to indicate aliases for given author strings.
        The authormap file needs to provide a dict named alias of the following
        form:
        alias = {
            "Name found in log": "Identification of author",
            "Other name in log": "Identification of author"
        }
    '''
    try:
        import authormap
        author_alias = authormap.alias
    except:
        author_alias = {}

    hg_out = subprocess.check_output(
                 ['hg', 'log', '-f',
                  '--template', r'''{date(date, '%Y')},{author}\n''', f])
    logreader = csv.reader(hg_out.decode('utf-8').split('\n'), delimiter=',')
    authors = {}
    rownum = 0
    for row in logreader:
        # Keep track of the row number to sort resulting authors list according
        # to log.
        rownum += 1
        if len(row) == 2:
            # Get the author name from the alias dictionary if it is present,
            # otherwise use the found string itself as author.
            authorname = author_alias.get(row[1], row[1])
            year = int(row[0])
            if authorname in authors:
                authors[authorname][0] = rownum
            else:
                authors[authorname] = [rownum, set()]
            authors[authorname][1].add(year)
    return authors


parser = argparse.ArgumentParser(description=
'''Finding author contributions to files.

The resulting list of authors will be sorted by first occurence in the
Mercurial history and prefixed by the string provided via the --prefix
option, which defaults to '! '.
You can provide an authormap.py file in the working directory to define
aliases for given authorstrings.
This module needs to define a dictionary name alias of the following form:
    alias = {
        "Name found in log": "Identification of author",
        "Other name in log": "Identification of author"
    }
''')

parser.add_argument('-l', '--list', action='store_true',
                    help=
'''Only list all found authors uniquely from all files without years and put
them into a template for an alias dictionary.
This is helpful to create an authormap.py file''')

parser.add_argument('-p', '--prefix', default='! ',
                    help=
''' Prefix to use for the copyright notices and license.
    Defaults to the '! ' Fortran comment.
''')

parser.add_argument('-L', '--license',
                    help=
''' File with the license to use when amending the given files.
    This file will be read, each line gets the prefix prepended and the message
    is written after the copyright lines.
''')

parser.add_argument('-M', '--modify', action='store_true',
                    help=
''' Modify the given files and put the copyright notices along with
    the license at the beginning of the file.
''')

parser.add_argument('files', metavar='FILE', nargs='+', help='files to process')
args = parser.parse_args()


all_authors = []
for fname in args.files:
    sources = glob.glob(fname)
    for f in sources:
        authors = authors_of_file(f)
        print(f)

        contributors = []

        for author,commitinfo in authors.items():
            years = commitinfo[1]
            if args.list:
                if author not in all_authors:
                     all_authors.append(author)
            else:
                ylist = list(years)
                year_order = sorted(ylist)
                year_str = "{0}Copyright (c) {1}".format(args.prefix,
                                                         year_order[0])
                prev_year = year_order[0]
                lastseq = 0
                for year in year_order[1:]:
                    if year == prev_year+1:
                        lastseq = year
                    elif lastseq > 0:
                        year_str = "{0}-{1}, {2}".format(year_str, lastseq,
                                                         year)
                        lastseq = 0
                    else:
                        year_str = "{0}, {1}".format(year_str, year)
                    prev_year = year
                if lastseq > 0:
                    year_str = "{0}-{1}".format(year_str, lastseq)
                year_str = "{0} {1}".format(year_str, author)
                contributors.append((commitinfo[0], year_str))

        if args.modify:
            with open(f, 'r') as sourcefile:
                origin = sourcefile.read().splitlines()
            import shutil
            shutil.copy(f, f+'.orig')
            target = open(f, 'w')

        # Print contributor lines sorted by their first appearance in the
        # history.
        for copyright in sorted(contributors, reverse=True):
            if args.modify:
                target.write(copyright[1].strip()+'\n')
            else:
                print(copyright[1].strip())

        if args.license:
            if args.modify:
                target.write(args.prefix.strip()+'\n')
                with open(args.license, 'r') as licfile:
                    license_text = licfile.read().splitlines()
                for licline in license_text:
                   target.write('{0}{1}'.format(args.prefix, licline).strip() + '\n')
            else:
                print(args.prefix.strip())
                with open(args.license, 'r') as licfile:
                    license_text = licfile.read().splitlines()
                for licline in license_text:
                   print('{0}{1}'.format(args.prefix, licline).strip())

        if args.modify:
            for line in origin:
                if line != '! See copyright notice in the COPYRIGHT file.':
                    target.write(line+'\n')


if args.list:
    print('alias = {')
    for author in sorted(all_authors):
        print('"{0}": "",'.format(author))
    print('}')
