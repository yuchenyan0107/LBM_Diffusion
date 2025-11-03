def options(opt):
	"""
	Provide the ``--showtests`` command-line option. To allow the display
        of testcases which were successful.
	"""
	opt.add_option('--showtests', action='store_true', default=False,
                       help='Show output of successful unit tests.',
                       dest='show_tests')
	opt.add_option('--mpicmd', action='store', default='mpiexec',
                       help='MPI command to execute a parallel application.',
                       dest='mpicmd')

def summary(bld):
	"""
	Get the test results from last line of output::

		Fortran applications can not return arbitrary return codes in
		a standarized way, instead we use the last line of output to
		decide the outcome of a test: It has to state "PASSED" to count
		as a successful test.

		Otherwise it is considered as a failed test. Non-Zero return codes
		that might still happen are also considered as failures.

		Display an execution summary:

		def build(bld):
			bld(features='cxx cxxprogram test', source='main.c', target='app')
			from waflib.extras import utest_results
			bld.add_post_fun(utest_results.summary)
	"""
	from waflib import Logs
	import sys

	lst = getattr(bld, 'utest_results', [])

	# Check for the PASSED keyword in the last line of stdout, to
	# decide on the actual success/failure of the test.
        nlst = []
	for (f, code, out, err) in lst:
		ncode = code
		if not code:
			if sys.version_info[0] > 2:
				lines = out.decode('ascii').splitlines()
			else:
				lines = out.splitlines()
			if lines:
				ncode = lines[-1].strip() != 'PASSED'
			else:
				ncode = True
		nlst.append([f, ncode, out, err])
	lst = nlst

	if lst:
		Logs.pprint('CYAN', 'Utests execution summary:')

		total = len(lst)
		tfail = len([x for x in lst if x[1]])

		if tfail > 0:
			failcolor = 'RED'
		else:
			failcolor = 'CYAN'

		if tfail == total:
			passcolor = 'RED'
		else:
			passcolor = 'CYAN'

		Logs.pprint(passcolor, '  tests that pass %d/%d' % (total-tfail, total))
		for (f, code, out, err) in lst:
			if not code:
				Logs.pprint('GREEN', '    %s' % f)
				if bld.options.show_tests:
					if out:
						if sys.version_info[0] > 2:
							Logs.pprint('NORMAL', out.decode('ascii'))
						else:
							Logs.pprint('NORMAL', out)
					if err:
						if sys.version_info[0] > 2:
							Logs.pprint('YELLOW', err.decode('ascii'))
						else:
							Logs.pprint('YELLOW', err)
		Logs.pprint('NORMAL', '')

		Logs.pprint(failcolor, '  tests that fail %d/%d' % (tfail, total))
		for (f, code, out, err) in lst:
			if code:
				Logs.pprint('RED', '    %s' % f)
				if out:
					if sys.version_info[0] > 2:
						Logs.pprint('NORMAL', out.decode('ascii'))
					else:
						Logs.pprint('NORMAL', out)
				if err:
					if sys.version_info[0] > 2:
						Logs.pprint('YELLOW', err.decode('ascii'))
					else:
						Logs.pprint('YELLOW', err)
		if tfail == 0:
			Logs.pprint('GREEN', 'Utests successful!')
		else:
			bld.fatal('Utests FAILED!')

def utests(bld, use, path='utests', preprocessor=None):
	"""
	Define the unit tests from the programs found in the utests directory.
	"""
	from waflib import Options
	import os
	preprocessor_string = ''
	if preprocessor is not None:
		preprocessor_string = '{0} '.format(preprocessor)
	ppsources = bld.path.ant_glob(path + '/*_test.fpp')
	for utest in ppsources + bld.path.ant_glob(path + '/*_test.f90'):
		nprocs = search_procs_in_file(utest.abspath())
		if int(nprocs) > 0:
			bld(
			    features = preprocessor_string+'fc fcprogram test',
			    source = utest,
			    use = use,
			    ut_cmd='{0} -n {1} %s'.format(Options.options.mpicmd,nprocs),
			    install_path = None,
			    target = os.path.basename(utest.change_ext('').abspath()))
		else:
			bld(
			    features = preprocessor_string+'fc fcprogram test',
			    source = utest,
			    use = use,
			    install_path = None,
			    target = os.path.basename(utest.change_ext('').abspath()))

def search_procs_in_file(fname):
	import re
	PROC_REGEX = r'\s*!MPI!\s*nprocs\s*=\s*(\w+)'
        re_proc = re.compile(PROC_REGEX, re.I)
	with open(fname) as f:
		for line in f:
			m = re_proc.match(line)
			if m:
				return m.group(1)
		return '0'
