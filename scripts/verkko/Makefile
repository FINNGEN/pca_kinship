# This is a Makefile for doing various tests related to verkko, like
# building documentation.  It is in a Makefile for concise scripting,
# but could be moved else where later.  Maybe it would better fit in a
# setup.py file.


.PHONY: default test docs coverage cron

default:
	@echo "Usage:"
	@echo "  make clean - remove all temporary files created by these commands"
	@echo "  make test - run unit tests (making docs/test-results.html if the)"
	@echo "              html-output module is installed)"
	@echo "  make docs - make HTML documentation in build/html/"
	@echo "  make coverage - makes test coverage results in ./docs/coverage/"

clean:
	rm -rf docs/build/
	rm -rf docs/coverage/
	rm -f docs/test-results.html

# Python wrapper which disables X in pylab/pyplot
PYTHON = python ./misc/python_no_pyplot.py

# For HTML output this must be installed:
#  https://github.com/cboylan/nose-html-output
NOSE_HTML = $(if $(shell nosetests -p | grep html-output), \
                   --with-html-output --html-out-file=docs/test-results.html)

# These files should NOT be imported under any circumstances, due to
# side-effects.  Using this variable should be a last-resort for
# modules which intrinsically function via side-effects.
EXCLUDE_FILES = misc/interactnow.py

# Run tests
test:
	$(PYTHON) `which nosetests` . $(NOSE_HTML)

# Generate all docs
docs:
	mkdir -p docs/
	$(PYTHON) docs/apidoc.py ../verkko/ -o ./docs/api/ -f -d 0 --separate \
                 $(EXCLUDE_FILES)     # exclude paths
#	PYTHONPATH=$PYTHONPATH:. sphinx-build -b html ./docs/ ./docs/build/html/
#	This is needed in order to handle 'import pylab' in scripts.
#	PYTHONPATH=..:$$PYTHONPATH python -c 'import matplotlib ; matplotlib.use("Agg"); import sphinx ; sphinx.main(argv="sphinx-build -E -a -b html ./docs/ ./docs/build/html/".split())'
	PYTHONPATH=..:$$PYTHONPATH $(PYTHON) `which sphinx-build` -E -a -b html ./docs/ ./docs/build/html/


# Make the coverage tests in ./docs/coverage/
coverage:
	$(PYTHON) `which nosetests` \
	--with-coverage ../verkko/ \
	--cover-erase --cover-package=verkko \
	--cover-html --cover-html-dir=docs/coverage/
#	--cover-inclusive

# Automatic script to fetch, update.  The weird usage of '-f
# $(MAKEFILE)' here is so that I can make a copy of the Makefile to
# use with cron.  Git updates won't automatically change the cron
# script.  This doesn't provide any *real* security, but until there
# is more isolation I want to go with it.  You should "cp Makefile
# Makefile.local" manually on each significant update.  The cron
# command should be:
#  ... cd /path/to/this/ && make -f Makefile.local cron > cron.output 2>&1
MAKEFILE = $(MAKEFILE_LIST)
cron:
	git fetch
	git checkout origin/master
	make -f $(MAKEFILE) clean || true
	make -f $(MAKEFILE) test || true
	make -f $(MAKEFILE) coverage || true
	make -f $(MAKEFILE) docs || true
	git checkout master

