# Python style and documentation guide for the group code-library

The reason for aiming to proper style and documentation is the following:

> Code is read much more often than it is written.
> Thus readability counts.

We try to follow the PEP
[style](http://legacy.python.org/dev/peps/pep-0008/) and
[docstring](http://legacy.python.org/dev/peps/pep-0257/) guidelines,
of which the following are modified excerpts with minor modifications.

Try to follow these guides, but don't let them come too much on your way.
Sharing the code is more important!

**Note that** for many texteditors there exists tools for automatically checking/correcting
code style, which can make your life a lot easier. See for example [autopep8](https://github.com/hhatto/autopep8).


### Use 4 spaces for intendation!
(Python 3 does not even allow mixed use of tabs and spaces)

### Limit lines to 79 characters
* Possible to have multiple source code windows open at the same time
* Easier comparison with code difference tools

### Spacing
Use spacing to make code easy to read and beautiful.  Don't cram
everything into one line:

    a=1+2+(5*6*7)                # bad
    a = 1 + 2 + (5*6*7)          # better - spaces are guide for eyes
    if a == b  or  c in (5,6):   # spacing helps parse
        ....

Parallelism is important, too, it makes it easier to read.  For
example:

    first  = 1   # extra space so that it lines up
    second = 2

    if a == 'const_1':     b    += 1
    if a == 'const_other': cval += 1

### Imports
Put imports on separate lines:

#### Yes

```python
import sys
import os
```

#### No

```python
import sys, os
```

Wildcard imports are **NOT** preferred:

* However, from `pylab import *` might not be always that dangerous.


### Comments
Please comment your code reasonably and smartly.  You **don't** need to make
something for every single line:

    for comm in part1:        # for each community
        size += len(comm)     # calculate total size

If you use good names and code structure, code is almost
self-documenting at the line level (when it's not, *then* have line
comments).  Instead, consider comments as your abstract or
introduction, explaining higher level logic:

    # Find sum of community size: different from total graph size because
    # of overlaps.
    for comm in part1:     # for each community
        size += len(comm)  # calculate total size

or

    # Max jaccard score, return (J, best_name) tuple
    max( (J(a, b), bname) for bname, b in b.iteritems() )

Use organization, function names, docstrings, and comments, and
spacing all to make your code easily understandable.


### Naming conventions

Spend a few seconds to come up with good names before writing, it will
go a long way to making code organized later.  One of the most clever
pieces of advice I have read was "program with an open thesaurus" to
find the perfect names for modules, functions, variables, and so on.
Another thing people say is to think up a metaphor for what you are
programming, and use that for picking good names.


Modules:

* Short lowercase, like `plots.py`

Class names:

* CapWords, such as `MyClass`

Function names:

* lowercase words separated with underscores, e.g., `my_great_func()`

Constants (on the module level):

* All capital letters with underscores, e.g., `VALUE_OF_PI`

### Documenting code

See notes about Sphinx documentation below.

#### Docstrings:

Use one-liners for really obvious cases:

```python
def sum(x, y):
	""" Return the sum of x and y."""
	return x + y
```

Multiline docstrings:
* First one line summary, then more elaborate description, if necessary.

```python
def complex(real=0.0, imag=0.0):
    """Form a complex number.

    Keyword arguments:
    real -- the real part (default 0.0)
    imag -- the imaginary part (default 0.0)

    Returns:
    complex -- a complex number
    """
    if imag == 0.0 and real == 0.0:
        return complex_zero
    ...
```

Here is a quick checklist to help you write a docstring quickly.  Not
everything listed is needed for every function, just use this as an
idea list, not a checklist.  Start with the most important things, and
if the code becomes popular others will improve it as they review and
use the code.

 * One line imperative sentence describing what the function does.
   "Find the best match community", "Add legend to plot", "Calculate
   burstiness of time sequence array" (most important)
 * Any assumptions made.
 * Meaning of each input argument
 * Output arguments.
 * Time and memory complexity, if nontrivial.


####``verkko``'s Sphinx documentation
``verkko`` uses [Sphinx](http://sphinx-doc.org/) for generating the html
 documentation.
 Especially our Sphinx build support the ``numpydoc`` format for docstrings.
 An example numpydoc, which shows you most of the features:
 [numpydoc example](https://github.com/numpy/numpy/blob/master/doc/example.py)
 or e.g. ``verkko.ptests.permute``

 Otherwise Sphinx uses reStructuredText to compile the documentation.
 See the [Sphinx reST-intro](http://sphinx-doc.org/rest.html) to get going.

 [More info on Shpinx rest references.](http://sphinx-doc.org/domains.html#signatures)

 To compile the documentation locally, use ``make docs`` (at verkko/)
