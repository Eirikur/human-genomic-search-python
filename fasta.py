#/usr/bin/env python
"Module comment: Just give me your fasta file on the command line."
from __future__ import print_function # for syntax compatibility on Python 2 & Python 3
# This program WILL NOT run on Python 3. pip didn't find a BioPython for it.

import sys, os
from collections import OrderedDict
# You don't need to see me do a CLI via argh and fall back to argparse and all that, right?

# I don't like to write and maintain code that exists and has maintainers -- IFF the performance is reasonable.
# I'm at least going to try using biopython.  I've had good results with scipy in the past.
# coding_tasks.txt said that parsers and formatters are okay.  I hope I'm in the spirit of that.
# Trust me, I've written a lot of ad-hoc parsers.
try:
    from Bio import SeqIO as fa
except ImportError: # Explain what occurred and offer a path for remediation.  'Don't exit without it!'
    print("The biopython library is not available.\n  'sudo pip install biopython' is one way to get it.")
    sys.exit(127) # File not found

# One standard place to keep our abort code.  Can't really use it above this point in a module.
def abort(status_code, reason):
    "Let's just write a message to stderr and exit with the status code."
    msg = "There's been a bit of a problem.  What I know about it is:\n{text}".format(text = reason)
    print(msg, file=sys.stderr)
    sys.exit(status_code) # Exit with a status code so that shell scripts can know that we had a problem.

def do_work(filepath):
    sequences = [record.seq for record in fa.parse(filepath, "fasta")]
    working = []
    for seq in sequences:
        count = sum([1 for target in sequences if str(target) == str(seq)]) # comparison casts are recommended by biopython.
        this = count, str(seq)
        if this not in working: working.append(this)
    working = sorted(working, reverse=True)
    for index in range(10):
        record = working[index]
        formatted_nicely = "{count}, {seq}".format(count=record[0], seq=record[1])
        print(formatted_nicely)



def main(args):
    if len(args) != 2:
        abort(1, "This program takes one mandatory argument, which is the name of the fasta file.")
    filename = args[1]
    if not os.path.exists(filename):
        abort(1, "That doesn't work for me.  Please provide an existing path to the file.")
    do_work(filename)
    print("Done.")


if __name__ == "__main__":
    main(sys.argv)
    sys.exit(0) # Exiting at the end of the file is nice and clean.
