#/usr/bin/env python
"Module comment: Just give me your position and GTF files on the command line."

# STRATEGY STATEMENT:
# I'm going to organize both the work list (positions), and the database (the GTF file) by chromosome.
# By doing that I stand a reasonable chance of having the 'current chromosome' in
# the L3/L2 processor caches.  The chromosome is represented by a dict in both places
# (the work list and the database).  It's those dicts that we want to have cached, if possible.
# I really don't have time to fully performance tune this, but I'd love to do it.
# Recriminations: Maybe if I thought in database-ese I'd have used sqlite here.  I could have gotten
# some algorithms for free there.  Again, this needs actual measurement before warping it to fit
# some performance hypothesis.

# This file is probably difficult to read in spots despite my care for clarity.  I use
# rainbow-delimiters-mode for colorization of nested [({<>})] and I wouldn't want to be without it.

# Your term 'annotation' doesn't match anything in the GTF (GFF V2) file format specs that I can find.
# I'm guessing that you want the tag/values from the attribute field.  There are no records with data in the comment field.
# I didn't send in a question about this because if I can find the right record, it really doesn't matter what field.

from __future__ import print_function # for syntax compatibility on Python 2 & Python 3

# This particular program DOES actually run correctly on Python3.  Tested on Python 2.7.4 and 3.3.1

import sys, os
from bisect import bisect_left
from timeit import timeit

from dotdict import DotDict # My own little Python language extension.  OO notation for dicts.
# Using dotdict greatly clarified the use of named fields in the GTF records.
# I didn't use it on the first pass and the nested brackets became very hard to read.

# One standard place to keep our abort code.  Can't really use it above this point in a module.
def abort(status_code, reason):
    "Let's just write a message to stderr and exit with the status code."
    msg = "There's been a bit of a problem.  What I know about it is:\n{text}".format(text = reason)
    print(msg, file=sys.stderr)
    sys.exit(status_code) # Exit with a status code so that shell scripts can know that we had a problem.

def do_work(filenames):
    positions_filename = filenames[0]; gtf_filename = filenames[1] # I depart from PEP-8 standard format only to trivially save vertical space.
    # I love method-chaining.
    position_data =  open(positions_filename, 'r').readlines() # Reference to file object is not kept; GC can take it; we are done.
    position_data = [record.strip().split('\t') for record in position_data] # 2X memory usage, but briefly.
    positions = {record[0]: [] for record in position_data if len(record) > 1}# Each chromosome owns a list of positions
    [positions[record[0]].append(int(record[1])) for record in position_data if len(record) > 1] # Rare non-use of explict result.
    # [print(position, positions[position]) for position in positions] # This actually helped in debugging.
    # The GTF format is unlikely to change much, but this is a good generic approach to records.
    gtf_fields = ['seqname', # Field names from the GTF (GFF V2) spec.
                  'source',
                  'feature',
                  'start',
                  'end',
                  'score',
                  'strand',
                  'frame',
                  'attribute',
                  'comment']
    # Dictionary comprehension.  This is a new feature in Python 2.7/3.x
    gtf = {key: index for (index, key) in enumerate(gtf_fields)}
    gtf = DotDict(gtf) # Enable OO dot-notation for field names. My own feature.
    gtf_records = open(gtf_filename, 'r').readlines()
    gtf_records = [record.strip().split('\t') for record in gtf_records] # 2X Memory, briefly.
    chromosomes = {record[gtf.seqname]: {} for record in gtf_records} # Each Ch is a dict in a master dict.
    def populate(record):
        chromosomes[record[gtf.seqname]] [int(record[gtf.start])] = [int(record[gtf.end]),
                                                                      record[gtf.attribute]]

    # The below was a very close horserace even at 100 iterations (15 minutes total execution time)
    #
    # def test1():
    #     [populate(record) for record in gtf_records] # Speed?  timeit()
    # def test2():
    #     map(populate, gtf_records) # Compare speed with the above implementation
    # def test3():
    #     [chromosomes[record[gtf.seqname]].update({int(record[gtf.start]) : [int(record[gtf.end]),
    #                                                                   record[gtf['attribute']]
    #                                                                  ]}) for record in gtf_records] # I think this will be faster
    # t1=timeit(test1, number=100)
    # t2=timeit(test2, number=100)
    # t3=timeit(test3, number=100)
    # print(t1, t2, t3)

    [populate(record) for record in gtf_records] # This was the fastest by a few milliseconds     # map(populate, gtf_records) # Compare speed with the above implementation

    # Sort the starting coordinates to allow later binary search.
    # I think this is my major insight for optimization here and I don't think
    # it is a premature optimization.  Any speed gain happens in the inner loop
    # at the cost of this setup.
    for ch in chromosomes:
        chromosomes[ch]['sorted_starts'] = sorted(list(set(chromosomes[ch]))) # TODO: Is this really the place to do de-duplication?

    for chromosome in positions: # It's time to get to the actual work now.
        target_coordinates = positions[chromosome] # List of coordinates on chromosome, from the input.
        target_coordinates = list(set(target_coordinates)) # TODO: Is this really the place to do de-duplication?
        if chromosome in chromosomes:  # There's a rogue in the input that's not in the GTF.
            target_starts = chromosomes[chromosome]['sorted_starts'] # GTF data, seq. starts are keys in a dict.
            for coordinate in target_coordinates:
                index = bisect_left(target_starts, coordinate) # Binary search from the standard library.  Fast.
                closest_sequence_start = target_starts[index]
                sequence_end = chromosomes[chromosome][closest_sequence_start][0] # First element is end of sequence.
                if closest_sequence_start <= coordinate <= sequence_end: # We have a winner!  Vanna will show you your prize after this message.
                    # I'm getting tired at this point.  Forgive the constant indeces, below.
                    gene_name = chromosomes[chromosome][closest_sequence_start][1].split(';')[0].split()[1].replace('"', '')
                    output = "{coordinate}\tfound in {gene} on {chromosome}".format(gene=gene_name, coordinate=coordinate,chromosome=chromosome)
                    print(output) # Output to sdtout, redirect or pipe to where you want it.

        else:
            print("<stderr> Hey! Chromosome {} has positions in the input, but it does not occur in the GTF.".format(chromosome),
                  file=sys.stderr)
        for coord in target_coordinates:
            pass

def main(args):
    if len(args) != 3:
        abort(1, "This program takes two mandatory arguments, a position file and a GTF file")
    filenames = args[1:]
    for filename in filenames:
        if not os.path.exists(filename):
            abort(1, "That doesn't work for me.  Please provide existing, acessible, files.")
    do_work(filenames)

if __name__ == "__main__":
    main(sys.argv)
    sys.exit(0) # Exiting at the end of the file is nice and clean.
