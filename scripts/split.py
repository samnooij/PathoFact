#!/usr/bin/env python

from Bio import SeqIO
import argparse
import os

def batch_iterator(iterator, batch_size):
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.__next__()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Split fasta')
    parser.add_argument('infile',  metavar='input', type=str, help='The input file (FASTA format)')
    parser.add_argument('split_size', metavar='split_size', type=int, help='Size of splitted fasta files')
    parser.add_argument('out', metavar="output_directory", type=str, help='Location of the output directory')
    args = parser.parse_args()

    try:
        os.stat(args.out)
    except:
        os.mkdir(args.out)

    record_iter = SeqIO.parse(open(args.infile),"fasta")
    for i, batch in enumerate(batch_iterator(record_iter, args.split_size)):
        save_to_path = args.out
        filename = "group_%i.fasta" % (i + 1)
        complete_name = os.path.join(save_to_path, filename)
        with open(complete_name, "w") as handle:
            count = SeqIO.write(batch, handle, "fasta")
        print("Wrote %i records to %s" % (count, filename))

