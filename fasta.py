"""
    fasta.py

    This module contains simple functions for parsing FASTA input. 

    Sam Vohr (svohr@soe.ucsc.edu)
"""

import re
import textwrap

class FastaError(Exception):
    """
    Exception for reporting specific problems reading FASTA input.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def read_seq ( fasta_in ):
    """ Generator for reading fasta files one sequence at a time. """
    
    seq_id = None   # The name of the current sequence.
    comment = None  # The comment data for the current sequence.
    seq = None      # what we have read for the sequence so far.
                
    for line in fasta_in:
        
        line = line.rstrip('\r\n')

        if line.startswith('>'):
            # Id line
            if seq_id is not None:
                # return the previous sequence
                yield (seq_id, comment, seq)
            
            match = re.match(r'^>([^\s,]*),?\s*(.*)$', line)
            if match is not None:
                seq_id  = match.group(1)
                comment = match.group(2)
                seq = ""
            else:
                raise FastaError('Bad ID line.')
        else:
            # Sequence line
            if seq_id is not None:
                # add the sequence data and remove the whitespace.
                seq = seq + ''.join(line.split())
            else:
                raise FastaError('Sequence data before ID line.')
    if seq_id is not None:
        yield (seq_id, comment, seq)

# Write a single sequence to a file.
def write_seq( fasta_out, name, comment, sequence ):
    """
    Write a single fasta entry (name, comment, sequence) to the output stream.
    Wraps sequences at 60 chars for readability.
    """
    fasta_out.write( '>%s, %s\n' % (name, comment) )
    for line in textwrap.wrap( sequence, 60 ):
        fasta_out.write( '%s\n' % line )


