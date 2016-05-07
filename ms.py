"""
ms.py

This module is for reading the simulation results from ms.

Sam Vohr (svohr@ucsc.edu)
October 21, 2012
"""

import collections

class MSResult:
    """ Wrapper class for a stream containing the results of a set of
        ms simulations. """
    cmdraw = ''
    nsamples = 0
    nreps = 0
    args = collections.defaultdict(list)
    rseed = ''
    ms_input = None
    curr_sim = None

    def _read_args(self):
        """ Parse the arguments used in the MS invocation into a nice dict """
        items = self.cmdraw.split('-')
        required = items[0].split()
        self.nsamples = int(required[1])
        self.nreps = int(required[2])
        for arg in items[1:]:
            flag, values = arg[0:arg.find(' ')], arg[arg.find(' ')+1:].rstrip()
            self.args[flag].append(values)

    def __init__(self, ms_input):
        """ creates a new MSResult based on the input stream. """
        self.ms_input = ms_input
        # First line is the command with arguments for this simulation.
        self.cmdraw = self.ms_input.readline().rstrip()
        self._read_args()
        self.rseed = self.ms_input.readline().rstrip()
        # Advance the file to the next simulation.
        line = self.ms_input.readline()
        while not line.startswith('//'):
            line = self.ms_input.readline()

    def next_sim(self):
        """ Reads the next simulation from the file. """
        lines = []
        line = self.ms_input.readline()
        if line == '':
            # end of file
            return None
        else:
            while not (line.startswith('//') or line == ''):
                if line != '\n':
                    # if line is not blank, add it to the list for this sim.
                    lines.append(line.rstrip())
                line = self.ms_input.readline()
            # hit the next simulation or end of file
            ms_current_sim = MSSim(self, lines)
            return ms_current_sim

    def pop_sizes( self ):
        """ Reads the population sizes from the ms arguments """
        if 'I' in self.args:
            pop_str = self.args['I'][0]
            items = pop_str.split()
            npops = int(items[0])
            pops = [ int(pop_size) for pop_size in items[1:] ]
            assert npops == len(pops)
            assert self.nsamples == sum(pops)
        else:
            pops = list( [ self.nsamples ] )
        return pops


class MSSim:
    """ Stores the results of a single simulation """
    nsites = 0
    pos = []
    base_pos = None
    base_len = None
    genotypes = []
    ms_parent = None
    def __init__(self, parent, ms_lines):
        self.ms_parent = parent
        self.genotypes = list()
        site_line = ms_lines[0]
        assert site_line.startswith('segsites:')
        self.nsites = int( site_line[site_line.find(' ') + 1:] )
        if self.nsites > 0:
            pos_line = ms_lines[1]
            assert pos_line.startswith('positions:')
            self.pos = pos_line.split()[1:]
            self.genotypes = ms_lines[2:]
        if 'r' in self.ms_parent.args:
            # Grab the locus lenth from the recombination arguments.
            self.base_len = int(self.ms_parent.args['r'][0].split()[1])
            self.base_pos = self.translate_pos_to_base()
        return

    def pop_genotypes(self):
        """ Return the genotypes partitioned by sub-population. """
        pop_size = self.ms_parent.pop_sizes()
        pop_geno = list()
        cur = 0
        for pop in pop_size:
            pop_geno.append( self.genotypes[cur:cur + pop] )
            cur += pop
        return pop_geno

    def translate_pos_to_base(self):
        """ Converts the relative positions into base positions given the
            length of the locus. """
        # if we have more sites than the locus's length, change the length to
        # the number of sites we have. This is an extreme case.
        if len(self.pos) > self.base_len:
            self.base_len = len(self.pos)
        # Make sure there are no double hits here. Shift sites that occur at the
        # same bases over by one.
        bases = [int(float(p) * self.base_len) for p in self.pos]
        base_count = collections.Counter(bases)
        while base_count.most_common(1)[0][1] > 1:
            for bsite in base_count.keys():
                i = 1
                # for each site over 1, move a site over to the next one.
                while base_count[bsite] > 1:
                    base_count[bsite] -= 1
                    base_count[(bsite + i) % self.base_len] += 1
                    i += 1
        base_pos = base_count.keys()
        base_pos.sort()
        return base_pos

