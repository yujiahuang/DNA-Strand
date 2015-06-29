from __future__ import print_function
import math
import random
from simanneal import Annealer

def complement(strand):
  comp = []
  for n in strand:
    if n == 'A': comp.append('T')
    elif n == 'T': comp.append('A')
    elif n == 'C': comp.append('G')
    elif n == 'G': comp.append('C')
  return comp

class GenerateStrand(Annealer):

    ###
    # Test annealer for generating strands
    ###

    def __init__(self, state):
        super(GenerateStrand, self).__init__(state)  # important! 

    def move(self):
      # randomly choose one strand, and randomly do one of the following moves:
      # 0: flip: flip 5' and 3'
      # 1: switch: change one nucleotide to another type
      # 2: shuffle: cut strand into pieces and randomize the order
      def flip(strand):
        strand.reverse()
      def switch(strand):
        # do we need to check CG content?
        i = random.randint(0, len(strand)-1)
        nucleotides = ['A', 'T', 'C', 'G']
        nucleotides.remove(strand[i])
        n = nucleotides[random.randint(0, 2)]
        strand[i] = n
      def shuffle(strand):
        groups = []
        for i, n in enumerate(strand):
          if i % 3 == 0:
            groups.append([n])
          else:
            groups[i/3].append(n)
        if len(groups[-1]) == 1:
          groups[-2] = groups[-2] + groups[-1]
          groups.pop(-1)

        random.shuffle(groups)
        i = 0
        for g in groups:
          for n in g:
            strand[i] = n
            i += 1

      # choose strand
      num_strands = len(self.state)
      strand = self.state[ random.randint(0, num_strands-1) ]['sequence']

      # choose move
      move = random.randint(0, 2)
      if move == 0:
        flip(strand)
      if move == 1:
        switch(strand)
      if move == 2:
        shuffle(strand)
      # else:
        # impossible

    def energy(self):
      def similarity(s1, s2):
        sim = 0
        length = len(s1) if len(s1) < len(s2) else len(s2)
        for i in range(length):
          if s1[i] == s2[i]: sim += 1
        return sim

      # Cost1: similarity of each other
      strands = self.state
      sim = []
      for i in range(len(strands)-1):
        for j in range(i+1, len(strands)):
          s1 = strands[i]['sequence']
          s2 = strands[j]['sequence']
          s2_star = complement(s2)
          sim.append( similarity(s1, s2) )
          sim.append( similarity(s1, s2_star) )

      # Cost2: slide similarity (TODO)
      # Cost3: weighted similarity (TODO)
          
      return sum(sim)


def init_random(length):
  
  init_strand = []
  
  # define CG content
  CG_content = random.uniform(0.4, 0.6)
  switch_points = [ 
    CG_content * length / 2 , 
    CG_content * length ,
    (CG_content + (1 - CG_content)/2) * length 
  ]

  # generate an ordered string
  for i in range(length):
    nucleotide = 'C'
    if i > switch_points[0]: nucleotide = 'G'
    if i > switch_points[1]: nucleotide = 'A'
    if i > switch_points[2]: nucleotide = 'T'
    init_strand.append(nucleotide)

  # randomize it
  random.shuffle(init_strand)

  return init_strand


if __name__ == '__main__':

    # Specify Needs
    strands = [
      { 'name': 'a', 'length': 10 },
      { 'name': 'b', 'length': 16 },
      { 'name': 'c', 'length': 10 },
    ]

    # initial state
    # 1: random
    for s in strands:
      init_strand = init_random(s['length'])
      s['sequence'] = init_strand
    print('initial strands:')
    print(strands)
    # 2: nupack (TODO)

    strand_annealer = GenerateStrand(strands)
    # since our state is just a list, slice is the fastest way to copy
    strand_annealer.copy_strategy = "slice"  
    new_strands, e = strand_annealer.anneal()

    print('after annealing:')
    print(new_strands)

