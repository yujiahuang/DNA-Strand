from __future__ import print_function
import math
import random
import sys
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
      strands = self.state

      # helper functions

      # similarity: hamming distance of 2 strands
      def similarity(s1, s2):
        length = len(s1) if len(s1) < len(s2) else len(s2)
        sim = 0
        for i in range(length):
          if s1[i] == s2[i]: sim += 1
        return sim

      # slide similarity: hamming distance of 2 strands, with slide (simple similarity excluded)
      def slide_similarity(s1, s2):
        slide_sim = 0
        slide_range = [ (-1) * (len(s1)-1), len(s2)-1 ]
        for i in range(slide_range[0], slide_range[1]): # slide
          sim = 0
          overlap = 0
          for j in range(len(s1)-1):
            if 0 <= j+i < len(s2):
              overlap += 1
              if s1[j] == s2[j+i]:
                sim += 1

          if overlap != 0 and i != 0:  # don't count simple similarity
            slide_sim += sim / overlap  # weighted by overlap: the more overlap the more important
        return slide_sim

      # weighted similarity: the closer to 5', the more weight is added
      def weighted_similarity(s1, s2):
        sim = 0
        length = len(s1) if len(s1) < len(s2) else len(s2)
        for i in range(length):
          if s1[i] == s2[i]: sim += math.pow(0.5, i)
        return sim

      # possible maximum of 3 similarity functions; used in normalization
      def possible_max():

        def calc_max(s1, s2):
          shorter_length = len(s1) if len(s1) < len(s2) else len(s2)
          length_diff = abs(len(s2)-len(s1))
          
          simple_sim = shorter_length
          slide_sim = shorter_length * (length_diff - 1) + shorter_length * shorter_length
          weighted_sim = (1 - math.pow(0.5, shorter_length)) * 2
  
          return [simple_sim, slide_sim, weighted_sim]

        maxes = [0, 0, 0]
        for i in range(len(strands)-1):
          for j in range(i+1, len(strands)):
            s1 = strands[i]['sequence']
            s2 = strands[j]['sequence']
            s2_star = complement(s2)
            maxes = [ x + y for x, y in zip( maxes, calc_max(s1, s2) ) ]
            maxes = [ x + y for x, y in zip( maxes, calc_max(s1, s2_star) ) ]            

        return maxes

      # calculate similarity of each other
      def calc_sim(sim_function):
        sim = 0
        for i in range(len(strands)-1):
          for j in range(i+1, len(strands)):
            s1 = strands[i]['sequence']
            s2 = strands[j]['sequence']
            s2_star = complement(s2)
            sim += sim_function(s1, s2)
            sim += sim_function(s1, s2_star)
        return sim

      # normalization
      def normalize(arr, max):
        # normalize to 0 ~ ceiling and sum up
        normalized = []
        ceiling = 100
        for i, e in enumerate(arr):
          normalized.append( ceiling * e / max[i] )
        return sum(normalized)

      # costs

      cost = []
      # Cost1: similarity of each other
      cost.append( calc_sim(similarity) )
      # Cost2: slide similarity
      cost.append( calc_sim(slide_similarity) )
      # Cost3: weighted similarity (TODO)
      cost.append( calc_sim(weighted_similarity) )

      return normalize(cost, possible_max())


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

# print the DNA strands specifically for amplifier
def print_amplifier(a, b, c):
  hairpin1 = a + b + complement(c) + complement(b)
  hairpin2 =   complement(b) + complement(a) + b + c   
  initiator1 =   complement(b) + complement(a)   
  initiator2 =   complement(c) + complement(b)   
  detect1 =   complement(b) + complement(a) + a + b + complement(c) + complement(b)   
  detect2 =   complement(b) + complement(a) + b + c + complement(c) + complement(b)   

  print('==================================')
  print('a:          ', ''.join(a) )
  print('b:          ', ''.join(b) )
  print('c:          ', ''.join(c) )
  print('hairpin1:   ', ''.join(hairpin1) )
  print('hairpin2:   ', ''.join(hairpin2) )
  print('initiator1: ', ''.join(initiator1) )
  print('initiator2: ', ''.join(initiator2) )
  print('detect1:    ', ''.join(detect1) )
  print('detect2:    ', ''.join(detect2) )
  print('==================================\n')


if __name__ == '__main__':

    # Specify Needs
    strands = [
      { 'name': 'a', 'length': 10 },
      { 'name': 'b', 'length': 16 },
      { 'name': 'c', 'length': 10 },
    ]

    # initial state
    # 1: random
    # 2: nupack
    init_type = 'nupack' if len(sys.argv) > 1 and sys.argv[1] == '-n' else 'random'
    if init_type == 'random':
      print('#################################')
      print('Randomly generate initial strands')
      print('#################################')
      for s in strands:
        init_strand = init_random(s['length'])
        s['sequence'] = init_strand
      print('\ninitial strands:')
      # print(strands)
      print_amplifier(
        strands[0]['sequence'], 
        strands[1]['sequence'], 
        strands[2]['sequence'])
    elif init_type == 'nupack':
      print('##################################')
      print('  Using Nupack as initial strands ')
      print('##################################')
      strands[0]['sequence'] = list('GTTTTTGGTT')
      strands[1]['sequence'] = list('GGGGTGGGGGGGGTGG')
      strands[2]['sequence'] = list('AGGAGAGAGG')
      print('\ninitial strands:')
      print_amplifier(
        strands[0]['sequence'], 
        strands[1]['sequence'], 
        strands[2]['sequence'])


    # anneal
    strand_annealer = GenerateStrand(strands)
    # since our state is just a list, slice is the fastest way to copy
    strand_annealer.copy_strategy = "slice"  
    # auto_schedule = strand_annealer.auto(minutes=1)
    # print(auto_schedule)
    # strand_annealer.set_schedule(auto_schedule)
    strand_annealer.steps = 60000
    # strand_annealer.Tmax = 200
    strand_annealer.Tmin = 1.0
    new_strands, e = strand_annealer.anneal()

    print('\nafter annealing:')
    # print(new_strands)
    print_amplifier(
      new_strands[0]['sequence'], 
      new_strands[1]['sequence'], 
      new_strands[2]['sequence'])

