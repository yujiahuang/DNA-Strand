import random

strands = [ 
  ['C', 'G', 'G', 'G', 'T', 'A', 'C', 'C', 'T', 'A'], 
  ['A', 'A', 'C', 'A', 'C', 'T', 'C', 'T', 'G', 'A', 'C', 'G', 'C', 'G', 'T', 'G']
]

def flip(strand):
  strand.reverse()
def switch(strand):
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

# output all
for s in strands:
  print('strand:')
  print(s)
  flip(s)
  print(s)
  switch(s)
  print(s)
  shuffle(s)
  print(s)
  print('\n')
