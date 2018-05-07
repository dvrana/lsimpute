# Generates a fake PED and MAP files for benchmarking
# Will generate two ped / map pairs, both containing nsnp SNPs
# One will have nsam samples, the other nref samples
# All SNPs will be on chromosome 1

from random import choice, randrange, random

nsnp = 10
nref = 10
nsam = 10
refname = "ref"
samname = "sam"

alleles = ['A','C','T','G']

# Generate map
snps = []
snpids = set([''])
bp = 0
mapdist = 0.0
for i in range(nsnp):
  a1 = choice(alleles)
  a2 = choice(alleles)
  while (a1 == a2):
    a2 = choice(alleles)
  bp += randrange(700)
  mapdist += random() / 100.0
  rsid = ''
  while rsid in snpids:
    rsid = 'rs' + str(randrange(1000000000))
  snps.append(['1', rsid, str(mapdist), str(bp) ,a1, a2])

# Generate sample ped
sped = []
for i in range(nsam):
  sped.append([str(i), str(i), '0', '0', '0', '0'])
  for j in range(nsnp):
    sped[i].append(choice(snps[j][4:]))
    sped[i].append(choice(snps[j][4:]))

# Generate reference ped
rped = []
for i in range(nsam):
  rped.append([str(i), str(i), '0', '0', '0', '0'])
  for j in range(nsnp):
    rped[i].append(choice(snps[j][4:]))
    rped[i].append(choice(snps[j][4:]))

# Write map files
with open(refname + '.map','w') as f:
  for line in snps:
    f.write(' '.join(line[:4]) + '\n')

with open(samname + '.map','w') as f:
  for line in snps:
    f.write(' '.join(line[:4]) + '\n')

# Write ped files
with open(refname + '.ped','w') as f:
  for line in rped:
    f.write(' '.join(line) + '\n')

with open(samname + '.ped','w') as f:
  for line in sped:
    f.write(' '.join(line) + '\n')
