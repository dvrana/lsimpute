# Generates a fake PED and MAP files for benchmarking
# Will generate two ped / map pairs, both containing nsnp SNPs
# One will have nsam samples, the other nref samples
# All SNPs will be on chromosome 1

from random import choice, randrange, random

nsnp = 60000
nref = 2500
nsam = 1
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

print("Generated map")

# Generate sample ped
sped = [[str(i), str(i), '0', '0', '0', '0'] for i in range(nsam)]
for i in range(nsam):
  for j in range(nsnp):
    sped[i].append(snps[j][4 + randrange(2)])
    sped[i].append(snps[j][4 + randrange(2)])

print("Generated sample ped")

# Generate reference ped
rped = [[str(i), str(i), '0', '0', '0', '0'] for i in range(nref)]
for i in range(nref):
  for j in range(nsnp):
    rped[i].append(snps[j][4 + randrange(2)])
    rped[i].append(snps[j][4 + randrange(2)])

print("Generated reference ped")

# Write map files
with open(refname + '.map','w') as f:
  for line in snps:
    f.write(' '.join(line[:4]) + '\n')

with open(samname + '.map','w') as f:
  for line in snps:
    f.write(' '.join(line[:4]) + '\n')

print("Written maps")

# Write ped files
with open(refname + '.ped','w') as f:
  for line in rped:
    f.write(' '.join(line) + '\n')

with open(samname + '.ped','w') as f:
  for line in sped:
    f.write(' '.join(line) + '\n')

print("Written peds")
