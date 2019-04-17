from itertools import product

alpha = ['A', 'T', 'C', 'G']
motifs = [a+b+c+d+e+f+g for a,b,c,d,e,f,g in product(alpha, repeat=5)]
with open('motifs5.txt', 'w') as f:
    for item in motifs:
        f.write("%s\n" % item)
