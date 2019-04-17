# aeffp

Andy's 815 Final Project: Using lasso-penalized logistic regression to explore the influence of local nucleotide sequence context on mutation rate heterogeneity.

## Preparing Data

1. Assuming variant calls are available in a vcf/bcf file, run a command similar to this to extract information about the singletons (note: not all the fields here need to be extracte, but following steps might assume these fields are present in the output)

```bash
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%DP\t%NS\n" /path/to/vcf >> /desired/output/path.txt
```

This should result in a file that looks like the following:

```
10      60524   G       T       34140   3764
10      60606   T       C       10843   3473
10      60765   G       A       29355   3753
10      60776   A       G       33805   3757
10      60828   A       G       47986   3761
```

2. Each site in the above output needs to be annotated with the local nucleotide context. This requires a reference genome and a script to pull the sequence from the reference for each singleton (an example can be found in the `inst/perl` directory -> `annotate.pl`)

3. Generate a count of the number of times each motif occurs in the reference genome. A script that does this is available in the `inst/python` directory (`motif_count.py`). This script requires a text file that lists the motifs you want to count (one per line); an additional script is included (`gen_motif.py`) to generate this file for you.

```bash
python gen_motif.py 
python motif_count.py -i /refGenome.fasta -m motifs5.txt -o output/path
```

You also need to match each motif with its reverse-complement and sum the two together.