# aeffp

Andy's 815 Final Project: Using lasso-penalized logistic regression to explore the influence of local nucleotide sequence context on mutation rate heterogeneity.

## Preparing Data

1. Assuming variant calls are available in a vcf/bcf file, run a command similar to this to extract information about the singletons (note: not all the fields here need to be extracte, but following steps might assume these fields are present in the output)

```bash
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%DP\t%NS\n" /path/to/vcf >> /desired/output/path.txt
```