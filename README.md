# PrimeGradeBeef
A new BEEF tailored for command line execution and integration with geneiousprime via a wrapper plugin. The script accurately filters putative TargetAID base editor targets down to the candidates that can induce premature stop codons within a coding sequence. 

The geneiousprime workflow should be ran on an extracted CDS annotation representing the gene intended for knockout. 

Alternatively, a Genbank formatted file containing a linear extracted sequence with a CDS annotation, and annotated CRISPR target sites can be fed into the script through the commandline using:

```./PGBv5 -input [inputFileName] -output [outputfilename]```
