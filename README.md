# PrimeGradeBeef
I originally made a Base Editor Enrichment Function (BEEF) to simplify the process of finding CRISPR target sites amenable to premature stop codons. That script was pretty unintuitive to use, and required many button clicks and tedious directory navigation to run. It also improperly handled several cases in guide filtering. 

This script, PrimeGradeBeef, is a new BEEF tailored for command line execution and integration with Geneious Prime via a wrapper plugin. The script accurately filters putative TargetAID base editor targets down to the candidates that can induce premature stop codons within a coding sequence. 

The Geneious Prime workflow should be ran on an extracted CDS annotation representing the gene intended for knockout. 

Alternatively, a Genbank formatted file containing a linear extracted sequence with a CDS annotation, and annotated CRISPR target sites can be fed into the script through the commandline using:

```./PGBv5 -input [inputFileName] -output [outputfilename]```

An example input file can be found as 'PrimeGradeBeef_ex.gb.'
