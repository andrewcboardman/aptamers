#!/usr/bin/bash
# Without partition function
RNAfold -i samples.fasta --noPS -P ../../../ViennaRNA-2.4.10/misc/dna_mathews2004.par --noconv > samples_struct.fasta
# With partition function
RNAfold -i samples.fasta --noPS -P ../../../ViennaRNA-2.4.10/misc/dna_mathews2004.par --noconv -p > samples_struct.fasta