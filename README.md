# Fast FFP
Fast Feature Frequency Profiling

This project contains the code needed to run feature frequency profile (FFP) comparisons of many species. This is useful for analyzing entire phylogenies. This is done through the calculation of pairwise distances between the FFPs of all species combinations. This FFP tool uses Jellyfish k-mer counts to build the FFP of each species.


Commmand Options:

    -f, --folder
Path to the folder containing the fasta sequences of each species. Currently require each fasta to be in a folder with the same name as the fasta file. Ex. /Acanthisitta_chloris/Acanthisitta_chloris.fq
  
    -a, --sequence1
 Full path to sequence 1 when using pairwise option. Pairwise option used if the folder path is not provided.
  
    -b, --sequence2
Full path to sequence 2 when using pairwise option. Pairwise option used if the folder path is not provided.
  
    -l, --length
k-mer length to use in Jellyfish k-mer counting. Default is 20.
  
    -o, --output
Optional path to output file.
  
    -d, --distance-method
Distance method used for FFP comparison. Avalible options: Jenson-Shannon (js), Euclidean(e), Euclidean Squared (e2).
  
    -r, --rerun
Recounts all k-mers with Jellyfish.
  
    -v, --versbose
Prints additional information to console while running. 
    
  
Future improvements:
1. Pre-calculate the FFP of each species. Then reuse the FFPs to decrease the running time.
2. Use jellyfish python bindings to decrease the running time.
3. Added an installation method.
