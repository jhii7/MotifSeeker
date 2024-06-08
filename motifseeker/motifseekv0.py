import argparse
from argparse import RawTextHelpFormatter
import os
import sys
import numpy as np
import pandas as pd
import scipy
from Bio import SeqIO
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 
from scipy.stats import fisher_exact
import math
import random
import timeit

nucs = {"A": 0, "C": 1, "G": 2, "T": 3}

def ERROR(msg):
	"""
	Print an error message and die

	Parameters
	----------
	msg : str
	   Error message to print
	"""
	sys.stderr.write("[ERROR]: " + "{msg}\n".format(msg=msg) )
	sys.exit(1)

def ExtractSequencesFromBed(bed_file, ref_genome_file):
    """
    Extract sequences from a reference genome based on BED file coordinates.

    Parameters
    ----------
    bed_file : str
        Path to the BED file containing coordinates.
    ref_genome_file : str
        Path to the reference genome in FASTA format.

    Returns
    -------
    sequences : list of str
        List of extracted sequences.
    """
    sequences = []

    # Read reference genome
    ref_genome = SeqIO.to_dict(SeqIO.parse(ref_genome_file, "fasta"))

    # Read BED file and extract sequences
    with open(bed_file, "r") as bed:
        for line in bed:
            if line.strip():
                fields = line.strip().split()
                chrom, start, end = fields[0], int(fields[1]) - 1, int(fields[2])
                if chrom in ref_genome:
                    seq_record = ref_genome[chrom]
                    seq = seq_record.seq[start:end]
                    sequences.append(str(seq))
                else:
                    print(f"Chromosome {chrom} not found in reference genome.")

    return sequences

def get_reads(sequences, length):
    """
    Generate reads of a specified length from a list of sequences.

    Parameters
    ----------
    sequences : list of str
        List of nucleotide sequences from which to generate reads.
    length : int
        Length of the reads to generate.

    Returns
    -------
    reads : list of str
        List of reads of the specified length extracted from the sequences.
    """
    # Generates reads of length from sequences, store in array
    reads = []
    for sequence in sequences:
        i = 0
        end = length
        while end < len(sequence) + 1:
            reads.append(sequence[i:end])
            i += 1
            end += 1
    return reads

def ComputeNucFreqs(sequences):
    """ Compute nucleotide frequencies of a list of sequences
    
    Parameters
    ----------
    sequences : list of str
       List of sequences
       
    Returns
    -------
    freqs : list of float
       Frequencies of A, C, G, T in the sequences
    """
    freqs = [0.25, 0.25, 0.25, 0.25]
    count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    total = 0

    # Count each nucleotide
    for seq in sequences:
        for nucleotide in seq:
            if nucleotide in count:
                count[nucleotide] += 1
                total += 1

    # Calculate frequencies
    if total > 0:
        freqs = [count[n] / total for n in ['A', 'C', 'G', 'T']]
    else:
        freqs = [0.25, 0.25, 0.25, 0.25]  # Default to equal distribution if no valid data

    return freqs

def RandomSequence(n, freqs=[0.25, 0.25, 0.25, 0.25]):
    """ Generate a random string of nucleotides of length n
    
    Use the given nucleotide frequences
    
    Parameters
    ----------
    n : int
       Length of random string to generate
    freqs : list of float
       List of frequencies of A, C, G, T
       
    Returns
    -------
    seq : str
       random sequence of length n with the specified allele frequencies
    """
    seq = "A"*n
    nucs = ['A', 'C', 'G', 'T']
    seq = ''.join(random.choices(nucs, weights=freqs, k=n))
    return seq

def GetPFM(sequences):
    """ Compute the PFM for a set of sequences
    
    Parameters
    ----------
    sequences : list of str
        List of sequences (e.g. binding_sites)
    
    Returns
    -------
        pfm : 2d np.array
    """
    nucs = {"A": 0, "C": 1, "G": 2, "T": 3}
    pfms = []
    for seq in sequences:
        pfm = np.zeros((4, len(seq)))
        rowcount = 0
        for nucleotide in seq:
            pfm[nucs[nucleotide]][rowcount] += 1
            rowcount += 1
        pfms.append(pfm)
    return pfms

def GetPWM(binding_sites, background_freqs=[0.25, 0.25, 0.25, 0.25]):
    """ Compute the PWM for a set of binding sites
    
    Parameters
    ----------
    binding_sites : list of str
        List of sequences 
    background_freqs: list of float
        Background frequency of A, C, G, T
    
    Returns
    -------
        pwm : 2d np.array
        
    Assumes all sequences have the same length
    """
    pwms = []
    nucs = {"A": 0, "C": 1, "G": 2, "T": 3}
    
    for seq in sequences:
        length = len(seq)
        pfm = np.zeros((4, length))
        for idx, nucleotide in enumerate(seq):
            if nucleotide in nucs:
                pfm[nucs[nucleotide], idx] += 1
        
        pfm += 0.01
        
        col_sums = np.sum(pfm, axis=0)
        pwm = np.zeros_like(pfm)
        
        for i in range(4):
            pwm[i, :] = np.log2((pfm[i, :] / col_sums) / background_freqs[i])
        
        pwms.append(pwm)
    
    return pwms

def ScoreSeq(pwm, sequence):
    """ Score a sequence using a PWM
    
    Parameters
    ----------
    pwm : 2d np.array
       Position weight matrix
    sequence : str
       Sequence of nucleotides to be scored
       
    Returns
    -------
    score : float
       PWM score of the sequence
    """
    score = 0
    for i, nucleotide in enumerate(sequence):
        row_index = nucs[nucleotide]
        score += pwm[row_index, i]

    return score

def GetThreshold(null_dist, pval):
    """ Find the threshold to achieve a desired p-value
    
    Given a null distribution (list of values),
    find the threshold to achieve a desired p-value
    
    Parameters
    ----------
    null_dist : list of float
       Null distribution of scores
    pval : float
       % of null_dist that should be above the threshold returned
       
    Returns
    -------
    thresh : float
       Threshold to achieve the desired p-value    
    """
    thresh = 0 
    thresh = np.percentile(null_dist, 100 * (1 - pval))
    return thresh

def ParseMotifsFile(filename):
    """
    Parse a file containing motif sequences in HOMER-style format.

    Parameters
    ----------
    filename : str
        Path to the file containing motifs.

    Returns
    -------
    motifs : list of str
        List of motif sequences parsed from the file.
    """
    # Outputs a list of lists of motif sequences from HOMER-style motif databases
    motifs = []
    with open(filename, "r") as f:
        for line in f:
            if line.startswith(">"):
                motifs.append(line.split()[0].strip(">"))
    return motifs
  
def FindExactMatches(sequences, motifs):
    """
    Find exact matches for motifs in a list of sequences.

    Parameters
    ----------
    sequences : list of str
        List of nucleotide sequences to search for motifs.
    motifs : list of str
        List of motif sequences to search for in the sequences.

    Returns
    -------
    found_motifs : list of str
        List of motifs found in the sequences.
    """
    # Find exact matches for motifs in database
    found_motifs = []
    for motif in motifs:
        for sequence in sequences:
            # Check if motif is present in the sequence using string matching
            if motif in sequence:
                found_motifs.append(motif)
    return found_motifs

def ReverseComplement(sequence):
    """ Get the reverse complement of a sequence
    
    Parameters
    ----------
    sequence : str
      Sequence of nucleotides
      
    Returns
    -------
    revcomp : str
      Reverse complement of sequence
    """
    revcomp = ""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    revcomp = ''.join(complement[nt] for nt in reversed(sequence))
    return revcomp

def FindMaxScore(pwm, sequence):
    """ Get highest PWM match for a sequence
    
    Scan a sequence with a pwm
    Compute the highest pwm score for a given sequence
    Be sure to check the forward and reverse strands!
    
    Parameters
    ----------
    pwm : 2d np.array
       PWM matrix
    sequence : str
       Sequence of nucleotides
       
    Returns
    -------
    max_score : float
       Score of top match to the PWM
    """
    max_score = -1*np.inf
    n = pwm.shape[1]
    max_score = -np.inf

    def scan_seq(seq):
        max_seq_score = -np.inf
        for i in range(len(seq) - n + 1):
            subseq = seq[i:i+n]
            score = ScoreSeq(pwm, subseq)
            max_seq_score = max(max_seq_score, score)
        return max_seq_score

    max_score = max(scan_seq(sequence), scan_seq(ReverseComplement(sequence)))

    return max_score

def GetThreshold(null_dist, pval):
    """
    Given a null distribution (list of values),
    find the threshold to achieve a desired p-value
    
    Parameters
    ----------
    null_dist : list of float
       Null distribution of scores
    pval : float
       pval% of null_dist should be above the threshold returned
       
    Returns
    -------
    thresh : float
       Threshold to achieve the desired p-value    
    """
    thresh = 0
    null_dist.sort()
    index = int((1 - pval) * len(null_dist))
    thresh = null_dist[index - 1] if index > 0 else null_dist[0]
    return thresh

def ComputeEnrichment(peak_total, peak_motif, bg_total, bg_motif):
    """ Compute fisher exact test to test whether motif enriched in bound sequences
    
    Parameters
    ----------
    peak_total : int
       Number of total peaks
    peak_motif : int
       Number of peaks matching the motif
    bg_total : int
       Number of background sequences
    bg_motif : int
       Number of background sequences matching the motif
       
    Returns
    -------
    pval : float
       Fisher Exact Test p-value    
    """
    table = [[peak_motif, peak_total - peak_motif],
             [bg_motif, bg_total - bg_motif]]

    _, pval = scipy.stats.fisher_exact(table)
    return pval

def MotifEnrichmentAnalysis(peak_seqs, pwm_list, pwm_thresholds, background_freqs, pval_threshold=0.05):
    """
    Perform motif enrichment analysis using sequences extracted from BED file.

    Parameters
    ----------
    bed_file : str
        Path to the BED file containing peak coordinates.
    ref_genome_file : str
        Path to the reference genome in FASTA format.
    pwm_list : list of np.array
        List of PWMs to be tested.
    pwm_thresholds : list of float
        List of thresholds corresponding to each PWM.
    background_freqs : list of float
        Background frequencies of A, C, G, T.

    Returns
    -------
    results : list of tuple
        List of tuples containing PWM index, threshold, number of peaks passing,
        number of background sequences passing, and p-value.
    """

    if (args.pval is not None) :
        pval_threshold = args.pval
    
    bg_seqs = [RandomSequence(len(peak_seqs[0]), background_freqs) for _ in range(len(peak_seqs))]

    results = []

    for i in range(len(pwm_list)):
        pwm = pwm_list[i]
        thresh = pwm_thresholds[i]

        num_peak_pass = np.sum([int(FindMaxScore(pwm, seq) > thresh) for seq in peak_seqs])
        num_bg_pass = np.sum([int(FindMaxScore(pwm, seq) > thresh) for seq in bg_seqs])

        pval = ComputeEnrichment(len(peak_seqs), num_peak_pass, len(bg_seqs), num_bg_pass)
        enriched = pval < pval_threshold

        results.append((i, motifs[i], thresh, num_peak_pass, num_bg_pass, pval, enriched))

    return peak_seqs, bg_seqs, results

def PrintResults(peak_seqs, bg_seqs, results):
    """
    Print the results of motif enrichment analysis.

    Parameters
    ----------
    peak_seqs : list of str
        List of sequences extracted from the peaks.
    bg_seqs : list of str
        List of background sequences.
    results : list of tuple
        List of tuples containing PWM index, threshold, number of peaks passing,
        number of background sequences passing, p-value, and enrichment status.
    """
    for result in results:
        pwm_index, motif_seq, thresh, num_peak_pass, num_bg_pass, pval, enriched = result
        enriched_str = "yes" if enriched else "no"
        print(f"{pwm_index}\t{motif_seq}\t{num_peak_pass}/{len(peak_seqs)}\t{num_bg_pass}/{len(bg_seqs)}\t{pval:.3e}\t{enriched_str}")


    
##################################

# Create parser
parser = argparse.ArgumentParser(prog="MotifSeeker", description= "Command-line tool to perform Motif Enrichment Analysis", formatter_class=RawTextHelpFormatter)

# Input file, we want it to work on .bed files #
parser.add_argument("inputfile", help=".bed file containing peaks", metavar ="FILE", type=str)

parser.add_argument("genome", help="fasta reference genome file", \
                        metavar="FILE", type=str)

# Other arguments
parser.add_argument("-o", "--out", help="Write output to file." \
                        "Default: stdout", metavar="FILE", type=str, required=False)
    
parser.add_argument("-p", "--pval", help="p value threshold for motif sequences" \
                         "Default: 0.05", metavar="NUMERICAL P VALUE", type=float, required=False)
    
parser.add_argument("-s", "--size", help="size of motif" \
                        "Default: 100", metavar="NUMERICAL SIZE VALUE", type=int, required=False)

# Parse arguments here
args = parser.parse_args()

if args.size is not None:
        size = args.size
else:
        size = 100
        
# Setup output file
if args.out is None:
        outf = sys.stdout
else: 
        outf = open(args.out, "w")


# Do something with input and genome file if they exist.
if ((args.inputfile is not None) and (args.genome is not None)):
       
        # Begin timer, use for runtime comparison against HOMER
        start = timeit.default_timer()



        # Get sequences, store in var sequences
        sequences = ExtractSequencesFromBed(args.inputfile, args.genome)
       

        # Get reads from sequences lengths 8 to 25
        reads = []
        i = 8 # Minimum motif length from HOMER database
        while i < 25: # Maximum motif length from HOMER database
            reads.append(get_reads(sequences, i))
            i += 1

        # GetPWM already runs pfms, don't need to run GetPFM twice.
        pwms = GetPWM(reads)


        # Store HOMER motifs in var motifs
        motifs = ParseMotifsFile("../motifs/custom.motifs")
        null_dist = null_dist = [
        0.1, 0.2, 0.15, 0.3, 0.25, 0.35, 0.4, 0.45, 0.5, 0.55,
        0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05,
        1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55,
        1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.05,
        2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55
        ]

        thresholds = [GetThreshold(null_dist, 0.05) for _ in pwms]

        # Find exact matches (TESTING)
        # print(FindExactMatches(reads, motifs))

        peak_seqs, bg_seqs, results = MotifEnrichmentAnalysis(
        peak_seqs=sequences,
        pwm_list=pwms,
        pwm_thresholds=thresholds,
        background_freqs=[0.25, 0.25, 0.25, 0.25])

        PrintResults(peak_seqs, bg_seqs, results)

        # End timer
        stop = timeit.default_timer()
        print("Time: ", stop - start, " seconds")