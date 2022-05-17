#############################################################################
# Functions to read sequences
# Gregg Thomas
#############################################################################

import sys
import os
import gzip
import lib.core as CORE
import multiprocessing as mp
from itertools import groupby

############################################################################# 

def readFasta(filename, seq_compression):
# Read a FASTA formatted sequence file
# Great iterator and groupby code from: https://www.biostars.org/p/710/ 

    if seq_compression == "gz":
        file_stream = gzip.open(filename);
        fa_iter = (x[1] for x in groupby(file_stream, lambda line: line.decode()[0] == ">"));
        readstr = lambda s : s.decode().strip();
    elif seq_compression == "none":
        file_stream = open(filename); 
        fa_iter = (x[1] for x in groupby(file_stream, lambda line: line[0] == ">"));
        readstr = lambda s : s.strip();
    # Read the lines of the file depending on the compression level
    # file_stream opens the file as an iterable
    # groupby takes an iterable (file_stream) and a function that indicates the key of the group. It iterates over
    # the iterable and when it encounters a key, it groups all following items to it (perfect for FASTA files).
    # fa_iter is a generator object holding all the grouped iterators.
    # readstr is a function that changes depending on compression level -- for compressed files we also need to decode
    # each string in the iterators below.

    seqdict = {};
    # A dictionary of sequences:
    # <sequence id/header> : <sequence>

    for header_obj in fa_iter:
        #print(header_obj)
        header = readstr(header_obj.__next__());
        # The header object is an iterator. This gets the string.

        curkey = header[1:];
        # This removes the ">" character from the header string to act as the key in seqdict

        seq = "".join(readstr(s) for s in fa_iter.__next__());
        # The current header should correspond to the current iterator in fa_iter. This gets all those
        # lines and combines them as a string.

        #print(header, len(seq));

        # if curkey in globs['st'].tips:
        seqdict[curkey] = seq;
        # Save the sequence in seqdict if it belongs to a tip branch in the input tree  

    return seqdict;

#############################################################################

def readSeq(seq_dir, globs):

    step = "Reading alignments";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    #step = "Getting FASTA files";
    #step_start_time = PC.report_step(globs, step, False, "In progress...");        
    aln_files = [ os.path.join(seq_dir, f) for f in os.listdir(seq_dir) if f.endswith((".fa", ".fa.gz", ".fasta", ".fasta.gz")) ];
    #step_start_time = PC.report_step(globs, step, step_start_time, "Success: " + str(len(aln_files)) + " FASTA files found");
    ## Read alignment file names
    
    #step = "Detecting compression";
    #step_start_time = PC.report_step(globs, step, False, "In progress...");
    seq_compression = CORE.detectCompression(aln_files[0]);
    # if globs['seq-compression'] == "none":
    #     step_start_time = PC.report_step(globs, step, step_start_time, "Success: No compression detected");
    # else:
    #     step_start_time = PC.report_step(globs, step, step_start_time, "Success: " + globs['seq-compression'] + " detected");
    ## Detect the compression of the input sequence files

    # step = "Reading FASTA files";
    # step_start_time = PC.report_step(globs, step, False, "In progress...");

    aln_dict = {};
    # The dictionary of individual alignments to return:
    # <locus id> : { <sequence id> : <sequence> }

    for f in aln_files:
        locus_id = os.path.splitext(os.path.basename(f))[0];
        # Get the locus ID from the file name

        cur_aln = readFasta(f, seq_compression);
        # Read the current file as a FASTA file
        
        aln_dict[locus_id] = cur_aln;
        # Add the current alignment to the main aln dict

    num_alns = len(aln_dict);
    # step_start_time = PC.report_step(globs, step, step_start_time, "Success: " + str(globs['num-alns']) + " files read");

    # Read sequences if input is a directory of alignment files
    #######################


    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(num_alns) + " alns read");
    # Status update

    return aln_dict, num_alns;

#############################################################################

def subsetAlns(globs, spec_list):
    step = "Subsetting alignments";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Stats update

    subset_alns = {};

    for aln_id in globs['alns']:
        subset_alns[aln_id] = { aln : globs['alns'][aln_id][aln] for aln in globs['alns'][aln_id] if aln_id not in spec_list };
        # Subset the current alignment

        if len(subset_alns[aln_id]) < 4:
            globs['warnings'] += 1;
            CORE.printWrite(globs['logfilename'], -2, "# WARNING: Alignment " + str(aln_id) + " has fewer than 4 sequences after subsetting for pruned species");


    step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
    # Status update

    return globs, subset_alns;

#############################################################################

def locusAlnStats(locus_item):
    locus, aln, skip_chars = locus_item;
    # Unpack the data for the current locus

    num_seqs = len(aln);
    aln_len = len(aln[list(aln.keys())[0]]);

    cur_stats = { 'num-seqs' : num_seqs, 'length' : aln_len, 'avg-nogap-seq-len' : 0, 'variable-sites' : 0, 'unique-seqs' : 0,
                                        'informative-sites' : 0, 'num-sites-w-gap' : 0, 'num-sites-half-gap' : 0,
                                        'num-seqs-half-gap' : 0 };
    # Initialize the stats dict for this locus

    half_aln_len = aln_len / 2;
    half_site_len = num_seqs / 2;
    # Compute half the alignment length and half the site length for the current locus.
    # Used to see if seqs and sites are half-gap or more

    ungapped_seq_lens = [];
    for seq in aln:
        ungapped_seq_lens.append(len(aln[seq].replace("-", "")));
    cur_stats['avg-nogap-seq-len'] = CORE.mean(ungapped_seq_lens);
    # Calculate the average sequence length without gaps for each sequence in the alignment

    for j in range(aln_len):
        site = "";
        for seq in aln:
            if aln[seq][j].count("-") >= half_aln_len:
                cur_stats['num-seqs-half-gap'] += 1;
            # Count whether the current sequence is more than half gap

            site += aln[seq][j];
        # Get each allele from each sequence as the site

        allele_counts = { allele : site.count(allele) for allele in site if allele not in skip_chars };
        # Count the occurrence of each allele in the site

        if len(allele_counts) > 1:
            cur_stats['variable-sites'] += 1;
            # If there is more than one allele in the site, it is variable

            multi_allele_counts = [ allele for allele in allele_counts if allele_counts[allele] >= 2 ];
            # Count the number of alleles present in at least 2 species

            if len(multi_allele_counts) >= 2:
                cur_stats['informative-sites'] += 1;
            # If 2 or more alleles are present in 2 or more species, this site is informative

        if "-" in site:
            cur_stats['num-sites-w-gap'] += 1;

            if site.count("-") >= half_site_len:
                cur_stats['num-sites-half-gap'] += 1;
        # Count whether this site contains a gap and, if so, whether more than half the sequences are a gap
    ## End site loop

    cur_stats['num-unique-seqs'] = len(set(aln.values()));
    # Count the number of unique sequences
    # Could also count number of times each sequence occurs: https://stackoverflow.com/questions/30692655/counting-number-of-strings-in-a-list-with-python#30692666

    # if cur_stats['num-sites-half-gap'] > half_site_len or cur_stats['num-seqs-half-gap'] > half_aln_len:
    #     cur_stats['low-qual'] = True;
    # Setting a flag for low quality sequence to be considered when estimating theta

    return locus, cur_stats;

#############################################################################

def alnStats(globs, alns, proc_pool):
    step = "Calculating alignment stats";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    aln_stats = {};

    with proc_pool as pool:
        for result in pool.imap(locusAlnStats, ((locus, globs['alns'][locus], globs['aln-skip-chars']) for locus in alns)):
        # Loop over every locus in parallel to calculate stats
        # Have to do it this way so it doesn't terminate the pool for sCF calculations

            aln, stats = result;
            aln_stats[aln] = stats;
            # Unpack the current result
            
            # if globs['aln-stats'][aln]['informative-sites'] == 0:
            #     globs['no-inf-sites-loci'].append(aln);
            # If the locus has no informative sites, add to the list here

    # sorted_aln_lens = sorted([ globs['aln-stats'][aln]['length'] for aln in globs['aln-stats'] ]);
    # globs['avg-aln-len'] = CORE.mean(sorted_aln_lens);
    # globs['med-aln-len'] = CORE.median(sorted_aln_lens);
    # Sort alignment lengths and calculate summary stats

    # sorted_avg_seq_lens = sorted([ globs['aln-stats'][aln]['avg-nogap-seq-len'] for aln in globs['aln-stats'] ]);
    # globs['avg-nogap-seq-len'] = CORE.mean(sorted_avg_seq_lens);
    # globs['med-nogap-seq-len'] = CORE.median(sorted_avg_seq_lens);
    # Sort average sequence lengths without gaps and calculate summary statistics

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(len(aln_stats)) + " alignments processed");
    #if globs['no-inf-sites-loci']:
    #    PC.printWrite(globs['logfilename'], globs['log-v'], "# INFO: " + str(len(globs['no-inf-sites-loci'])) + " loci have 0 informative sites and will be removed from the analysis.");
    # Status update

    return aln_stats;

#############################################################################

def writeAlnStats(globs, aln_stats, aln_stat_file):
# Writes alignment stats to output file

    step = "Writing alignment stats";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status updated

    loci_sorted = sorted(aln_stats);
    # Sorts the loci

    with open(aln_stat_file, "w") as outfile:
    # Open the alignment stats file for writing

        first = True;
        # First flace

        for locus in loci_sorted:
        # Write every locus

            if first:
            # For the first locus, extract and write the headers

                keys = list(aln_stats[locus].keys())
                # The headers will be the keys in the locus dict

                headers = ['locus'] + keys;
                # Add 'locus' to the headers for the locus ID

                outfile.write(",".join(headers) + "\n");
                first = False;
                # Write the headers and set the first flag to False
            ###

            outline = [locus] + [ str(aln_stats[locus][key]) for key in keys ];
            outfile.write(",".join(outline) + "\n");
            # Extract the stats for the current locus and write to the file
        ## End locus loop
    ## Close file

    # if globs['no-inf-sites-loci']:
    #     globs['no-inf-loci-file'] = os.path.join(globs['outdir'], globs['no-inf-loci-file']);
    #     with open(globs['no-inf-loci-file'], "w") as outfile:
    #         for locus in globs['no-inf-sites-loci']:
    #             outfile.write(locus + "\n");
    # Write the loci with no informative sites to a file

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: align stats written");
    #globs['aln-stats-written'] = True;
    # Status update

    return globs;

#############################################################################

def writeAlns(globs, alns, aln_dir):
    step = "Writing alignments";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Stats update

    subset_alns = {};
    num_written = 0;

    for aln_id in alns:
        outfilename = os.path.join(aln_dir, aln_id + ".fa");
        with open(outfilename, "w") as outfile:
            for aln in alns[aln_id]:
                outfile.write(">" + aln + "\n");
                outfile.write(alns[aln_id][aln] + "\n");
        num_written += 1;

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(num_written) + " alns written");
    # Status update

    return globs;

#############################################################################