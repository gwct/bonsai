#############################################################################
# Calculates concordance factors
# Gregg Thomas
#############################################################################

import sys
import time
import copy
import math
import lib.core as CORE
import lib.tree as TREE

#############################################################################

def gcf(st, gts):
# Calculates gene concordance factor for every node in the given species tree based on quartet
# presence in the gene trees

    gcf_counts = { node : { "concordant" : 0.0, "disco1" : 0.0, "disco2" : 0.0, "para" : 0.0, "decisive" : 0.0, "nondecisive" : 0.0 } for node in st.internals[:-1] };
    # Counts of the concordant/discordant topologies of each node in the species tree

    total_trees = 0.0;
    # Total gene trees coutned

    st_quartets = st.getQuartets(root=False);
    # Get the quartets around each node in the species tree

    for gt_id in gts:
    # Count topologies in every gene tree

        if gts[gt_id] == "NA":
            continue;

        gt_clades = gts[gt_id].getClades();
        gt_splits = gts[gt_id].getSplits();
        gt_tips = set(gts[gt_id].tips);
        # Get the clades, splits, and tips of the current gene tree

        cur_st_quartets = { node : {} for node in st_quartets };
        # The species tree quartets for the current gene tree, which may be different
        # from the species tree because of missing taxa in the gene tree

        for node in st_quartets:
        # Check every species tree node in the current gene tree

            for q in ["d1", "d2", "s", "q4"]:
                cur_st_quartets[node][q] = st_quartets[node][q] & gt_tips;
            # For each clade in the quartet of the curretn node, get the intersect between it and
            # the tips of the current gene tree

            if not all(cur_st_quartets[node][q] for q in ["d1", "d2", "s", "q4"]):
                gcf_counts[node]["nondecisive"] += 1;
                continue;
            # Make sure the gene tree is decisive for this node -- all current quartet clades are defined
            # If not, skip

            possible_topos = { "concordant" : cur_st_quartets[node]["d1"] | cur_st_quartets[node]["d2"],
                                "disco1" : cur_st_quartets[node]["d1"] | cur_st_quartets[node]["s"],
                                "disco2" : cur_st_quartets[node]["d1"] | cur_st_quartets[node]["q4"] };          
            # The poosible topologies for this node in the species tree, given the current quartet
            # The one defined in the species tree (concordant), and the two nearest neighbor interchanges (disco1, disco2)

            if gts[gt_id].findSplits(possible_topos["concordant"], gt_clades, gt_splits):
                gcf_counts[node]["concordant"] += 1;
            elif gts[gt_id].findSplits(possible_topos["disco1"], gt_clades, gt_splits):
                gcf_counts[node]["disco1"] += 1;
            elif gts[gt_id].findSplits(possible_topos["disco2"], gt_clades, gt_splits):
                gcf_counts[node]["disco2"] += 1;
            else:
                gcf_counts[node]["para"] += 1;
            # Check which topology for this quartet is found in the species tree            

            gcf_counts[node]["decisive"] += 1;
            # Add to the count of decisive gene trees
        ## End st node loop

        total_trees += 1.0;
        # Iterate the total tree count
    ## End gene tree loop

    return gcf_counts, total_trees;

#############################################################################

def countTopos(globs, st, gts):

    step = "Counting topologies";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update  

    num_trees, num_mismatch = 0, 0;

    #st_topo = set([ frozenset(st.getClade(node)) for node in st.internals ]);
    # For rooted trees?
    st_topo = set([ frozenset({frozenset(st.getClade(node)), frozenset(st.getSplit(node))}) for node in st.internals ]);
    # A tree is defined by splits at each branch, so we retrieve them for the species tree here to calculate RF

    topo_counts = { "topo" : [], "count" : [], "tree" : [], "rf" : [], "splits" : [] };
    # The topology count lists

    for gt_id in gts:
    # Count the topology of every gene tree
    
        if gts[gt_id] == "NA":
            continue;

        num_trees += 1;
        if set(st.tips) != set(gts[gt_id].tips):
            CORE.printWrite(globs['logfilename'], -2, "# WARNING: Gene tree on line " + str(gt_id) + " of gene tree file (-g) is missing taxa relative to the reference tree! Skipping topology counting.");
            num_mismatch += 1;
            continue;

        #gt_topo = set([ frozenset(gts[gt_id].getClade(node)) for node in gts[gt_id].internals ]);
        # For rooted trees?
        gt_topo = set([ frozenset({frozenset(gts[gt_id].getClade(node)), frozenset(gts[gt_id].getSplit(node))}) for node in gts[gt_id].internals ]);
        # Get the splits that define each branch in the gene tree

        if gt_topo in topo_counts["topo"]:
        # If the gene tree splits set has been found already, increment the count
            topo_ind = topo_counts["topo"].index(gt_topo);
            topo_counts["count"][topo_ind] += 1;
        
        else:
        # If this is the first time seeing this gene tree, add it and its info to the lists
            topo_counts["topo"].append(gt_topo);
            topo_counts["count"].append(1);
            topo_counts["tree"].append(gts[gt_id].topo_str);

            num_diffs = len((gt_topo - st_topo)) + len((st_topo - gt_topo));
            topo_counts["rf"].append(num_diffs);
            topo_counts["splits"].append(len(st_topo) + len(gt_topo));
            # The number of differences between the set of splits is the Robinson-Foulds distance
            # For normalization, save the total number of splits as well
    ## End gene tree loop

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: " + str(len(topo_counts["topo"])) + " unique topologies found");
    if num_mismatch == num_trees:
        CORE.printWrite(globs['logfilename'], globs['log-v'], "# WARNING: All gene trees are missing taxa relative to the reference tree (from -st). No topologies counted.");
    elif num_mismatch:
        CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: " + str(num_mismatch) + " gene trees are missing taxa relative to the reference tree (from -st) and were skipped during topology counting. See log for more info.");
    # Status updates

    return topo_counts, num_trees, num_mismatch;

#############################################################################

def locusSCF(locus_item):
    locus, aln, quartets, st, aln_skip_chars = locus_item
    # Unpack the data for the current locus

    aln_len = len(aln[list(aln.keys())[0]]);
    # Get the alignment length from the first sequence

    locus_scf = {};
    # The dictionary to calculate average sCF across all nodes in the current locus
    # <locus id> : { <node id> : <scf sum>, <num quartets>, <avg scf> }

    quartet_scores = {};
    # Dictionary to hold all counts for each quartet sampled

    # node_scf = [];
    # The list of sCFs in this locus

    for node in st.nodes:
    # Calculate sCF for every eligible node in the tree

        if node not in quartets:
            continue;
        # Cannot calculate sCF for tips, the root, or node descendant from the root

        # locus_scf[node] = { 'scf-sum' : 0, 'num-quartets' : 0, 'avg-scf' : "NA" };
        # Initialize the scf dict for the current node

        quartet_scores[node] = {};
        # Dictionary to hold all counts for each quartet sampled

        for quartet in quartets[node]:
        # quartet is a tuple of tuples:
        # ((split1-spec1, split1-spec2), (split2-spec1, split2-spec2))

            quartet_scores[node][quartet] = { 'variable-sites' : 0, 'decisive-sites' : 0, 'concordant-sites' : 0, 'disco1-sites' : 0, 'disco2-sites' : 0, 'scf' : "NA" };
            # Initialize the current scf dict for the current quartet

            split1 = quartet[0];
            split2 = quartet[1];
            quartet_spec_list = split1 + split2;
            # Parse out the species for easier access

            if not all(spec in list(aln.keys()) for spec in quartet_spec_list):
                continue;
            # Skip any alignments that don't have all 4 species from the current quartet
            # Not sure if this is the best way to deal with missing data...

            #sub_aln = { spec : aln[spec] for spec in quartet_spec_list };
            sub_aln = [ aln[spec] for spec in quartet_spec_list ];
            # Get the sub-alignment for the current quartet

            sites = zip(*sub_aln);

            for site in sites:
                if any(char in site for char in aln_skip_chars):
                    continue;
                # We only care about sites with full information for the quartet, so skip others here

                uniq_site = set(site);
                # The alleles at the current site

                if len(uniq_site) > 1:
                    quartet_scores[node][quartet]['variable-sites'] += 1;
                # If there is more than one allele, the site is variable                

                    if all(site.count(allele) == 2 for allele in uniq_site):
                    #if len(uniq_site) == 2:
                        quartet_scores[node][quartet]['decisive-sites'] += 1;
                    # If all alleles have a count of 2 at this site, it is decisive

                        if site[0] == site[1] and site[2] == site[3]:# and split1_site[0] != split2_site[0]:
                        #if split1_site.count(split1_site[0]) == len(split1_site) and split2_site.count(split2_site[0]) == len(split2_site):
                            quartet_scores[node][quartet]['concordant-sites'] += 1;
                        # If the alleles from split1 match and the alleles from split2 match, the site is concordant
                        elif site[0] == site[2] and site[1] == site[3]:
                            quartet_scores[node][quartet]['disco1-sites'] += 1;
                        elif site[0] == site[3] and site[1] == site[2]:
                            quartet_scores[node][quartet]['disco2-sites'] += 1;

            # for j in range(aln_len):
            # # Count every site in the sub-alignment

            #     site, split1_site, split2_site = "", "", "";
            #     # The full site, the site for split1 and the site for split2

            #     for seq in sub_aln:
            #         site += sub_aln[seq][j];
            #         if seq in split1:
            #             split1_site += sub_aln[seq][j];
            #         if seq in split2:
            #             split2_site += sub_aln[seq][j];
            #     # Get the allele for every sequence in the sub-alignment and add it to the appropriate sites                    

            #     if any(state in site for state in aln_skip_chars):
            #         continue;
            #     # We only care about sites with full information for the quartet, so skip others here

            #     uniq_site = set(site);
            #     # The alleles at the current site

            #     if len(uniq_site) > 1:
            #         quartet_scores[node][quartet]['variable-sites'] += 1;
            #     # If there is more than one allele, the site is variable

            #         if all(site.count(allele) == 2 for allele in uniq_site):
            #         #if len(uniq_site) == 2:
            #             quartet_scores[node][quartet]['decisive-sites'] += 1;
            #         # If all alleles have a count of 2 at this site, it is decisive

            #             #print(j, site);

            #             if split1_site[0] == split1_site[1] and split2_site[0] == split2_site[1]:# and split1_site[0] != split2_site[0]:
            #             #if split1_site.count(split1_site[0]) == len(split1_site) and split2_site.count(split2_site[0]) == len(split2_site):
            #                 quartet_scores[node][quartet]['concordant-sites'] += 1;
            #             # If the alleles from split1 match and the alleles from split2 match, the site is concordant
            #             elif split1_site[0] == split2_site[0] and split1_site[1] == split2_site[1]:
            #                 quartet_scores[node][quartet]['disco1-sites'] += 1;
            #             elif split1_site[0] == split2_site[1] and split1_site[1] == split2_site[0]:
            #                 quartet_scores[node][quartet]['disco2-sites'] += 1;

            ## End site loop

            # if quartet_scores[node][quartet]['decisive-sites'] != 0:
            #     quartet_scores[node][quartet]['scf'] = quartet_scores[node][quartet]['concordant-sites'] / quartet_scores[node][quartet]['decisive-sites'];
            # Calculate the scf for the current quartet
        ## End quartet loop

        # for quartet in quartet_scores[node]:
        #     if quartet_scores[node][quartet]['scf'] != "NA":
        #         locus_scf[node]['scf-sum'] += quartet_scores[node][quartet]['scf'];
        #         locus_scf[node]['num-quartets'] += 1;
        # # For all the quartets with sCF calculated, sum the sCF for this locus to be averaged later

        # if locus_scf[node]['num-quartets'] != 0:
        #     locus_scf[node]['avg-scf'] = locus_scf[node]['scf-sum'] / locus_scf[node]['num-quartets'];
        #     node_scf.append(locus_scf[node]['scf-sum'] / locus_scf[node]['num-quartets']);
        # Average all sCF from each quartet for this locus
    ## End node loop

    return locus, quartet_scores;

#############################################################################

def scf(globs, st, alns, proc_pool):
# A function to calculate site concordance factors for each input locus
# sCFs are calculate in two ways:
# 1. Average sCF across all nodes per LOCUS (written to aln stats file)
# 2. Average sCF across all loci per NODE (written to scf stats file)
# In both cases, a number of quartets (100) are sampled for each branch in each tree and sites are counted
# and averages are taken across quartets. See: https://doi.org/10.1093/molbev/msaa106

    step = "Sampling quartets";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    sampled_quartets = st.sampleQuartets();
    step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
    # Sample quartets for all nodes

    num_alns = len(alns);

    step = "Calculating per-locus sCF";
    step_start_time = CORE.report_step(globs, step, False, "Processed 0 / " + str(num_alns) + " alns...", full_update=True);
    # Status update

    scf_dict = {};
    # The dictionart of averages over quartets

    quartet_counts = {};
    # The dictionary of site counts per quartet

    for node in st.nodes:
        if node in st.tips or node == st.root:
            continue;
        # Cannot calculate sCF for tips, the root, or node descendant from the root

        scf_dict[node] = { 'variable-sites' : 0, 'decisive-sites' : 0, 
                                'concordant-sites' : 0, 'disco1-sites' : 0, 'disco2-sites' : 0,
                                 'quartet-scf-sum' : 0, 'total-quartets' : 0, 'avg-quartet-scf' : "NA",
                                 'scf' : "NA", 'sdf1' : "NA", 'sdf2' : "NA" };
        # The main scf dictionary with AVERAGES across quartets per node

        quartet_counts[node] = {};
        for quartet in sampled_quartets[node]:
            quartet_counts[node][quartet] = { 'variable-sites' : 0, 'decisive-sites' : 0, 
                                'concordant-sites' : 0, 'disco1-sites' : 0, 'disco2-sites' : 0 };
        # The dictionary of site COUNTS per quartet per node
    ## Initialize the dictionaries to calculate average sCF per node across all loci

    # if globs['qstats']:
    #     qstats_file = os.path.join(globs['outdir'], "quartet-stats.csv");
    #     qfile = open(qstats_file, "w");
    #     headers = ["locus","node","quartet","variable-sites","decisive-sites","concordant-sites"];
    #     qfile.write(",".join(headers) + "\n");
    # For --qstats, creates and opens a file to write site counts for each quartet within the pool loop

    with proc_pool as pool:
        counter = 0;
        # A counter to keep track of how many loci have been completed
        #for locus in globs['alns']:
        #    result = locusSCF((locus, globs['alns'][locus], sampled_quartets, st, globs['aln-skip-chars']));

        for result in pool.imap_unordered(locusSCF, ((locus, alns[locus], sampled_quartets, st, globs['aln-skip-chars']) for locus in alns)):
            # Loop over every locus in parallel to calculate sCF per node

            locus, quartet_scores = result;
            # Unpack the current result

            # globs['aln-stats'][locus]['node-scf-sum'] = sum([ n for n in node_scf ]);
            # globs['aln-stats'][locus]['num-nodes'] = len(node_scf);
            # globs['aln-stats'][locus]['node-scf-avg'] = "NA";
            # globs['aln-stats'][locus]['perc-low-scf-nodes'] = "NA";            
            # globs['aln-stats'][locus]['low-scf-nodes'] = len([ n for n in node_scf if n < globs['min-scf'] ]);
            # if globs['aln-stats'][locus]['num-nodes'] != 0:
            #     globs['aln-stats'][locus]['node-scf-avg'] = globs['aln-stats'][locus]['node-scf-sum'] / globs['aln-stats'][locus]['num-nodes'];
            #     globs['aln-stats'][locus]['perc-low-scf-nodes'] = globs['aln-stats'][locus]['low-scf-nodes'] / globs['aln-stats'][locus]['num-nodes'];
            # # Average sCF per locus stored with the other alignment stats

            for node in quartet_scores:
                for q in sampled_quartets[node]:

                    # if globs['qstats']:
                    #     q_str = ";".join(q[0]) + ";" + ";".join(q[1]);
                    #     outline = [locus, node, q_str] + [str(quartet_scores[node][q][col]) for col in headers[3:]];
                    #     qfile.write(",".join(outline) + "\n");
                    # For --qstats, writes the quartet stats to a file for the current quartet

                    quartet_counts[node][q]['variable-sites'] += quartet_scores[node][q]['variable-sites'];
                    quartet_counts[node][q]['decisive-sites'] += quartet_scores[node][q]['decisive-sites'];
                    quartet_counts[node][q]['concordant-sites'] += quartet_scores[node][q]['concordant-sites'];
                    quartet_counts[node][q]['disco1-sites'] += quartet_scores[node][q]['disco1-sites'];
                    quartet_counts[node][q]['disco2-sites'] += quartet_scores[node][q]['disco2-sites'];
                    # Sum sites for this quartet for this node

                    # if quartet_scores[node][q]['decisive-sites'] != 0:
                    #     globs['scf'][node]['quartet-scf-sum'] += quartet_scores[node][q]['concordant-sites'] / quartet_scores[node][q]['decisive-sites'];
                    #     globs['scf'][node]['total-quartets'] += 1;
                    # If there are decisive sites, calculate sCF for this quartet
            # For every quartet in every node, calculate sCF and sum up the number of sites in order to average
            # across nodes later

            counter += 1;
            if counter % 100 == 0:
                cur_scf_time = CORE.report_step(globs, step, step_start_time, "Processed " + str(counter) + " / " + str(num_alns) + " alns...", full_update=True);
            # A counter and a status update every 100 loci
        ## End imap locus loop
    ## End pool

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success", full_update=True);
    # Status update
    
    step = "Averaging site counts";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    for node in scf_dict:
    # Average site counts for every node

        node_concordant = [];
        node_scfs = [];
        node_d1 = [];
        node_sdf1 = [];
        node_d2 = [];
        node_sdf2 = [];
        node_decisive = [];
        # Lists of site counts and cfs to build for each quartet

        for q in sampled_quartets[node]:
        # Append counts from every quartet on this node
            node_scfs.append(quartet_counts[node][q]['concordant-sites'] / quartet_counts[node][q]['decisive-sites']);
            node_sdf1.append(quartet_counts[node][q]['disco1-sites'] / quartet_counts[node][q]['decisive-sites']);
            node_sdf2.append(quartet_counts[node][q]['disco2-sites'] / quartet_counts[node][q]['decisive-sites']);
            node_concordant.append(quartet_counts[node][q]['concordant-sites']);
            node_d1.append(quartet_counts[node][q]['disco1-sites']);
            node_d2.append(quartet_counts[node][q]['disco2-sites']);
            node_decisive.append(quartet_counts[node][q]['decisive-sites']);

            if quartet_counts[node][q]['decisive-sites'] != 0:
                scf_dict[node]['total-quartets'] += 1;
        ## End quartet loop

        scf_dict[node]['scf'] = CORE.mean(node_scfs);
        scf_dict[node]['concordant-sites'] = CORE.mean(node_concordant);
        scf_dict[node]['sdf1'] = CORE.mean(node_sdf1);
        scf_dict[node]['disco1-sites'] = CORE.mean(node_d1);
        scf_dict[node]['sdf2'] = CORE.mean(node_sdf2);
        scf_dict[node]['disco2-sites'] = CORE.mean(node_d2);

        scf_dict[node]['decisive-sites'] = CORE.mean(node_decisive);
        # Calculate averages across quartets of sites and concordance factors

        # globs['scf'][node]['scf'] = globs['scf'][node]['concordant-sites'] / globs['scf'][node]['decisive-sites'];
        #globs['scf'][node]['avg-quartet-scf'] = globs['scf'][node]['quartet-scf-sum'] / globs['scf'][node]['total-quartets'];
    ## End site count average loop

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
    # For every node, average the sCFs per quartet across all loci

    # if globs['qstats']:
    #     qfile.close()
    # For --qstats, closes the quartet stats file.

    return scf_dict;
    
#############################################################################