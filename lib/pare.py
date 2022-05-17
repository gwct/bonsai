#############################################################################
# The main tree paring functions
# Gregg Thomas
#############################################################################

import sys
import re
import random
import numpy as np
import lib.core as CORE
import lib.tree as TREE

#############################################################################

def pare(globs, cur_tree, iteration, debug=False):

    if debug:
        cur_tree.showAttrib("length", "label");

    cur_tree.subtree = cur_tree.genSubtrees();
    cur_tree.tree_str = cur_tree.subtree[cur_tree.root];

    cur_tree.bl = { n : float(cur_tree.bl[n]) if cur_tree.bl[n] != "NA" else 0.0 for n in cur_tree.bl };
    cur_tree.label = { n : float(cur_tree.label[n]) if cur_tree.label[n] != "NA" else 0.0 for n in cur_tree.label };

    bls = [ cur_tree.bl[node] for node in cur_tree.internals if node != cur_tree.root ];

    ####################

    cur_bl_threshold = np.percentile(bls, globs['bl-percentile']);
    if debug:
        print("BL THRESHOLD: ", cur_bl_threshold);
        print("--------");
    # Calculate the branch length threshold from the distribution of bls
    # Branches with length lower than this will be considered for paring
    ## TODO: Should this be updated every time or not?

    ####################

    branches_to_pare = [];
    spec_to_prune = [];
    # Keep track of both the branches that will be pared and the species
    # descending from those branches that will be pruned

    criteria = "max-gcf";
    if criteria not in ["deterministic", "random", "min-spec", "max-gcf"]:
        sys.exit("\n\nINTERNAL: Invalid clade selection criteria specified: " + criteria + "\n\n");
    # One of: deterministic, random, min-spec, max-gcf

    for n in cur_tree.internals:
        if n == cur_tree.root:
            continue;
        # Skip tips and the root when getting branches to pare
        ## TODO: Consider tips for non-Astral trees?

        if cur_tree.bl[n] < cur_bl_threshold and cur_tree.label[n] < globs['gcf-threshold']:
        # If the branch length of the current node is lower than both the bl threshold and the gcf threshold, select
        # species to pare

            clade1 = cur_tree.getClade(cur_tree.desc[n][0]);
            clade2 = cur_tree.getClade(cur_tree.desc[n][1]);
            # Get all species in the clades descending from the current descendants

            if criteria == "deterministic":
                cur_prune_clade = clade1;
                cur_prune_node = cur_tree.desc[n][0];
                prune_index = 0;
                # Variables to keep track of the pruned branch
            # Deterministic pruning (good for debugging)

            ##########

            elif criteria == "random":
                cur_prune_clade = random.choice([clade1, clade2]);
                
                if cur_prune_clade == clade1:
                    cur_prune_node = cur_tree.desc[n][0];
                    prune_index = 0;
                else:
                    cur_prune_node = cur_tree.desc[n][1];
                    prune_index = 1;
                # Variables to keep track of the pruned branch
            # Randomly choose one of the clades to pare

            ##########

            elif criteria == "min-spec":
                if len(clade1) <= len(clade2):
                    cur_prune_clade = clade1;
                    cur_prune_node = cur_tree.desc[n][0];
                    prune_index = 0;
                else:
                    cur_prune_clade = clade2;
                    cur_prune_node = cur_tree.desc[n][1];
                    prune_index = 1;
                # Variables to keep track of the pruned branch
            # Minimize species removed

            ##########

            elif criteria == "max-gcf":
                gcfs = [];
                # We will generate two lists of gcfs to average: one from each descendant node

                for d in [0,1]:
                # Loop through both descendant nodes

                    cur_desc = cur_tree.desc[n][d];

                    if debug:
                        print("CUR DESC: ", cur_desc);
                    # Debug statements

                    if cur_tree.type[cur_desc] != 'tip':
                    # If the tree isn't a tip, get gcfs

                        cur_gcfs = [];
                        # The list of subtrees for the current descendant node

                        cur_clade = cur_tree.getClade(cur_desc, full=True);
                        for subtree_node in cur_clade:
                            if cur_tree.type[subtree_node] != "tip":
                                cur_gcfs.append(cur_tree.label[subtree_node]);

                        # cur_subtree_str = cur_tree.subtree[cur_desc];
                        # cur_subtree_str = re.sub("<[\d]+>_", "", cur_subtree_str);
                        # if debug:
                        #     print("CUR SUBTREE: ", cur_subtree_str);
                        # # Get the subtree descending from the current descendant node and remove
                        # # the previous labels

                        # cur_subtree = TREE.Tree(cur_subtree_str);
                        # # Parse the subtree

                        # for sub_n in cur_subtree.internals:
                        #     if cur_subtree.label[sub_n] not in ["NA", ""]:
                        #         cur_gcfs.append(float(cur_subtree.label[sub_n]));
                        # # For every node in the subtree, if there is a gcf add it to the current list of gcfs

                        cur_gcfs.append(cur_tree.label[cur_desc]);
                        # Also add the gcf of the current descendant itself to the list of gcfs

                        gcfs.append(cur_gcfs);

                        if debug:
                            print("CUR GCFS: ", cur_gcfs);
                        # Add the current list of gcfs to the list of lists of gcfs
                    else:
                    # If the descendant node is a tip, there is no gcf so make sure it isn't chosen by making
                    # it 0.0
                        gcfs.append([0.0]);
                    # gcf gathering statements
                ## Descendant loop

                avg_gcfs = [];
                for gcf_list in gcfs:
                    avg_gcfs.append(CORE.mean(gcf_list));
                if debug:
                    print("AVG GCFS: ", avg_gcfs);
                # Average the gcfs from both subtrees

                if avg_gcfs[0] <= avg_gcfs[1]:
                    if debug:
                        print("HI");
                    cur_prune_clade = clade1;
                    cur_prune_node = cur_tree.desc[n][0];
                    prune_index = 0;
                else:
                    cur_prune_clade = clade2;
                    cur_prune_node = cur_tree.desc[n][1];
                    prune_index = 1;
                # Variables to keep track of the pruned branch
                if debug:
                    print("--------");
                # Pick the subclade with the lowest average gcf to pare
            # Maximize gCF of clade retained

            ##########

            if any(clade.issubset(set(cur_prune_clade)) for clade in globs['exempt-clades']):
            # If the selected clade to prune matches exactly one of the exempt clades, switch to its sister here

                CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: Branch " + n + " chosen for paring. Branch " + cur_prune_node + " chosen for pruning, but is exempt. Switching to sister branch.");
                # Info update

                if prune_index == 0:
                    cur_prune_clade = clade2;
                    cur_prune_node = cur_desc[1];
                    prune_index = 1;
                else:
                    cur_prune_clade = clade1;
                    cur_prune_node = cur_desc[0];
                    prune_index = 0;
                # Switch to the sister clade depending on the prune index
            ## End first check for exempt clade                  

            if any(clade.issubset(set(cur_prune_clade)) for clade in globs['exempt-clades']):
            # If the sister clade also exactly matches an exempt clade, do not prune anything

                CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: Branch " + n + " chosen for paring. Sister branch " + cur_prune_node + " is also exempt. Skipping paring.");
                # Info update

                cur_prune_clade = [];
                cur_prune_node = False;
                prune_index = False;
                # Switch the clade pruning variables to empty/false here
            ## End second check for exempt clade

            if cur_prune_clade:
            # If there is a clade to prune
            
                num_spec_to_pare = len(cur_prune_clade);
                # if debug:
                #     print(num_spec_to_pare);
                # The number of species to pare

                if num_spec_to_pare <= globs['branch-max-spec']:
                # If the number of species to prune doesn't exceed the max per branch, add them
                    CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: Branch " + n + " chosen for paring. Branch " + cur_prune_node + " will be pruned and remove " + str(num_spec_to_pare) + " species");
                    # Info update

                    branches_to_pare.append(n);
                    # Add this branch to the list of branches to pare

                    spec_to_prune += cur_prune_clade;
                    # Add all species in this clade to the to_pare list
                
                else:
                    CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: Branch " + n + " chosen for paring. Branch " + cur_prune_node + " was to be pruned but will be retained because it exceeds the max number of species to remove in a given branch (" + str(num_spec_to_pare) + ">" + str(globs['branch-max-spec']) + ")");
                    # Info update
                # Otherwise print an update that we are skipping this pruning

                #print(cur_prune_node);
                #print(cur_prune_clade);
                #print("-----")

            ## End clade prune block.
    ## End loop to select clades to pare

    spec_to_prune = list(set(spec_to_prune));
    # Get a uniqe list of species to pare

    if len(spec_to_prune) + globs['total-pruned-tips'] > globs['total-max-spec']:
        return globs, cur_bl_threshold, "NA", branches_to_pare, spec_to_prune, True;
    # If pruning tips this iteration would put us over the max number of species to prune, truncate
    # the iteration here and exit to main

    if debug:
        print("BRANCHES TO PRUNE: ", branches_to_pare);
        print("SPECIES TO PRUNE: ", spec_to_prune);
        print("# SPECIES TO PRUNE: ", len(spec_to_prune));
        print("--------");

    # This block handles selection of species to prune/clades to pare
    ####################

    for spec in spec_to_prune:
        for i in range(len(globs['exempt-clades'])):
            if spec in globs['exempt-clades'][i]:
                globs['exempt-clades'][i].remove(spec);
    #print(globs['exempt-clades']);
    # For every species/tip we are pruning, we need to remove that tip from the exempt clades as well, otherwise
    # we risk removing the whole clade in a later iteration

    ####################

    pruned_tree_str = cur_tree.Prune(spec_to_prune);

    return globs, cur_bl_threshold, pruned_tree_str, branches_to_pare, spec_to_prune, False;
    

#############################################################################

def pruneGT(globs, gts, tips_to_prune):

    step = "Pruning gene trees";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Stats update

    pruned_gts = {};
    gt_pruned, gt_skipped = 0, 0;
    for gt_id in gts:
        #print(gt_id);
        # print("orig st: ", cur_st.labeled_topo_str);
        cur_gt = gts[gt_id];
        # Get trees

        if cur_gt.num_tips < 4 or len(cur_gt.tips) - len(tips_to_prune) < 4:
            gt_skipped += 1;
            globs['warnings'] += 1;
            CORE.printWrite(globs['logfilename'], -2, "# WARNING: Gene tree on line " + str(gt_id) + " has too few tips after pruning to match the species tree. Removing from this iteration.");
            pruned_gts[gt_id] = "NA";
            continue;
        
        # Throw a warning if there aren't enough tips left after pruning the gene tree and skip

        cur_gt = cur_gt.Prune(tips_to_prune);
        cur_gt = TREE.Tree(cur_gt);
        #print("prune gt:", cur_gt.tree);
        #CORE.printWrite(globs['logfilename'], -2, "# INFO: The gene tree on line " + str(gt_id) + " was pruned of the following " + str(len(gt_tips_to_prune)) + " species because they were not in the species tree: " + str(gt_tips_to_prune));
        pruned_gts[gt_id] = cur_gt;
        gt_pruned += 1;
        # Prune gene tree to match tips in species tree if necessary

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
    # Status update

    return globs, pruned_gts;

#############################################################################