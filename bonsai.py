#!/usr/bin/env python3
#############################################################################
# bonsai is a script to remove branches from a tree to maximize concordance/support
# of the underlying loci. This is the main interface.
#
# Gregg Thomas
# Fall 2021
#############################################################################

# pare.py -t full_coding_iqtree_astral.cf.rooted.tree -o test --overwrite -i 30

# PRUNING: Removing a branch from a tree
# PARING: Pruning one of the two descending clades from a branch to increase its length/support

#############################################################################

import sys
import os
import re
import lib.core as CORE
import lib.params as params
import lib.opt_parse as OP
import lib.seq as SEQ
import lib.tree as TREE
import lib.treeio as TREEIO
import lib.cf as CF
import lib.pare as PARE

#############################################################################

if __name__ == '__main__':
# Main is necessary for multiprocessing to work on Windows.

    globs = params.init();
    ## Get the global params as a dictionary.
    
    print("\n" + " ".join(sys.argv) + "\n");

    if any(v in sys.argv for v in ["--version", "-version", "--v", "-v"]):
        print("# bonsai version " + globs['version'] + " released on " + globs['releasedate'])
        sys.exit(0);
    # The version option to simply print the version and exit.
    # Need to get actual PhyloAcc version for this, and not just the interface version.

    print("#");
    print("# " + "=" * 125);
    #print(CORE.welcome());
    if "-h" not in sys.argv:
        print("            Heuristic tree paring\n");
    ## A welcome banner.

    globs = OP.optParse(globs);
    ## Getting the input parameters from optParse.

    if globs['info']:
        print("# --info SET. EXITING AFTER PRINTING PROGRAM INFO...\n#")
        sys.exit(0);
    if globs['norun']:
        print("# --norun SET. EXITING AFTER PRINTING OPTIONS INFO...\n#")
        sys.exit(0);
    ## Early exit options

    step_start_time = CORE.report_step(globs, "", "", "", start=True);
    ## Initialize the step headers

    # Initializaiton
    ####################

    if not globs['debug-tree']:
        globs = TREEIO.readST(globs);
    else:
        TREE.debugTree(globs);
        sys.exit(0);

    if globs['label-tree']:
        globs['st'].subtree = globs['st'].genSubtrees();
        globs['st'].tree_str = globs['st'].subtree[globs['st'].root];
        globs['st'].showAttrib("type", "label", "desc");
        print("\n" + globs['st'].tree_str +"\n");
        print();
        sys.exit(0);
    #else:
    #    CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: Original tree with node labels:\t" + globs['st'].tree_str);
    # Print the tree and exit if --labeltree is set

    ## Read species tree
    ####################

    if globs['gt-input']:
        globs, num_rooted, num_unrooted = TREEIO.readGT(globs);

    ## Read gene trees
    ####################

    if globs['gt-input'] and not globs['use-labels']:
        step = "Calculating gCF";
        step_start_time = CORE.report_step(globs, step, False, "In progress...");  
        globs['gcf-dict'], total_trees = CF.gcf(globs['st'], globs['gts']);
        step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
        # gCF

        topo_counts, num_trees, num_mismatch = CF.countTopos(globs, globs['st'], globs['gts']);
        TREEIO.writeTopoCounts(globs, topo_counts, num_trees, globs['topo-count-file']);
        # Topology counts and RF

    ## gCF
    ####################
  
    if globs['aln-dir']:
        globs['alns'], globs['num-alns'] = SEQ.readSeq(globs['aln-dir'], globs);
    
    ## Read alignments    
    ####################

    if globs['aln-stat-file']:
            globs['aln-stats'] = SEQ.alnStats(globs, globs['alns'], globs['aln-pool']);
            # Calculate some basic alignment stats

            globs = SEQ.writeAlnStats(globs, globs['aln-stats'], globs['aln-stat-file']);
            # Write out the alignment summary stats

    ## Calculate alignment stats
    ####################

    if globs['scf-quartets']:
        globs['scf-dict'] = CF.scf(globs, globs['st'], globs['alns'], globs['scf-pool']);
        # Sample quartets and calculate sCF

    ## Calculate sCF
    ####################

    if (globs['gt-input'] and not globs['use-labels']) or globs['scf-quartets']:
        globs['st'] = TREEIO.writeCF(globs['st'], globs['gcf-dict'], globs['scf-dict'], globs['st-cf-stat-file'], globs['st-cf-tree-file'], globs)

    ## CF output
    ####################

    if globs['cf-only']:
        CORE.endProg(globs, "");

    ## Exit here if --cf is set
    ####################

    if globs['prune-file']:
    # If a file with a list of branches to prune form the species tree is given, do not do the iterative
    # paring, rather just pruned those branches here

        step = "Reading branches to prune";
        step_start_time = CORE.report_step(globs, step, False, "In progress...", full_update=True);
        # Status update

        globs, prune_branches, prune_clades = TREEIO.readBranches(globs, globs['prune-file']);
        # Read the file with branches to prune

        globs['prune-branches'] = prune_branches;
        globs['prune-clades'] = prune_clades;
        globs['pruned-tips'] = list(set().union(*globs['prune-clades']));
        # Save the sets of branches, clades, and tips to prune

        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(len(globs['prune-branches'])) + " branches will be pruned", full_update=True);
        # Status update

        ####################

        step = "Pruning species tree";
        step_start_time = CORE.report_step(globs, step, False, "In progress...");
        # Stats update
        pruned_tree_str = globs['st'].Prune(globs['prune-branches']);
        globs['st'] = TREE.Tree(pruned_tree_str);
        step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
        # Status update

        ## Prune the species tree
        ####################

        if globs['gt-input']:
            globs, pruned_gts = PARE.pruneGT(globs, globs['gts'], globs['pruned-tips']);
            globs['gts'] = pruned_gts

            step = "Recalculating gCF";
            step_start_time = CORE.report_step(globs, step, False, "In progress...");  
            pruned_gcf_dict, total_trees = CF.gcf(globs['st'], globs['gts']);
            step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
            # gCF

        ## Prune the gene trees of the same tips if they are given and recalculate gCF
        ####################

        final_msg = "# Pruned " + str(len(globs['prune-branches'])) + " total branches and removed " + str(len(globs['pruned-tips'])) + " total tips.";
        # Info message

    ## End -prune block
    ####################

    else:
    # If a file isn't given with -p, run the iterative pruning algorithm

        if globs['exempt-file']:
        # If a file is given with branches to exempt from pruning, read it here

            step = "Reading exempt branches";
            step_start_time = CORE.report_step(globs, step, False, "In progress...", full_update=True);
            # Status update

            globs, exempt_branches, exempt_clades = TREEIO.readBranches(globs, globs['exempt-file']);
            # Read the file

            globs['exempt-branches'] = exempt_branches;
            globs['exempt-clades'] = exempt_clades;
            # Save sets of branches and clades to exempt

            step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(len(globs['exempt-branches'])) + " branches will be exempt from paring.", full_update=True);
            # Status update

        ## Reading the file with branches exempt from pruning
        ####################

        pare = True;
        # Bool for the paring loop -- it will run until one of the stopping conditions is met:
        # 1. The number of species pruned is over the max allowed with -m
        # 2. No more branches are pruned
        # 3. The number of iterations reaches the max allowed with -i

        iteration = 0;
        # Iteration counter

        while pare:
        # The main paring loop

            iteration += 1;
            # Iterate the iteration counter

            ####################

            CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 50);
            step = "Paring iteration " + str(iteration);
            step_start_time = CORE.report_step(globs, step, False, "In progress...", full_update=True);
            # Status update

            globs, bl_threshold, pruned_tree_str, pared_branches, pruned_tips, over_max_spec = PARE.pare(globs, globs['st'], iteration);
            # Call the paring algorithm for the current tree

            pruned_tree = TREE.Tree(pruned_tree_str);
            #pruned_tree.showAttrib("length", "label");
            # Read the returned pruned tree as a Tree

            ####################

            if over_max_spec:
                step_start_time = CORE.report_step(globs, step, step_start_time, "Truncated", full_update=True);
                # Status update      

                CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: Branch length threshold for iteration:\t" + str(bl_threshold));
                CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: This iteration would remove " + str(len(pruned_tips)) + " tips, which puts the total number pruned over the maximum limit (" + str(globs['max-spec']) +  "). Paring complete.");
                # Info update

                pare = False;
            ## Stopping condition 1: the number of species pruned is over the max allowed with -m

            ####################

            elif len(pruned_tips) == 0:
                step_start_time = CORE.report_step(globs, step, step_start_time, "Success: No branches pared this iteration. Exiting.", full_update=True);
                # Status update              

                CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: Branch length threshold for iteration:\t" + str(bl_threshold));
                # Info update

                pare = False;
            ## Stopping condition 2: No species were pruned in the last iteration

            ####################

            else:
            # If the first 2 stopping conditions aren't met, process the current iteration

                step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(len(pruned_tips)) + " tips removed, " + str(len(pared_branches)) + " branches pared.", full_update=True);
                # Status update

                CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: Branch length threshold for iteration:\t" + str(bl_threshold));
                # Info update

                step = "Writing iteration " + str(iteration) + " files";
                step_start_time = CORE.report_step(globs, step, False, "In progress...");
                # Stats update

                with(open(os.path.join(globs['outdir'], "iter-" + str(iteration) + "-pruned-spec.txt"), "w")) as pruned_file:
                    for tip in pruned_tips:
                        pruned_file.write(tip + "\n");
                # Write the species pruned in this iteration to a file

                with(open(os.path.join(globs['outdir'], "iter-" + str(iteration) + "-pared-branches.txt"), "w")) as pared_file:
                    for n in pared_branches:
                        pared_file.write(n + "\n");
                # Write the branches pruned in this iteration to a file

                step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
                # Status update

                ## Writing the branches and tips pruned in this iteration to files
                ####################

                if globs['gt-input']:
                # If gene trees are provided, prune them here and recalculate gCF

                    globs, pruned_gt = PARE.pruneGT(globs, globs['gts'], pruned_tips);
                    globs['gts'] = pruned_gt;     

                    ## Prune gene trees
                    ####################

                    step = "Recalculating gCF";
                    step_start_time = CORE.report_step(globs, step, False, "In progress...");
                    # Status update
                    
                    pruned_gcf_dict, total_trees = CF.gcf(pruned_tree, globs['gts']);
                    # Count quartet topos

                    gcf_labels = {};
                    for node in pruned_gcf_dict:
                        gcf = round(pruned_gcf_dict[node]["concordant"] / pruned_gcf_dict[node]["decisive"], 3);
                        gcf_labels[node] = str(gcf);
                    # Calculate gCF

                    pruned_tree_str = pruned_tree.addLabel(gcf_labels);
                    # Add gCF labels to the branches of the pruned tree

                    globs['st'] = TREE.Tree(pruned_tree_str);
                    # Read the tree with gCF values as a Tree and set it as the species tree

                    step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
                    # Status update

                ## End gene tree block
                ####################

                else:
                    globs['st'] = pruned_tree;
                # If no gene trees are given, just set the pruned tree as the species tree

                ####################

                with(open(os.path.join(globs['outdir'], "iter-" + str(iteration) + "-pared.tre"), "w")) as pared_tree_file:
                    pared_tree_file.write(pruned_tree_str);
                # Write the pared tree from this iteration to a file  

                ####################

                globs['bl-thresholds'].append(bl_threshold);
                globs['pared-branches'].append(pared_branches);
                globs['total-pared-branches'] += len(pared_branches);
                globs['pruned-tips'].append(pruned_tips);
                globs['total-pruned-tips'] += len(pruned_tips);
                globs['iter-trees'].append(pruned_tree_str);
                # Add the pared branches and pruned tips from this iteration to the global lists for all iterations

                ####################

                if iteration == globs['max-iterations']:
                    iteration += 1;
                    pare = False;
                ## Stopping condition 3: The number of iterations is over the max allowed by -i

                ####################

        # The paring loop
        ####################

        CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 50);
        step = "Writing summary stats to log";
        step_start_time = CORE.report_step(globs, step, False, "In progress...");
        # Status update    

        headers = ["iteration", "branch length threshold", "branches pared", "tips pruned", "tree"];
        CORE.printWrite(globs['logfilename'], 3, "\t".join(headers));
        # Headers for the iteration output in the logfile

        for i in range(iteration-1):
        # Write output for every iteration

            iter_tree = globs['iter-trees'][i] + ";";
            # The tree from the current iteration

            if i == 0:
                outline = [ "0", "NA", "NA", "NA", iter_tree ];
            # Write the original tree as the 0th iteration
            else:
                iter_str = str(i);
                outline = [ iter_str, str(globs['bl-thresholds'][i-1]), str(len(globs['pared-branches'][i-1])), str(len(globs['pruned-tips'][i-1])), iter_tree ];
            CORE.printWrite(globs['logfilename'], 3, "\t".join(outline));

        tips_file = os.path.join(globs['outdir'], "all-pruned-tips.txt");
        with open(tips_file, "w") as tipsfile:
            for iter_tips in globs['pruned-tips']:
                for tip in iter_tips:
                    tipsfile.write(tip + "\n");
        # Write the pruned tips to a file

        step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
        # Status update

        ####################

        final_msg = "# Pared " + str(globs['total-pared-branches']) + " total branches and removed " + str(globs['total-pruned-tips']) + " total tips.";
        # Info message

    ## End paring block
    ####################
    
    if globs['gt-input']:
    # If gene trees are given as input, this block handles their related output
        globs = TREEIO.writeGT(globs, globs['gts'], globs['prune-gt-file']);
        # Write the pruned gene trees to a file

        topo_counts, num_trees, num_mismatch = CF.countTopos(globs, globs['st'], globs['gts']);
        TREEIO.writeTopoCounts(globs, topo_counts, num_trees, globs['topo-count-pruned-file']);
        # Topology counts and RF
    
    ## End gene tree output block
    ####################    

    pruned_scf_dict = {};
    if globs['aln-dir']:
    # If alignments are given as input, this block handles their related output

        globs, pruned_alns = SEQ.subsetAlns(globs, globs['pruned-tips']);
        globs = SEQ.writeAlns(globs, pruned_alns, globs['prune-aln-dir']);
        # Subset and write the alignments based on the pruned species

        ####################

        if globs['aln-stat-file']:
            globs['aln-stats-pruned'] = SEQ.alnStats(globs, pruned_alns, globs['aln-pool-prune']);
            globs = SEQ.writeAlnStats(globs, globs['aln-stats-pruned'], globs['aln-stat-prune-file']);
            # Calculate and write some basic alignment stats for the pruned alignments

        if globs['scf-quartets']:
            pruned_scf_dict = CF.scf(globs, globs['st'], pruned_alns, globs['scf-pool-prune']);
        # Sample quartets and calculate sCF if specified in input params

    ## End alignment output block
    ####################

    if (globs['gt-input'] and not globs['use-labels']) or globs['scf-quartets']:
        globs['st'] = TREEIO.writeCF(globs['st'], pruned_gcf_dict, pruned_scf_dict, globs['st-cf-stat-prune-file'], globs['st-final-file'], globs)
    # If CF has been calculated, add gCF labels to the branches of the pruned tree
    else:
        with open(globs['st-final-file'], "w") as treefile:
            final_tree_str = globs['st'].addLabel(globs['st'].label);
    # Otherwise, add the original labels back to the tree and write it to a file            

    if globs['prune-file']:
        print("\n# Pruned species tree: ")
        print("\n" + pruned_tree_str + "\n");

    CORE.endProg(globs, final_msg);
    # A nice way to end the program

#############################################################################

