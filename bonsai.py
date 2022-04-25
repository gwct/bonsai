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
import lib.tree as TREE
import lib.treeio as TREEIO
import lib.pare as PARE

#############################################################################

if __name__ == '__main__':
# Main is necessary for multiprocessing to work on Windows.

    globs = params.init();
    # Get the global params as a dictionary.
    
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
    # A welcome banner.

    globs = OP.optParse(globs);
    # Getting the input parameters from optParse.

    if globs['info']:
        print("# --info SET. EXITING AFTER PRINTING PROGRAM INFO...\n#")
        sys.exit(0);
    if globs['norun']:
        print("# --norun SET. EXITING AFTER PRINTING OPTIONS INFO...\n#")
        sys.exit(0);
    # Early exit options

    ####################

    step_start_time = CORE.report_step(globs, "", "", "", start=True);
    # Initialize the step headers

    ####################

    step = "Parsing input species tree";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    globs = TREEIO.readST(globs);

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: tree with " + str(globs['num-tips']) + " tips and " + str(globs['num-internals']) + " internal nodes read.");
    # Status update
    
    if globs['label-tree']:
        print("\n" + globs['labeled-st-str'] +"\n");
        sys.exit(0);
    else:
        CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: Original tree with node labels:\t" + globs['labeled-st-str']);
    # Print the tree and exit if --labeltree is set

    ####################

    step = "Parsing input gene trees";
    step_start_time = CORE.report_step(globs, step, False, "In progress...", full_update=True);
    # Status update    

    globs = TREEIO.readGT(globs);

    if globs['gt-init']:
        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(len(globs['gt-init'])) + " gene trees read", full_update=True);
        if globs['gt-input-empty']:
            CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: Skipped " + str(globs['gt-input-empty']) + " lines because they were empty");
        if globs['gt-input-skipped']:
            CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: Skipped " + str(globs['gt-input-skipped']) + " lines because they couldn't be read as trees");
    else:
        step_start_time = CORE.report_step(globs, step, step_start_time, "Failure: no gene trees read", full_update=True);
        print("\n\n");
        CORE.errorOut("2", "Error reading gene trees! Couldn't parse any trees from file. Make sure they are formatted as rooted, Newick trees.", globs);
    # Status update

    ####################

    if globs['exempt-file']:
        step = "Reading exempt branches";
        step_start_time = CORE.report_step(globs, step, False, "In progress...", full_update=True);
        # Status update

        globs = TREEIO.readExempt(globs);

        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(len(globs['exempt-branches'])) + " branches will be exempt from paring.", full_update=True);
        # Status update
    ## Reading the file with branches exempt from pruning

    ####################

    pare = True;
    # Boolean for the paring loop

    iteration = 0;
    # Iteration counter

    cur_tree = globs['orig-st-str'];
    # For the first iteration, the tree will be the input tree.

    while pare:
        iteration += 1;
        # Iterate the iteration counter

        ####################

        CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 50);
        step = "Paring iteration " + str(iteration);
        step_start_time = CORE.report_step(globs, step, False, "In progress...", full_update=True);
        # Status update

        globs, bl_threshold, pared_tree, pared_branches, pruned_tips, over_max_spec = PARE.pare(globs, cur_tree, iteration);
        # Call the paring algorithm for the current tree

        if over_max_spec:
            step_start_time = CORE.report_step(globs, step, step_start_time, "Truncated", full_update=True);
            # Status update      

            CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: Branch length threshold for iteration:\t" + str(bl_threshold));
            CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: This iteration would remove " + str(len(pruned_tips)) + " tips, which puts the total number pruned over the maximum limit (" + str(globs['max-spec']) +  "). Paring complete.");
            # Info update

            pare = False;

        elif len(pruned_tips) == 0:
            step_start_time = CORE.report_step(globs, step, step_start_time, "Success: No branches pared this iteration. Exiting.", full_update=True);
            # Status update              

            CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: Branch length threshold for iteration:\t" + str(bl_threshold));
            # Info update

            pare = False;

        else: 
            step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(len(pruned_tips)) + " tips removed, " + str(len(pared_branches)) + " branches pared.", full_update=True);
            # Status update

            ####################

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

            with(open(os.path.join(globs['outdir'], "iter-" + str(iteration) + "-pared.tre"), "w")) as pared_tree_file:
                pared_tree_file.write(pared_tree);
            # Write the pared tree from this iteration to a file  

            step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
            # Status update

            ####################

            cur_tree = pared_tree;
            # Set the tree for the next iteration to be the pared tree from this iteration

            globs['bl-thresholds'].append(bl_threshold);
            globs['pared-branches'].append(pared_branches);
            globs['total-pared-branches'] += len(pared_branches);
            globs['pruned-tips'].append(pruned_tips);
            globs['total-pruned-tips'] += len(pruned_tips);
            globs['iter-trees'].append(pared_tree);
            # Add the pared branches and pruned tips from this iteration to the global lists for all iterations

            if iteration == globs['max-iterations']:
                iteration += 1;
                pare = False;
            # Check for stopping conditions to set pare to False
    # The paring loop

    ####################

    CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 50);
    step = "Writing summary stats to log";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update    

    headers = ["iteration", "branch length threshold", "branches pared", "tips pruned", "tree"];
    CORE.printWrite(globs['logfilename'], 3, "\t".join(headers));
    for i in range(iteration):
        
        iter_tree = globs['iter-trees'][i] + ";";
        
        #iter_tree = re.sub("<[\d]+>_", "", globs['iter-trees'][i]);
        #iter_tree = re.sub("<[\d]+>", "", iter_tree) + ";";
        # Removes the internal node labels from the tree

        if i == 0:
            outline = [ "0", "NA", "NA", "NA", iter_tree ];
        else:
            iter_str = str(i);
            outline = [ iter_str, str(globs['bl-thresholds'][i-1]), str(len(globs['pared-branches'][i-1])), str(len(globs['pruned-tips'][i-1])), iter_tree ];
        CORE.printWrite(globs['logfilename'], 3, "\t".join(outline));

    tips_file = os.path.join(globs['outdir'], "all-pruned-tips.txt");
    with open(tips_file, "w") as tipsfile:
        for iter_tips in globs['pruned-tips']:
            for tip in iter_tips:
                tipsfile.write(tip + "\n");

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
    # Status update

    ####################

    final_msg = "# Pared " + str(globs['total-pared-branches']) + " total branches and removed " + str(globs['total-pruned-tips']) + " total tips.";
    CORE.endProg(globs, final_msg);
    # A nice way to end the program

#############################################################################

