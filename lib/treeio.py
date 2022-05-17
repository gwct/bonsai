#############################################################################
# Handles the reading of the input trees
# Gregg Thomas
#############################################################################

import sys
import lib.core as CORE
import lib.tree as TREE

#############################################################################

def readST(globs):

    step = "Reading input species tree";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    if globs['st-input-type'] == "file":
        globs['orig-st-str'] = open(globs['st-input'], "r").read().strip();
    else:
        globs['orig-st-str'] = globs['st-input']
    # If the input type is a file, read the file here, otherwise take the string as the tree input

    globs['st'] = TREE.Tree(globs['orig-st-str']);
    # Read the tree

    root_str = "unrooted";
    add_root_node = 0;
    if globs['st'].rooted:
        root_str = "rooted";
        add_root_node = 1;
    # Node counting for log

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + root_str + " tree read");
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: Tree has "+ str(globs['st'].num_tips) + " tips and " + str(globs['st'].num_internals + add_root_node) + " internal nodes");
    # Status updates

    return globs;

#############################################################################

def readGT(globs):

    step = "Reading input gene trees";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update    

    line_num, tree_num, num_rooted, num_unrooted = 0, 0, 0, 0;
    for line in open(globs['gt-input']):
        line = line.strip();
        line_num += 1;

        if not line:
            globs['gt-input-empty'] += 1;
            continue;

        try:
            tree_num += 1;
            globs['gts'][tree_num] = TREE.Tree(line);
        except:
           globs['warnings'] += 1;
           globs['gt-input-skipped'] += 1;
           CORE.printWrite(globs['logfilename'], -2, "# WARNING: Could not read line " + str(line_num) + " of gene tree file (-g) as a tree! Skipping.");
           # Warning statement if a line can't be read as a tree

        if globs['gts'][tree_num].rooted:
            num_rooted += 1;
        else:
            num_unrooted += 1;

    if globs['gts']:
        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(len(globs['gts'])) + " gene trees read");
        CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: " + str(num_rooted) + " trees are rooted and " + str(num_unrooted) + " are unrooted");
        if globs['gt-input-empty']:
            CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: Skipped " + str(globs['gt-input-empty']) + " lines because they were empty. See log for more info.");
        if globs['gt-input-skipped']:
            CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: Skipped " + str(globs['gt-input-skipped']) + " lines because they couldn't be read as trees. See log for more info.");
    else:
        step_start_time = CORE.report_step(globs, step, step_start_time, "Failure: no gene trees read");
        print("\n\n");
        CORE.errorOut("2", "Error reading gene trees! Couldn't parse any trees from file. Make sure they are formatted as rooted, Newick trees.", globs);
    # Status update

    return globs, num_rooted, num_unrooted;

#############################################################################

def readBranches(globs, branch_file):

    branches, clades = [], [];

    for line in open(branch_file):
        # Every line in the file corresponds to one branch

            if line.startswith("#"):
                continue;
            # Skip lines that are commented out

            line = line.strip();
            # Remove trailing whitespace

            if " " in line:
            # If there is a space in the line, assume it is defining a branch with 2 tip labels

                specs = line.split(" ");
                # Split the line by the two tips

                if not all(s in globs['st'].tips for s in specs):
                    globs['warnings'] += 1;
                    CORE.printWrite(globs['logfilename'], globs['log-v'], "# WARNING: Label not found in tree. Skipping line: " + line);
                # Throw a warning if both of the tips aren't found in the tree and skip

                else:
                    exempt_node = globs['st'].LCA(specs);
                    # Get the least common ancestor of the given tips

                    CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: " + line + " -> " + exempt_node);
                    # Info statement to show branch assignment

                    branches.append(exempt_node);
                    # Add the internal branch to the list of exempt branches

                    clades.append(set(globs['st'].getClade(exempt_node)));
                    # Get and add the full clade descending from the internal branch to the list of exempt clades
                # If both tips are found, get the full clade
            ## End tip block

            else:
            # If there is no space, assume it is defining an internal branch by label

                if line not in globs['st'].nodes:
                    globs['warnings'] += 1;
                    CORE.printWrite(globs['logfilename'], globs['log-v'], "# WARNING: Label in exempt file not found in tree. Skipping line: " + line);
                # Print a warning if the branch label isn't found in the tree

                else:
                    branches.append(line);
                    clades.append(set(globs['st'].getClade(line)));
                # Otherwise, add the branch and its clade to their global lists
            ## End label block
        ## End line/branch loop   

    return globs, branches, clades;

#############################################################################

def writeCF(cur_tree, gcf_dict, scf_dict, outfilename, treefilename, globs):

    step = "Writing CF files"
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    gcf_headers = ["ID", "gCF", "gCF_N", "gDF1", "gDF1_N", "gDF2", "gDF2_N", "gDFP", "gDFP_N", "gN", "gM" ];
    scf_headers = [ "sCF", "sCF_N", "sDF1", "sDF1_N", "sDF2", "sDF2_N", "sN" ];

    cf_headers = gcf_headers;
    if scf_dict:
        cf_headers = gcf_headers + scf_headers;
    cf_headers += [ "label", "length" ];

    gcf_labels = {};
    scf_labels = {};
    #print();

    with open(outfilename, "w") as outfile:
        outfile.write("# Gene concordance factor statistics\n");
        outfile.write("# ID:        Branch ID in the species tree\n");
        outfile.write("# gCF:       Gene concordance factor (=gCF_N/gN)\n");
        outfile.write("# gCF_N:     The number of trees concordant with this branch\n");
        outfile.write("# gDF1:      Gene discordance factor for NNI-1 (=gDF1_N/gN)\n");
        outfile.write("# gDF1_N:    Number of trees concordant with NNI-1 topology\n");
        outfile.write("# gDF2:      Gene discordance factor for NNI-2 (=gDF2_N/gN)\n");
        outfile.write("# gDF2_N:    Number of trees concordant with NNI-2 topology\n");
        outfile.write("# gDFP:      Gene discordance factor due to polyphyly (=gDFP_N/gN)\n");
        outfile.write("# gDFP_N:    Number of trees decisive but discordant due to polyphyly\n");
        outfile.write("# gN:         The number of trees decisive for this branch\n");
        outfile.write("# gM:         The number of trees with one or more groups missing from this branch\n");

        if scf_dict:
            outfile.write("# sCF:       Site concordance factor averaged over qN quartets (=sCF_N/sN)\n");
            outfile.write("# sCF_N:     sCF in absolute number of sites\n");
            outfile.write("# sDF1:      Site discordance factor for alternative quartet 1 (=sDF1_N/sN)\n");
            outfile.write("# sDF1_N:    sDF1 in absolute number of sites\n");
            outfile.write("# sDF2:      Site discordance factor for alternative quartet 2 (=sDF2_N/sN)\n");
            outfile.write("# sDF2_N:    sDF2 in absolute number of sites\n");
            outfile.write("# sN:        Number of decisive sites averaged over 100 quartets\n");
            outfile.write("# qN:        Number of quartets sampled with at least 1 decisive site on this branch\n");

        outfile.write("# label:     The label in the original species tree string\n");
        outfile.write("# length:    The branch length in the original species tree string\n");
        outfile.write("# NOTE: (gCF+gDF1+gDF2+gDFP) = 1.0 and (gCF_N+gDF1_N+gDF2_N+gDFP_N) = N\n");
        outfile.write("# NOTE: (N+M) = total number of trees\n");

        outfile.write("\t".join(cf_headers) + "\n");
        for node in gcf_dict:
            outline = [];

            if gcf_dict[node]["decisive"] == 0:
                outline = [node] + [ "NA" for h in gcf_headers ];
                outline[-2] = str(gcf_dict[node]["decisive"]);
                outline[-1] = str(gcf_dict[node]["nondecisive"])
                gcf_labels[node] = "NA";
            
            else:
                gcf = round(gcf_dict[node]["concordant"] / gcf_dict[node]["decisive"], 3);
                outline.append(gcf);
                outline.append(gcf_dict[node]["concordant"]);
                gcf_labels[node] = str(gcf);

                outline.append(round(gcf_dict[node]["disco1"] / gcf_dict[node]["decisive"], 3));
                outline.append(gcf_dict[node]["disco1"]);

                outline.append(round(gcf_dict[node]["disco2"] / gcf_dict[node]["decisive"], 3));
                outline.append(gcf_dict[node]["disco2"]);

                outline.append(round(gcf_dict[node]["para"] / gcf_dict[node]["decisive"], 3));
                outline.append(gcf_dict[node]["para"]);

                outline.append(gcf_dict[node]["decisive"]);
                outline.append(gcf_dict[node]["nondecisive"]);

            if scf_dict:
                if node in scf_dict:

                    outline.append(round(scf_dict[node]['scf'], 3));
                    outline.append(round(scf_dict[node]['concordant-sites'], 2));

                    outline.append(round(scf_dict[node]['sdf1'], 3));
                    outline.append(round(scf_dict[node]['disco1-sites'], 2));

                    outline.append(round(scf_dict[node]['sdf2'], 3));
                    outline.append(round(scf_dict[node]['disco2-sites'], 2));                

                    outline.append(round(scf_dict[node]['decisive-sites'], 2));

                    scf_labels[node] = str(round(scf_dict[node]['scf'], 3));

                    outline.append(scf_dict[node]['total-quartets']);

                else:
                    outline += ["NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"];
                    scf_labels[node] = "NA";

            outline = [node] + [ str(v) for v in outline ] + [cur_tree.label[node], cur_tree.bl[node]];
            outfile.write("\t".join(outline) + "\n");
            #print("\t".join(outline));
            
    if scf_dict:
        new_labels = { node : gcf_labels[node] + "/" + scf_labels[node] for node in gcf_dict };
        cf_tree_str = cur_tree.addLabel(new_labels);
        gcf_tree_str = cur_tree.addLabel(gcf_labels);
    else:
        cf_tree_str = cur_tree.addLabel(gcf_labels);
        gcf_tree_str = cf_tree_str;

    with open(treefilename, "w") as cf_tree_file:
        cf_tree_file.write(cf_tree_str);

    gcf_tree = TREE.Tree(gcf_tree_str);

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
    # Status updated

    return gcf_tree;

#############################################################################

def writeTopoCounts(globs, topo_counts, num_trees, topo_count_file):

    step = "Writing topo counts"
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update


    with open(topo_count_file, "w") as outfile:
        outfile.write("# Topology counts\n");
        outfile.write("# " + str(num_trees) + " trees counted\n");
        outfile.write("# " + str(len(topo_counts["topo"])) + " unique topologies found\n");
        outfile.write("# Reference tree (from -st): " + globs['st'].topo_str + "\n");
        outfile.write("# count:             Number of times this topology was found in the input set\n");
        outfile.write("# rf.to.ref.tree:    The Robinson-Foulds distance from this topology to the reference topology\n");
        outfile.write("# total.splits:      The sum of the total possible splits in this topology and the reference topology\n");
        outfile.write("# topology:          The Newick formatted string of the topology\n");
        outfile.write("# NOTE: rf.to.ref.tree / total.splits = normalized Robinson-Foulds distance\n");

        headers = ["count", "rf.to.ref.tree", "total.splits", "topology", ];
        outfile.write("\t".join(headers) + "\n");

        topo_lists = zip(topo_counts["count"], topo_counts["rf"], topo_counts["splits"], topo_counts["tree"]);
        topo_lists = sorted(topo_lists, reverse=True);
        # Also reverse sorts rf...

        #topo_lists = zip(*sorted(zip(topo_counts["count"], topo_counts["rf"], topo_counts["splits"], topo_counts["tree"]),reverse=True))

        for row in topo_lists:
            outline = [ str(row[0]), str(row[1]), str(row[2]), row[3] ];
            outfile.write("\t".join(outline) + "\n")

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
    # Status updated

#############################################################################

def writeGT(globs, gts, gt_file):

    step = "Writing gene trees"
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    with open(gt_file, "w") as outfile:
        for gt_id in sorted(list(gts.keys())):
            if gts[gt_id] == "NA":
                outfile.write("NA\n");
            else:
                outfile.write(gts[gt_id].orig_tree_str + "\n");

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
    # Status updated

    return globs;

#############################################################################