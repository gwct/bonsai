#############################################################################
# Functions to read and handle phylogenetic trees
# Gregg Thomas
#############################################################################

import sys
import os
import re
import random
import numpy as np
import lib.core as CORE

#############################################################################
def treeParse(tree, debug=0):
# The treeParse function takes as input a rooted phylogenetic tree with branch lengths and returns the tree with node labels and a
# dictionary with usable info about the tree in the following format:
# New (current) format:
# node:[branch length (if present), ancestral node, node type, node label (if present)]

    tree = tree.strip();
    if tree[-1] != ";":
        tree += ";";
    # Some string handling

    nodes, bl, supports, ancs = {}, {}, {}, {};
    # Initialization of all the tracker dicts

    topology = remBranchLength(tree);

    if debug == 1:
        print("TOPOLOGY:", topology);

    nodes = {};
    for n in topology.replace("(","").replace(")","").replace(";","").split(","):
        nodes[n] = 'tip';
    # nodes = { n : 'tip' for n in topology.replace("(","").replace(")","").replace(";","").split(",") };
    # Retrieval of the tip labels

    if debug == 1:
        print("NODES:", nodes);

    new_tree = "";
    z = 0;
    numnodes = 1;
    while z < (len(tree)-1):
        new_tree += tree[z];
        if tree[z] == ")":
            node_label = "<" + str(numnodes) + ">";
            new_tree += node_label;
            nodes[node_label] = 'internal';
            numnodes += 1;
        z += 1;
    nodes[node_label] = 'root';
    rootnode = node_label;
    # This labels the original tree as new_tree and stores the nodes and their types in the nodes dict

    if debug == 1:
        print("NEW TREE:", new_tree);
        print("TREE:", tree);
        print("NODES:", nodes);
        print("ROOTNODE:", rootnode);
        #sys.exit();
    topo = "";
    z = 0;
    numnodes = 1;
    while z < (len(topology)-1):
        topo += topology[z];
        if topology[z] == ")":
            node_label = "<" + str(numnodes) + ">";
            topo += node_label;
            numnodes += 1;
        z += 1;
    # This labels the topology with the same internal labels

    if debug == 1:
        print("TOPO:", topo);
        print("----------");
        print("TOPOLOGY:", topo);

    for node in nodes:
        if node + node in new_tree:
            new_tree = new_tree.replace(node + node, node);

    # if debug == 1:
    #     print new_tree;
    #     sys.exit();

    for node in nodes:
    # One loop through the nodes to retrieve all other info
        if debug == 1:
            print("NODE:", node);

        if nodes[node] == 'tip':
            supports[node] = "NA";
            if node + ":" in tree:
                cur_bl = re.findall(node + ":[\d.Ee-]+", new_tree);
                cur_bl = cur_bl[0].replace(node + ":", "");
                if debug == 1:
                    print("FOUND BL:", cur_bl);
                bl[node] = cur_bl;                
            else:
                bl[node] = "NA";

        elif nodes[node] == 'internal':
            if node + node in new_tree:
                new_tree = new_tree.replace(node + node, node);

            if node + "(" in new_tree or node + "," in new_tree or node + ")" in new_tree:
                if debug == 1:
                    print("NO BL OR LABEL");
                supports[node] = "NA";
                bl[node] = "NA";

            elif node + ":" in new_tree:
                supports[node] = "NA";
                cur_bl = re.findall(node + ":[\d.Ee-]+", new_tree);
                cur_bl = cur_bl[0].replace(node + ":", "");
                if debug == 1:
                    print("FOUND BL:", cur_bl);
                bl[node] = cur_bl;                                

            else:
                cur_bsl = re.findall(node + "[\d\w<>_*+.Ee/-]+:[\d.Ee-]+", new_tree);
                if cur_bsl:
                # If the pattern above is found then the node has both support and branch length
                    cur_bs = cur_bsl[0].replace(node, "");
                    cur_bs = cur_bs[:cur_bs.index(":")];
                    cur_bl = cur_bsl[0].replace(node, "").replace(cur_bs, "").replace(":", "");
                    if debug == 1:
                        print("FOUND BL AND LABEL:", cur_bl, cur_bs);
                    supports[node] = cur_bs;
                    bl[node] = cur_bl;
                    #new_tree = new_tree.replace(cur_bs, "");
                else:
                # If it is not found then the branch only has a label
                    cur_bs = re.findall(node + "[\w*+.<> -]+", new_tree);
                    cur_bs = cur_bs[0].replace(node, "");
                    if debug == 1:
                        print("FOUND LABEL:", cur_bs);
                    supports[node] = cur_bs;
                    bl[node] = "NA";
                    #new_tree = new_tree.replace(cur_bs, "");

        elif nodes[node] == 'root':
            bl[node] = "NA";
            supports[node] = new_tree[new_tree.index(node)+len(node):];
            ancs[node] = "NA";
            continue;

        # Next we get the ancestral nodes. If the node is the root this is set to NA.
        anc_match = re.findall('[(),]' + node, new_tree);

        #if nodes[node] == 'internal':
        #    sys.exit();
        # anc_match = re.findall(node + '[\d:(),]+', new_tree);
        anc_match = re.findall(node, topo);
        if debug == 1:
            print("ANC MATCH:", anc_match);

        anc_tree = new_tree[new_tree.index(anc_match[0]):][1:];
        # Ancestral labels are always to the right of the node label in the text of the tree, so we start our scan from the node label

        if debug == 1:
            print("NODE:", node);
            print("ANC_MATCH:", anc_match);
            print("ANC_TREE:", anc_tree);
            
        cpar_count = 0;
        cpar_need = 1;

        for i in range(len(anc_tree)):
        # We find the ancestral label by finding the ) which matches the nesting of the number of ('s found
            if anc_tree[i] == "(":
                cpar_need = cpar_need + 1;
            if anc_tree[i] == ")" and cpar_need != cpar_count:
                cpar_count = cpar_count + 1;
            if anc_tree[i] == ")" and cpar_need == cpar_count:
                anc_tree = anc_tree[i+1:];
                ancs[node] = anc_tree[:anc_tree.index(">")+1];
                break;

        if debug == 1:
            print("FOUND ANC:", ancs[node]);
            print("---");
    nofo = {};
    for node in nodes:
        nofo[node] = [bl[node], ancs[node], nodes[node], supports[node]];
    # Now we just restructure everything to the old format for legacy support

    if debug == 1:
    # Debugging options to print things out
        print(("\ntree:\n" + tree + "\n"));
        print(("new_tree:\n" + new_tree + "\n"));
        print(("topology:\n" + topo + "\n"));
        print("nodes:");
        print(nodes);
        print()
        print("bl:");
        print(bl);
        print()
        print("supports:");
        print(supports);
        print()
        print("ancs:");
        print(ancs);
        print()
        print("-----------------------------------");
        print()
        print("nofo:");
        print(nofo);
        print()

    return nofo, topo, rootnode;

#############################################################################

def remBranchLength(treestring):
# Removes branch lengths from a tree.

    treestring = re.sub('[)][\d\w<>/.eE_:-]+', ')', treestring);
    treestring = re.sub(':[\d.eE-]+', '', treestring);

    return treestring;

#############################################################################

def getDesc(d_node, d_treedict):
# This function takes a species in the current tree and the dictionary of the current tree
# returned by treeparse and finds the direct descendants of the species.
    d_list = [];
    for node in d_treedict:
        if d_treedict[node][1] == d_node:
            d_list.append(node);

    if d_list == []:
        return [d_node];
    else:
        return d_list;

#############################################################################

def getClade(c_node, c_treedict):
# This function takes a species in the current tree and the dictionary of the current tree
# returned by treeparse and finds all tip labels that are descendants of the current node.
# This is done by getting the direct descendants of the current node with getDesc and then
# recursively calling itself on those descendants.
    clade = [];
    c_desc = getDesc(c_node, c_treedict);
    for d in c_desc:
        if c_treedict[d][2] != 'tip':
            clade.append(getClade(d, c_treedict));
        else:
            clade.append(d);

    r_clade = [];
    for c in clade:
        if type(c) == list:
            for cc in c:
                r_clade.append(cc);
        else:
            r_clade.append(c);

    return r_clade;

#############################################################################

def addBranchLength(tree, treedict):
# Re-writes the branch lengths onto a topology parsed by treeParse.

    for node in treedict:
        if treedict[node][2] == 'root':
            continue;
        # No branch length on the root, so skip

        if treedict[node][0] != "NA":
            if treedict[node][3] != "NA":
                tree = tree.replace(node, node + "_" + str(treedict[node][3]) + ":" + str(treedict[node][0]));
            else:
                tree = tree.replace(node, node + ":" + str(treedict[node][0]));
        # If the node has a branch length, check if it has a label and add both or just the bl

        elif treedict[node][3] != "NA":
            tree = tree.replace(node, node + "_" + str(treedict[node][3]));
        # If the node has no branch length, check if it has a label and add it back

    return tree;

#############################################################################

def LCA(spec_list, treedict):
# Given a list of nodes, this function finds the least common ancestor of them,
# and tells whether the nodes provided form a monophyletic clade.

	ancs = {};
	for spec in spec_list:
		ancs[spec] = [spec];

	for spec in spec_list:
		if treedict[spec][2] == 'root':
			continue;
		curanc = treedict[spec][1];
		ancs[spec].append(curanc);
		while treedict[curanc][2] != 'root':
			curanc = treedict[curanc][1];
			ancs[spec].append(curanc);

	intersect_anc = set.intersection(*list(map(set, list(ancs.values()))))
	lcp = [t for t in list(ancs.values())[0] if t in intersect_anc]

	#lcp = sorted(set.intersection(*map(set, ancs.values())), key=lambda x: ancs.values()[0].index(x))
	monophyletic = False;
	if set(getClade(lcp[0],treedict)) == set(spec_list):
		monophyletic = True;

	return lcp[0], monophyletic;

#############################################################################

def getSubtree(node, tree):
# Gets the subtree string at a given node from a labeled tree string from treeParse

    subtree = "";
    # Initialize the subtree string

    partree = tree[:tree.index(node)][::-1]
    # Slice the tree string at the index of the node label and reverse it

    cp = 0;
    op = 0;
    # Counts of closing an opening parentheses. The subtree will be complete when they
    # are equal

    for c in partree:
        if c == ")":
            cp = cp + 1;
        if c == "(":
            op = op + 1;
        subtree = subtree + c;
        if cp == op:
            break;
    # Loop through every character in the sliced and reversed tree, counting parentheses and adding
    # charcters to the subtree string one at a time. Stop when the number of closing parentheses equals
    # the number of opening

    return subtree[::-1];
    # Return the reverse of the subtree, which is the original orientation

#############################################################################

def adjustTreeDict(tree_dict, root):
    
    num_no_supp = 0;
    # Number of nodes with no support

    bls = [];
    # A list of all branch lengths in the tree

    na_bl_spec = [];
    # A list of species that have no branch length (i.e. tips from Astral)

    for n in tree_dict:
        if n == root:
            continue;
        # Skip the root since it won't have a length or support

        if tree_dict[n][0] == "NA":
            tree_dict[n][0] = 0.0;
            na_bl_spec.append(n);
        # If the current node has no branch length, set it as 0.0 in the dictionary for later and
        # add to the na_bl_spec dict to add the NA back in later. Also don't add this to the list
        # of bls to calculate the percentile from
        else:
            tree_dict[n][0] = float(tree_dict[n][0]);
            bls.append(tree_dict[n][0]);
            # Add the floated bl to the bl list to calculate percentile. The NA branch lengths should
            # be left out
        # Convert the branch lengths in the tree dict to floats

        if tree_dict[n][3] == 'NA':
            num_no_supp += 1;
            continue;
        cur_supp = tree_dict[n][3];
        if "/" in cur_supp:
            cur_supp = cur_supp.split("/")[1];
        cur_supp = cur_supp.replace("_", "");
        tree_dict[n][3] = float(cur_supp);
        # Parse out the gCF from the supports, skip if there is no support
        # TODO: This is still very specific to my input format

    return tree_dict, bls, na_bl_spec;

#############################################################################

def pare(globs, cur_tree, iteration, debug=False):
    #debug=True;
    if debug:
        print();
        print("CUR TREE: ", cur_tree);

    cur_tree_dict, cur_labeled_tree, cur_root = treeParse(cur_tree);
    if debug:
        print("CUR TREE DICT: ", cur_tree_dict)
    cur_tree_dict, bls, na_bl_spec = adjustTreeDict(cur_tree_dict, cur_root);
    if debug:
        print("CUR TREE DICT: ", cur_tree_dict)
    # Parse the tree string and do some extra manipulation of the tree dict

    cur_bl_tree = addBranchLength(cur_labeled_tree, cur_tree_dict);
    if debug:
        print("CUR BL TREE: ", cur_bl_tree);
        print("--------");
    # Add branch lengths back to the labeled tree string

    if iteration == 1:
        no_tip_bl_tree = cur_bl_tree;
        for n in cur_tree_dict:
            if cur_tree_dict[n][2] != "tip":
                continue;
            no_tip_bl_tree = re.sub(n + ":0.0[\)]", n + ")", no_tip_bl_tree);
            no_tip_bl_tree = re.sub(n + ":0.0,", n + ",", no_tip_bl_tree);
        globs['iter-trees'].append(no_tip_bl_tree);
    # For the first iteration, add the original tree to the global list here so it
    # is the same format as the others

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

    criteria = "min-spec";
    if criteria not in ["deterministic", "random", "min-spec", "max-gcf"]:
        sys.exit("\n\nINTERNAL: Invalid clade selection criteria specified: " + criteria + "\n\n");
    # One of: deterministic, random, min-spec, max-gcf

    for n in cur_tree_dict:
        if cur_tree_dict[n][2] in ['root', 'tip']:
            continue;
        # Skip tips and the root when getting branches to pare
        ## TODO: Consider tips for non-Astral trees?

        if cur_tree_dict[n][0] < cur_bl_threshold and cur_tree_dict[n][3] < globs['gcf-threshold']:
        # If the branch length of the current node is lower than both the bl threshold and the gcf threshold, select
        # species to pare

            cur_desc = getDesc(n, cur_tree_dict);
            # Get the two descendant nodes of the node to pare

            clade1 = getClade(cur_desc[0], cur_tree_dict);
            clade2 = getClade(cur_desc[1], cur_tree_dict);
            # Get all species in the clades descending from the current descendants

            if criteria == "deterministic":
                cur_prune_clade = clade1;
                cur_prune_node = cur_desc[0];
                prune_index = 0;
                # Variables to keep track of the pruned branch
            # Deterministic pruning (good for debugging)

            ##########

            elif criteria == "random":
                cur_prune_clade = random.choice([clade1, clade2]);
                
                if cur_prune_clade == clade1:
                    cur_prune_node = cur_desc[0];
                    prune_index = 0;
                else:
                    cur_prune_index = cur_desc[1];
                    prune_index = 1;
                # Variables to keep track of the pruned branch
            # Randomly choose one of the clades to pare

            ##########

            elif criteria == "min-spec":
                if len(clade1) <= len(clade2):
                    cur_prune_clade = clade1;
                    cur_prune_node = cur_desc[0];
                    prune_index = 0;
                else:
                    cur_prune_clade = clade2;
                    cur_prune_node = cur_desc[1];
                    prune_index = 1;
                # Variables to keep track of the pruned branch
            # Minimize species removed

            ##########

            elif criteria == "max-gcf":
                gcfs = [];
                # We will generate two lists of gcfs to average: one from each descendant node

                for d in [0,1]:
                # Loop through both descendant nodes

                    if debug:
                        print("CUR DESC: ", cur_desc[d]);
                        print(cur_tree_dict[cur_desc[d]]);
                    # Debug statements

                    if cur_tree_dict[cur_desc[d]][2] != 'tip':
                    # If the tree isn't a tip, get gcfs

                        cur_gcfs = [];
                        # The list of subtrees for the current descendant node

                        cur_subtree = getSubtree(cur_desc[d], cur_bl_tree);
                        cur_subtree = re.sub("<[\d]+>_", "", cur_subtree);
                        if debug:
                            print("CUR SUBTREE: ", cur_subtree);
                        # Get the subtree descending from the current descendant node and remove
                        # the previous labels

                        subtree_dict, labeled_subtree, subtree_root = treeParse(cur_subtree);
                        # Parse the subtree

                        for sub_n in subtree_dict:
                            if subtree_dict[sub_n][3] not in ["NA", ""]:
                                cur_gcfs.append(float(subtree_dict[sub_n][3]));
                        # For every node in the subtree, if there is a gcf add it to the current list of gcfs

                        cur_gcfs.append(cur_tree_dict[cur_desc[d]][3]);
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
                for gcf in gcfs:
                    avg_gcfs.append(sum(gcf) / len(gcf));
                if debug:
                    print("AVG GCFS: ", avg_gcfs);
                # Average the gcfs from both subtrees

                if avg_gcfs[0] <= avg_gcfs[1]:
                    if debug:
                        print("HI")
                    cur_prune_clade = clade1;
                    cur_prune_node = cur_desc[0];
                    prune_index = 0;
                else:
                    cur_prune_clade = clade2;
                    cur_prune_node = cur_desc[1];
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
                if debug:
                    print(num_spec_to_pare);
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

    num_pruned = 0;
    #to_pare = ["Mallomys_rothschildi_ABTC47402", "Abeomelomys_sevia_KUM161018"]
    for n in spec_to_prune:
        if debug:
            #n = "Zyzomys_pedunculatus_Z34925";
            #n = "Mylomys_dybowskii_MNHN1997072";
            #n = "Abeomelomys_sevia_KUM161018";
            # if iteration == 2:
            #     n = "Bunomys_chrysocomus_JAE4867";
            print("SPEC TO PRUNE:", n);
            print(cur_tree_dict[n])
            print("--------");
        
        if n not in cur_tree_dict:
            print(n);
            print(cur_labeled_tree);
            print(cur_tree_dict);
            sys.exit("SPECIES NOT FOUND IN DICT");
        # A check... need to make better

        anc = cur_tree_dict[n][1];
        anc_supp = cur_tree_dict[anc][3];
        anc_bl = cur_tree_dict[anc][0];
        if debug:
            print("ANCESTOR OF NODE TO PARE: ", anc);
            print(cur_tree_dict[anc]);
            print("--------");
        # Get info about the ancestor of the species to prune

        anc_desc = getDesc(anc, cur_tree_dict);
        sis = [ d for d in anc_desc if d != n ][0];
        if debug:
            print("SISTER OF NODE TO PARE: ", sis);
            print(cur_tree_dict[sis]);
            print("--------");
        sis_bl = cur_tree_dict[sis][0];
        if sis_bl == "NA":
            sis_bl = 0.0;
        sis_supp = cur_tree_dict[sis][3];
        # Get info about the sister of the species to prune by checking the
        # other descendant of the ancestor

        cur_subtree = getSubtree(anc, cur_bl_tree);
        cur_subtree += anc + "_" + str(anc_supp) + ":" + str(anc_bl);
        if debug:
            print("SUBTREE AT ANCESTOR: ", cur_subtree);
            print("--------");
        # Get the subtree of the ancestor and ad the ancestral node to it

        sis_subtree = getSubtree(sis, cur_bl_tree);
        if debug:
            print("SUBTREE AT SISTER: ", sis_subtree);
            print("--------");
        # Get the subtree of the sister

        new_bl = anc_bl + sis_bl;
        # After removing the current species, the new branch length of the combined anc and sis
        # branch will be their sum

        if cur_tree_dict[sis][2] == 'tip':
            new_subtree = sis + ":" + str(new_bl);
        else:
            new_subtree = sis_subtree + sis + "_" + str(sis_supp) + ":" + str(new_bl);
        if debug:
            print("NEW SUBTREE: ", new_subtree);
            print("--------");
        # Make the new subtree to replace the current ancestral one. If the sister is
        # a tip, then this is just the sister label with the new branch length.
        # If the sister is a clade, then this is the sister subtree with the sister added
        # on with the new branch length

        new_tree = cur_bl_tree.replace(cur_subtree, new_subtree);
        new_tree = re.sub("<[\d]+>_", "", new_tree);
        new_tree = re.sub("<[\d]+>", "", new_tree);
        if debug:
            print("NEW TREE: " + new_tree);
            print("--------");
        # Replace the old subtree with the new one and remove the internal node labels from treeParse

        cur_tree_dict, cur_labeled_tree, cur_root = treeParse(new_tree);
        cur_tree_dict, bls, na = adjustTreeDict(cur_tree_dict, cur_root);
        cur_bl_tree = addBranchLength(cur_labeled_tree, cur_tree_dict);
        if debug:
            print("NEW PARSED BL TREE: ", cur_bl_tree);
            print("--------");
        # Parse the new tree and add the branch lengths back on

        num_pruned += 1;
        # Increment the number of tips pruned
    # Tip pruning block
    ####################

    if debug:
        print("PARED TREE: ", cur_bl_tree);
        print("PARED TREE DICT: ", cur_tree_dict);
        print("--------");

        print("TIPS IN ORIG TERE: ", len(globs['tips']));
        print("TIPS EXPECTED TO PRUNE: ", len(spec_to_prune));
        print("ACTUAL TIPS PRUNED: ", num_pruned);
        final_tips = [ n for n in cur_tree_dict if cur_tree_dict[n][2] == 'tip' ];
        print("TIPS IN FINAL TREE: ", len(final_tips));
        print("--------");

        for n in globs['tips']:
            if n not in spec_to_prune and n not in final_tips:
                print("Should NOT have been pared: " + n);

            if n in spec_to_prune and n in final_tips:
                print("Should have been pared: " + n);
    # Some debug statements
    ####################
        #print(spec_to_prune);

    for spec in na_bl_spec:
        cur_bl_tree = cur_bl_tree.replace(spec + ":0.0,", spec + ",");
        cur_bl_tree = cur_bl_tree.replace(spec + ":0.0)", spec + ")");
    if debug:
        print("CUR BL TREE: ", cur_bl_tree);
    # For the next iteration, we need to remove those 0.0 branch lengths added in so
    # they aren't counted in the distribution when calculating percentiles
    ####################

    return globs, cur_bl_threshold, cur_bl_tree, branches_to_pare, spec_to_prune, False;
    
#############################################################################