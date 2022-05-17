            #iter_tree = re.sub("<[\d]+>_", "", globs['iter-trees'][i]);
            #iter_tree = re.sub("<[\d]+>", "", iter_tree) + ";";
            # Removes the internal node labels from the tree

#############################################################################

def gcf(globs, st, gts):

    node_counts = { node : { "concordant" : 0.0, "disco1" : 0.0, "disco2" : 0.0, "para" : 0.0, "total" : 0.0 } for node in st.internals };
    gt_pruned, st_pruned = 0, 0;
    gt_skipped = 0;
    total_trees = 0.0;
    # Count vars

    gt_used = {};
    print();
    # if count_tops:
    #     tops, top_counts, top_trees = [], [], [];

    st_clades = st.getClades();
    st_splits = st.getSplits();
    st_quartets = st.getQuartets();
    # for node in st.internals:
    #     #print(node);
    #     st_quartets[node] = st.getQuartet(node);

    for gt_id in gts:
        print(gt_id);
        # print("orig st: ", cur_st.labeled_topo_str);
        cur_gt = gts[gt_id];
        # Get trees

        # if count_tops:
        #     gclade = set([frozenset(tp.getClade(node, ginfo)) for node in ginfo if ginfo[node][2] != 'tip']);
        #     if gclade in tops:
        #         topind = tops.index(gclade);
        #         top_counts[topind] += 1;
        #     else:
        #         tops.append(gclade);
        #         top_counts.append(1);
        #         top_trees.append(gtree);


        gt_tips_to_prune = set(cur_gt.tips) - set(st.tips);
        if cur_gt.num_tips < 4 or len(cur_gt.tips) - len(gt_tips_to_prune) < 4:
            gt_skipped += 1;
            globs['warnings'] += 1;
            CORE.printWrite(globs['logfilename'], -2, "# WARNING: Gene tree on line " + str(gt_id) + " has too few tips after pruning to match the species tree. Skipping.");
            continue;
        
        # Throw a warning if there aren't enough tips left after pruning the gene tree and skip

        if gt_tips_to_prune:
            #print("ya", gt_id, gt_tips_to_prune);
            gt_nodes_to_prune = cur_gt.findClades(gt_tips_to_prune);
            cur_gt = cur_gt.Prune2(gt_nodes_to_prune);
            cur_gt = TREE.Tree(cur_gt);
            #print("prune gt:", cur_gt.tree);
            CORE.printWrite(globs['logfilename'], -2, "# INFO: The gene tree on line " + str(gt_id) + " was pruned of the following " + str(len(gt_tips_to_prune)) + " species because they were not in the species tree: " + str(gt_tips_to_prune));
            gt_pruned += 1;
        # Prune gene tree to match tips in species tree if necessary



        st_tips_to_prune = set(st.tips) - set(cur_gt.tips);
        if st_tips_to_prune:
            # cur_st.showDesc();
            # print("+++++++++++++")
            # cur_st.showAnc();
            # print(st.labeled_topo_str);
            # print(st.num_tips);
            # # print(cur_st.desc["<30>"]);
            # print("==============================");
            # print("hey", gt_id, len(st_tips_to_prune), st_tips_to_prune);
            # print("==============================");
            # st_nodes_to_prune = list(cur_st.findClades(st_tips_to_prune));
            # print("hey agi ", gt_id, len(st_nodes_to_prune), st_nodes_to_prune);

            pruned_st_str = st.Prune2(list(st_tips_to_prune), debug=False);

            pruned_st = TREE.Tree(pruned_st_str);

            # cur_st.showDesc();
            # print("+++++++++++++")
            # cur_st.showAnc();
            # print("prune st:", pruned_st.labeled_topo_str);
            # print(pruned_st.num_tips);
            #print("----");

            assert pruned_st.num_tips == st.num_tips - len(st_tips_to_prune), "ERROR: INCORRECT NUMBER OF TIPS AFTER PRUNING (" + str(gt_id) + ")"; 

            CORE.printWrite(globs['logfilename'], -2, "# INFO: To compare with the gene tree on line " + str(gt_id) + " the species tree was pruned of the following " + str(len(st_tips_to_prune)) + " species: " + str(st_tips_to_prune));
            st_pruned += 1;


            pruned_st_clades = pruned_st.getClades();
            pruned_st_splits = pruned_st.getSplits();

            node_map, rev_map = TREE.mapNodes(pruned_st, pruned_st_clades, pruned_st_splits, st, st_clades, st_splits);

            # print(node_map);
            # print("==============================");
            # print(cur_gt.labeled_topo_str);
            # print("==============================");

            quartets = pruned_st.getQuartets();

        else:
            pruned_st = st;
            quartets = st_quartets;

            
            #print(node_map);
        # print("===================")
        # Prune species tree to match tips in gene tree if necessary

        # print(cur_gt.tree_str);
        # ts = "(((a,b),(c,d)),(e,f),(g,h))";
        # d1 = "((a,(c,d)),((b,f),e),(g,h))";
        # d2 = "((b,(c,d)),((a,f),e),(g,h))"
        # pruned_st = TREE.Tree(ts);
        # print(pruned_st.labeled_topo_str)
        # pruned_st.showQuartet();
        # # cur_gt = TREE.Tree(ts);
        # print("-----");
        # print(cur_gt.labeled_topo_str)
        # cur_gt.showQuartet();
        # print("-----");
        # ##print(TREE.mapNodes(cur_gt, pruned_st))

        # print(pruned_st.labeled_topo_str);
        # pruned_st.showDesc();
        # print("-----");
        # pruned_st.showAnc();
        # print("-----");
        # pruned_st.showSplit();
        # print("-----");
        # pruned_st.showQuartet();
        # for node in quartets:
        #     print(node, quartets[node]);
        # print("-----");
        # print(pruned_st.labeled_topo_str);
        # print(pruned_st.rooted);
        # print(pruned_st.root);
        
        # print("============================================================")
        # sys.exit();

        gt_clades = cur_gt.getClades();
        gt_splits = cur_gt.getSplits();
        gt_tips = set(cur_gt.tips);


        for node in pruned_st.internals:
            st_node = node;
            if st_tips_to_prune:
                st_node = node_map[node];

            if not all(st_quartets[st_node][q] & gt_tips for q in ["d2", "s", "q4"]):
                continue;

            # if not st_quartets[st_node]["d1"] & gt_tips:
            #     continue;
            # if not st_quartets[st_node]["d2"] & gt_tips:
            #     continue;            
            # if not st_quartets[st_node]["s"] & gt_tips:
            #     continue;
            # if not st_quartets[st_node]["q4"] & gt_tips:
            #     continue;

            possible_topos = { "concordant" : quartets[node]["d1"] | quartets[node]["d2"],
                                "disco1" : quartets[node]["d1"] | quartets[node]["s"],
                                "disco2" : quartets[node]["d1"] | quartets[node]["q4"] }            
            # print("gcf topos ", time.time() - start2);

            classification = "NA";

            # start2 = time.time()
            if cur_gt.findSplits(possible_topos["concordant"], gt_clades, gt_splits):
                classification = "concordant";
                node_counts[st_node]["concordant"] += 1;
            elif cur_gt.findSplits(possible_topos["disco1"], gt_clades, gt_splits):
                classification = "disco1";
                node_counts[st_node]["disco1"] += 1;
            elif cur_gt.findSplits(possible_topos["disco2"], gt_clades, gt_splits):
                classification = "disco2";
                node_counts[st_node]["disco2"] += 1;
            else:
                classification = "para";
                node_counts[st_node]["para"] += 1;

            node_counts[st_node]["total"] += 1;

            ### TODO: Condition for missing clade
                ## Check topos...?

        # For every clade in the species tree, check if it exists in the gene tree

        gt_used[gt_id] = cur_gt;

        total_trees += 1.0;
        
        # Iterate the total tree count

    ## End gene tree loop
    

    return globs, gt_used, node_counts, total_trees, gt_skipped, gt_pruned, st_pruned;

    # print("\n" + core.getTime() + " Done!");

    # if count_tops:
    #     print("\n----Topology counts----");
    #     tops_dict = {};
    #     for x in range(len(tops)):
    #         tops_dict[top_trees[x]] = top_counts[x];
    #     for item in sorted(list(tops_dict.items()), key=lambda x: x[1], reverse=True):
    #         print(item[0], item[1]);
    #     print(len(tops_dict), "total topologies found");

    # print("\n----Concordance factor nodes----");
    # stree = tp.addBranchLength(stree, sinfo);
    # for node in node_counts:
    #     cf = round(node_counts[node]/total_trees,2)
    #     print(node, cf);
    #     if sinfo[node][3] not in  ["NA",'']:
    #         stree = stree.replace(node, node + "_" + str(cf));
    #     else:
    #         stree = stree.replace(node, node+ "_" + str(cf));

    # if stree[-1] == "_":
    #     stree = stree[:-1];

    # print("\n----Concordance factor tree----");
    # print(stree);
    # print()

    # print("-----");
    # print(str(num_lines) + " total lines in gene tree file.");
    # if tre_skip != []:
    #     print("The following " + str(len(tre_skip)) + " lines couldn't be read as trees and were skipped: " + ",".join(tre_skip));
    # if sc_skip != []:
    #     print("The following " + str(len(sc_skip)) + " lines were skipped because they didn't have the same number of nodes as the species tree: " + ",".join(sc_skip));
    # print(str(total_trees) + " trees read.");
    # print("=======================================================================");





def pare(globs, cur_tree, iteration, debug=False):
    #debug=True;
    if debug:
        print();
        print("CUR TREE: ", cur_tree);

    cur_tree_dict, cur_labeled_tree, cur_root = TREE.treeParse(cur_tree);
    if debug:
        print("CUR TREE DICT: ", cur_tree_dict)
    cur_tree_dict, bls, na_bl_spec = TREE.adjustTreeDict(cur_tree_dict, cur_root);
    if debug:
        print("CUR TREE DICT: ", cur_tree_dict)
    # Parse the tree string and do some extra manipulation of the tree dict

    cur_bl_tree = TREE.addBranchLength(cur_labeled_tree, cur_tree_dict);
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

            cur_desc = TREE.getDesc(n, cur_tree_dict);
            # Get the two descendant nodes of the node to pare

            clade1 = TREE.getClade(cur_desc[0], cur_tree_dict);
            clade2 = TREE.getClade(cur_desc[1], cur_tree_dict);
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

                        cur_subtree = TREE.getSubtree(cur_desc[d], cur_bl_tree);
                        cur_subtree = re.sub("<[\d]+>_", "", cur_subtree);
                        if debug:
                            print("CUR SUBTREE: ", cur_subtree);
                        # Get the subtree descending from the current descendant node and remove
                        # the previous labels

                        subtree_dict, labeled_subtree, subtree_root = TREE.treeParse(cur_subtree);
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

        anc_desc = TREE.getDesc(anc, cur_tree_dict);
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

        cur_subtree = TREE.getSubtree(anc, cur_bl_tree);
        cur_subtree += anc + "_" + str(anc_supp) + ":" + str(anc_bl);
        if debug:
            print("SUBTREE AT ANCESTOR: ", cur_subtree);
            print("--------");
        # Get the subtree of the ancestor and ad the ancestral node to it

        sis_subtree = TREE.getSubtree(sis, cur_bl_tree);
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

        cur_tree_dict, cur_labeled_tree, cur_root = TREE.treeParse(new_tree);
        cur_tree_dict, bls, na = TREE.adjustTreeDict(cur_tree_dict, cur_root);
        cur_bl_tree = TREE.addBranchLength(cur_labeled_tree, cur_tree_dict);
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