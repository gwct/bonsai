#############################################################################
# Functions to read and handle phylogenetic trees
# Gregg Thomas
#############################################################################

import sys
import re
import copy
from itertools import chain

#############################################################################
class Tree:
# The treeParse function takes as input a rooted phylogenetic tree with branch lengths and returns the tree with node labels and a
# dictionary with usable info about the tree in the following format:
# node:[branch length (if present), ancestral node, node type, node label (if present)]

    def __init__(self, tree_string, debug=False):

        tree = tree_string.strip();
        if tree[-1] != ";":
            tree += ";";
        # Some string handling

        self.orig_tree_str = tree;
        self.topo_str = "";
        self.labeled_topo_str = "";
        self.tree_str = "";
        # Tree string types

        self.type = {};
        self.bl = {};
        self.has_bl = "";
        self.label = {};
        self.has_label = "";
        self.anc = {};
        self.desc = {};
        self.sis = {};
        self.clade = {};
        self.split = {};
        self.quartet = {};
        self.subtree = {};
        # Node attributes

        self.nodes = [];
        self.num_nodes = 0;

        self.tips = [];
        self.num_tips = 0;

        self.internals = [];
        self.num_internals = 0;

        self.rooted = "NA";
        self.root = "";
        # Node lists

        ## Class attributes
        #####

        self.topo_str = remBranchLength(tree);  

        ## Remove the branch lengths and node labels from the input tree string
        #####        

        for tip_label in self.topo_str.replace("(","").replace(")","").replace(";","").split(","):
            self.nodes.append(tip_label);
            self.tips.append(tip_label);
            self.type[tip_label] = 'tip';
        
        ## Retrieval of the tip labels
        #####

        labeled_tree = "";
        z = 0;
        numnodes = 1;
        while z < (len(tree)-1):
            labeled_tree += tree[z];
            if tree[z] == ")":
                node_label = "<" + str(numnodes) + ">";
                labeled_tree += node_label;
                self.nodes.append(node_label);
                self.internals.append(node_label);
                self.type[node_label] = 'internal';
                numnodes += 1;
            z += 1;

        self.internals = self.internals[:-1];
        self.root = node_label;
        # Get the root node as the last node read

        ## Generate internal node labels and add them to the original tree string
        #####

        self.num_tips = len(self.tips);
        self.num_internals = len(self.internals);
        self.num_nodes = len(self.nodes);

        self.rooted = self.checkRooted();
        
        if debug:
            print();
            print("TREE:", self.orig);
            print("TOPOLOGY:", self.topo_str);
            print("NODES:", self.nodes);
            print("ROOTED:", self.rooted);
            print("ROOTNODE:", self.root);
        ## Node counting and root checking
        #####

        self.bl = { node : "NA" for node in self.nodes };
        self.label = { node : "" for node in self.internals + [self.root] };
        self.anc = { node : "NA" for node in self.nodes };

        ## With the node lists set, initialize the other node attributes
        #####

        z = 0;
        numnodes = 1;
        while z < (len(self.topo_str)-1):
            self.labeled_topo_str += self.topo_str[z];
            if self.topo_str[z] == ")":
                node_label = "<" + str(numnodes) + ">";
                self.labeled_topo_str += node_label;
                numnodes += 1;
            z += 1;

        if debug:
            print("LABELED TOPOLOGY:", self.labeled_topo_str);
            print("----------");

        ## Add the generated internal node labels onto the topology string
        #####

        for node in self.nodes:
        # One loop through the nodes to retrieve all other info
            if debug:
                print("NODE:", node);

            if node in self.tips:
                if node + ":" in tree:
                    cur_bl = re.findall(node + ":[\d.Ee-]+", labeled_tree);
                    cur_bl = cur_bl[0].replace(node + ":", "");
                    if debug:
                        print("FOUND BL:", cur_bl);
                    self.bl[node] = cur_bl;
                ## If there are branch lengths, parse with regex and add to bl dict
            ## Parse tips
             
            elif node in self.internals:
                if node + "(" in labeled_tree or node + "," in labeled_tree or node + ")" in labeled_tree:
                    if debug:
                        print("NO BL OR LABEL");
                ## If there's no labels or branch lengths, keep the NAs

                elif node + ":" in labeled_tree:
                    cur_bl = re.findall(node + ":[\d.Ee-]+", labeled_tree);
                    cur_bl = cur_bl[0].replace(node + ":", "");
                    if debug:
                        print("FOUND BL:", cur_bl);
                    self.bl[node] = cur_bl;
                ## If there's only branch lengths, add them to the bl dict                               

                else:
                    cur_label_str = re.findall(node + "[\d\w<>_*+.Ee/-]+:[\d.Ee-]+", labeled_tree);
                    # If this pattern is found there is both a branch length and a label for the node

                    if cur_label_str:
                        cur_label = cur_label_str[0].replace(node, "");
                        cur_label = cur_label[:cur_label.index(":")];
                        cur_bl = cur_label_str[0].replace(node, "").replace(cur_label, "").replace(":", "");
                        if debug:
                            print("FOUND BL AND LABEL:", cur_bl, cur_label);
                        self.label[node] = cur_label;
                        self.bl[node] = cur_bl;
                    ## Parse if there is both bl and label
                        
                    else:
                        cur_label = re.findall(node + "[\w*+.<> -]+", labeled_tree);
                        cur_label = cur_label[0].replace(node, "");
                        if debug:
                            print("FOUND LABEL:", cur_label);
                        self.label[node] = cur_label;
                    ## Parse if there is only a label
                ## Check if there are branch lenghts and labels or just labels
            ## Parse internal nodes

            elif node == self.root:
                possible_label = labeled_tree[labeled_tree.index(node)+len(node):];
                if possible_label:
                    self.label[node] = possible_label;
                continue;
            ## Parse the root and continue since there is no ancestor to parse in the block below

            #####

            anc_tree = labeled_tree[labeled_tree.index(node):][1:];
            # Get the tree string starting at the current node label
            # Ancestral labels are always to the right of the node label in the text of the tree, 
            # so we start our scan from the node label

            if debug:
                print("NODE:", node);
                print("ANC_TREE:", anc_tree);
                
            cpar_count = 0;
            cpar_need = 1;

            for i in range(len(anc_tree)):
                if anc_tree[i] == "(":
                    cpar_need = cpar_need + 1;
                if anc_tree[i] == ")" and cpar_need != cpar_count:
                    cpar_count = cpar_count + 1;
                if anc_tree[i] == ")" and cpar_need == cpar_count:
                    anc_tree = anc_tree[i+1:];
                    self.anc[node] = anc_tree[:anc_tree.index(">")+1];
                    break;
                # When the parentheses counts match, the ancestor will start at the next position and
                # end at the >
            # We find the ancestral label by finding the ) which matches the nesting of the number of ('s found

            if debug:
                print("FOUND ANC:", self.anc[node]);
                print("---");

            #####
        ## End node loop

        if debug:
            for node in self.anc:
                print(node, self.anc[node]);

        #####

        if all(self.bl[n] == "NA" for n in self.bl):
            self.has_bl = False;
        else:
            self.has_bl = True;
            self.bl = { n : self.bl[n] if self.bl[n] != "NA" else "0.0" for n in self.bl };

        if all(self.label[n] == "" for n in self.label):
            self.has_label = False;
        else:
            self.has_label = True;
            self.label = { n : self.label[n] if self.label[n] != "" else "NA" for n in self.label };

        for node in self.nodes:
            self.desc[node] = self.getDesc(node);
            self.sis[node] = self.getSister(node);
            self.clade[node] = set(self.getClade(node));
            self.split[node] = self.getSplit(node);

        # self.showAnc();
        # print("--------------");
        # self.showClade();
        # print("--------------");
        # self.showSplit();
        # print("--------------");



        self.subtree = self.genSubtrees();
        # A few more useful node attributes to store

        if not self.rooted:
            self.internals += [self.root];
        # For unrooted trees, add the "root" node as an internal node
        # In this context, the "root" node is simply the last node read

        self.tree_str = self.subtree[self.root];
        # Full tree with node names, labels, and branch lengths

        # print(self.tree_str);

        for node in self.internals:
            self.quartet[node] = self.getQuartet(node);

    ##########

    def checkRooted(self):
        if (self.num_internals + 1) != (self.num_tips - 1):
            return False;
        elif (self.num_internals + 1) == (self.num_tips - 1):
            return True;
        else:
            return -1;

    ##########

    def getDesc(self, node):
        # This function takes a node in the current tree object
        # and returns a list of the two direct descendant nodes of it.

        if node in self.tips:
            return [node];
        else:
            return [ n for n in self.nodes if self.anc[n] == node ];

    ##########

    def getSister(self, node):
        # This function takes a node in the current tree object
        # and returns the other direct descendant of its ancestral node

        if node == self.root:
            return "NA";
        anc_desc = self.getDesc(self.anc[node]);
        return [ n for n in anc_desc if n != node ][0];

    ##########

    def getClade(self, node, full=False):
    # This function takes a node in the current tree object
    # and finds all tip labels that are descendants of it.
    # This is done by getting the direct descendants of the node with getDesc and then
    # recursively calling itself on those descendants.

        clade = [];
        desc = self.getDesc(node);
        for d in desc:
            if d not in self.tips:
                clade += self.getClade(d, full);
                if full:
                    clade.append(d);
                # If full is true, the function will also return all internal nodes
                # descending from the given node
            else:
                clade.append(d);

        return clade;

    ##########

    def getSplit(self, node):
    # Returns the tips from a node that do not descend from it to define a split (with clade)

        if node == self.root:
            return "NA";
        return set(self.tips) - self.clade[node];

    ##########

    def getQuartet(self, node):
    # Returns the tips from a node that do not descend from it to define a split (with clade)

        # print(node);
        if self.rooted:
            if node == self.root or self.type[node] == "tip":
                return "NA";
            # elif self.anc[node] == self.root and self.type[self.sis[node]] == "tip":
            #     return "NA";
            else:
                if self.anc[node] == self.root:
                    d1 = self.desc[node][0];
                    d2 = self.desc[node][1];

                    q1 = set(self.clade[d1]);
                    q2 = set(self.clade[d2]);

                    if self.type[self.sis[node]] == "tip":
                        q3 = set(self.sis[node]);
                        q4 = set(self.sis[node]);

                    else:
                        sis_desc = self.desc[self.sis[node]];
                        q3 = set(self.clade[sis_desc[0]]);
                        q4 = set(self.clade[sis_desc[1]]);              

                else:
                    d1 = self.desc[node][0];
                    d2 = self.desc[node][1];
                    sis = self.sis[node];

                    q1 = set(self.clade[d1]);
                    q2 = set(self.clade[d2]);
                    q3 = set(self.clade[sis]);
                    q4 = set(self.tips) - q1 - q2 - q3;

        else:
            if node == self.root or self.type[node] == "tip":
                return "NA";
            else:
                d1 = self.desc[node][0];
                d2 = self.desc[node][1];
                q1 = set(self.clade[d1]);
                q2 = set(self.clade[d2]);

                if len(self.desc[self.anc[node]]) > 2:
                    other = [ n for n in self.desc[self.anc[node]] if n != node ];
                    q3 = set(self.clade[other[0]]);
                    q4 = set(self.clade[other[1]]);    
                else:
                    q3 = set(self.clade[self.sis[node]]);
                    q4 = set(self.tips) - q1 - q2 - q3;                


        return { 'd1' : q1, "d2" : q2, "s" : q3, "q4" : q4 };

    ##########

    def genSubtrees(self):
        subtrees = {};
        
        for node in self.tips:
            if self.has_bl:
                subtrees[node] = node + ":" + str(self.bl[node]);
            else:
                subtrees[node] = node;

        for node in self.internals + [self.root]:
            subtrees[node] = "(";
            for d in self.desc[node]:
                subtrees[node] += subtrees[d] + ",";
            subtrees[node] = subtrees[node][:-1] + ")" + node;

            if self.has_label:
                subtrees[node] += str(self.label[node]);
            if self.has_bl and node != self.root:
                subtrees[node] += ":" + str(self.bl[node]);

        return subtrees;

    ##########

    def findClades(self, tip_set):
    # A function that takes a set of tips and finds all the monophyletic
    # clades within that set, and returns the LCA of each clade found

        clade_nodes = [];
        for node in self.internals:
            if self.clade[node] <= tip_set:
                clade_nodes.append(node);
        # In the input tree, find all nodes that have clades that are
        # a subset of the given tip set

        nodes_to_rm = [];
        for node in clade_nodes:
            for node_check in clade_nodes:
                if node == node_check:
                    continue;

                if node in self.getClade(node_check, full=True):
                #if node in self.desc[node_check]:
                    nodes_to_rm.append(node);
        clade_set = [ n for n in clade_nodes if n not in nodes_to_rm ];
        # We only want the deepest node from each possible clade, so remove
        # nodes from those found that are descendants of another node.

        nodes_to_rm = [];
        for node in tip_set:
            if any(node in self.getClade(n, full=True) for n in clade_set):
                nodes_to_rm.append(node);
        clade_set += [ n for n in tip_set if n not in nodes_to_rm ];
        # From the original input set, remove any nodes that are now found
        # in a clade

        return set(clade_set);

    ##########

    def findSplits(self, tip_set):
    # A function that takes a set of tips and finds whether any branches
    # are defined on either side by them

        split_nodes = [];
        for node in self.nodes:
            if self.clade[node] == tip_set or self.split[node] == tip_set:
                split_nodes.append(node);

        return set(split_nodes);

    ##########

    def LCA(self, node_list):
    # Given a list of nodes, this function finds the least common ancestor of them,
    # and tells whether the nodes provided form a monophyletic clade.

        ancs = {};
        for node in node_list:
            ancs[node] = [node];
        # For each node in the input list, we make a list of the path from that node 
        # to the root of the tree, including that node

        for node in node_list:
            if node == self.root:
                continue;

            cur_anc = self.anc[node];
            ancs[node].append(cur_anc);

            while cur_anc != self.root:
                cur_anc = self.anc[node];
                ancs[node].append(cur_anc);
        # For each node, add every node between it and the root to its ancs list


        intersect_anc = set.intersection(*list(map(set, list(ancs.values()))));
        # Gets the intersect of all lists of paths to the root, unordered

        lcp = [t for t in list(ancs.values())[0] if t in intersect_anc];
        # Orders the nodes in the intersect of all paths based on their order in
        # the path of an arbitrary node (since it should be identical for all of them)

        return lcp[0];
        # Returns the first node in the common path as the LCA  

    ##########

    def Monophyletic(self, node_list):
        monophyletic = False;
        if set(self.getClade(self.LCA(node_list))) == set(node_list):
            monophyletic = True;
        return monophyletic;

    ##########

    def addBranchLength(self):
    # Re-writes the branch lengths onto a tree object topology

        tree = self.labeled_topo_str;

        for node in self.nodes:
            new_node_str = node;
            if self.has_label and node in self.internals:
                new_node_str += str(self.label[node]);

            if self.has_bl and node != self.root:
                new_node_str += ":" + str(self.bl[node]);

            tree = tree.replace(node, new_node_str);

        return tree;

    ##########

    def addLabel(self, label_dict):
    # Given a dictionary with { node : label } format, adds those labels and the branch
    # lengths onto the given tree's topology
    # Returns: tree string

        new_tree_str = self.labeled_topo_str;
        # Extract labeled tree topology to add labels to

        for node in self.nodes:
            new_label = "";
            if self.type[node] == "tip":
                new_label = node;
            # If the node is a tip, always add the node as a label

            if node in label_dict:
                new_label += str(label_dict[node]);
            # If the node is in the given dictionary, add the new label

            if self.has_bl:
                new_label += ":" + str(self.bl[node]);
            # If the tree has branch lengths, add the bl

            if new_label:
                new_tree_str = new_tree_str.replace(node, new_label);
            # If a new label was created, replace the old node label in the
            # tree with it

        return new_tree_str;

    ##########

    def Prune(self, node_list, debug=False):
    # Given a node in the tree, this prunes that node out of the tree

        working_tree = self;
        
        for node in node_list:

            forward_map, prune_map = mapNodes(working_tree, self);

            if self.type[node] == 'tip':
                prune_node = node;
            else:
                prune_node = prune_map[node];
                #prune_node = [ k for k,v in prune_map.items() if v == node ][0];

            if debug:
                #print(working_tree.label);
                #print(working_tree.bl);
                print("TREE TO PRUNE: ", working_tree.tree_str);
                #print(working_tree.tree);
                #print(working_tree.topo);
                #print(working_tree.anc);
                #print(working_tree.tree);
                print("ORIG NODE:", node);
                print("NODE TO PRUNE:", prune_node);
                print("--------");

            ## Debug info
            #####

            anc = working_tree.anc[prune_node];
            if debug:
                print("ANCESTOR OF NODE TO PRUNE: ", anc);
                print("--------");
            
            ## Lookup the ancestor of the current node 
            #####

            sis = working_tree.sis[prune_node];
            if debug:
                print("SISTER OF NODE TO PRUNE: ", sis);
                print("--------");

            ## Look up the sister of the current node
            #####

            #cur_subtree = getSubtree(anc, working_tree.bltree);
            #cur_subtree += anc;
            cur_subtree = working_tree.subtree[anc];
            # if working_tree.has_label:
            #     cur_subtree += "_" + str(working_tree.label[anc]);
            # if working_tree.has_bl:
            #     cur_subtree += ":" + str(working_tree.bl[anc]);           

            if debug:
                print("SUBTREE AT ANCESTOR: ", cur_subtree);
                print("--------");
            
            ## Get the subtree of the ancestor and add the ancestral node to it
            #####

            # if working_tree.has_bl and sis in working_tree.tips:
            #     sis_subtree = sis + new_bl;

            # else:
                # sis_subtree = getSubtree(sis, working_tree.bltree);
            sis_subtree = working_tree.subtree[sis];
            #new_subtree = sis_subtree + sis;
            # if working_tree.has_label:
            #     new_subtree += "_" + str(working_tree.label[sis]);
            # new_subtree += new_bl;
            if debug:
                print("SUBTREE AT SISTER: ", sis_subtree);
                print("--------");
            # Get the subtree of the sister

            if working_tree.has_bl and anc != working_tree.root:
                new_bl = ":" + str(float(working_tree.bl[anc]) + float(working_tree.bl[sis]));
                sis_subtree = sis_subtree[:sis_subtree.rfind(":")] + new_bl;
            # After removing the current species, the new branch length of the combined anc and sis
            # branch will be their sum

            if debug:
                print("NEW SUBTREE: ", sis_subtree);
                print("--------");
            # Make the new subtree to replace the current ancestral one. If the sister is
            # a tip, then this is just the sister label with the new branch length.
            # If the sister is a clade, then this is the sister subtree with the sister added
            # on with the new branch length

            #####

            if anc == working_tree.root:
                new_tree = sis_subtree;
            else:
                new_tree = working_tree.tree_str.replace(cur_subtree, sis_subtree);
            new_tree = re.sub("<[\d]+>", "", new_tree);

            # Replace the old subtree with the new one and remove the internal node labels from treeParse 

            #####

            working_tree = Tree(new_tree);
            if debug:
                print("NEW TREE: " + working_tree.tree_str);
                print("--------");
            #####

            #cur_tree_dict, cur_labeled_tree, cur_root = TREE.treeParse(new_tree);
            #cur_tree_dict, bls, na = TREE.adjustTreeDict(cur_tree_dict, cur_root);
            #cur_bl_tree = TREE.addBranchLength(cur_labeled_tree, cur_tree_dict);
            # Parse the new tree and add the branch lengths back on

        return working_tree;

    ##########

    def Prune2(self, node_list, debug=False):
        subtrees = {};
        # A dict of subtrees for each node in the tree

        pruned_desc_dict = copy.deepcopy(self.desc);
        # A local copy of the descendants dict to adjust for pruning

        node_list = list(self.findClades(set(node_list)));
        # Make sure all the nodes to prune represent monophyletic groups (else this algorithm will error)

        for node in self.tips:
            if self.has_bl:
                subtrees[node] = node + ":" + str(self.bl[node]);
            else:
                subtrees[node] = node;
        ## Add tips to subtree dict
        #####

        for node in self.internals:
            # print("NODE:", node);
            desc_to_prune = [d for d in pruned_desc_dict[node] if d in node_list];
            desc_to_keep = [d for d in pruned_desc_dict[node] if d not in node_list];
            # print(desc_to_prune);
            # print(desc_to_keep);
            # For the current node, get which descendants to prune and which to keep

            if len(desc_to_prune) in [0, 2]:
                subtrees[node] = "(";
                for d in pruned_desc_dict[node]:
                    subtrees[node] += subtrees[d] + ",";
                subtrees[node] = subtrees[node][:-1] + ")" + node;

                if self.has_label:
                    subtrees[node] += str(self.label[node]);
                if self.has_bl and node != self.root:
                    subtrees[node] += ":" + str(self.bl[node]);
            # If neight or both descendants are pruned, build the subtree for the current node

            elif desc_to_prune:
            # If one descendant is pruned, make adjustments
                if not self.rooted and node == self.root:
                    continue;

                anc = self.anc[node];
                # print("ANC:", anc);
                # print("ANC DESC:", self.desc[anc]);
                desc_to_keep = desc_to_keep[0];
                # Unpack the node to keep and its ancestor

                if self.has_bl:
                    new_bl = ":" + str(float(self.bl[desc_to_keep]) + float(self.bl[node]));
                    subtrees[desc_to_keep] = subtrees[desc_to_keep][:subtrees[desc_to_keep].rfind(":")] + new_bl;
                # If the tree has branch lengths, we need to add the bl of the current node with the bl
                # of its ancestor and adjust the subtree to match
                            
                anc_desc_ind = pruned_desc_dict[anc].index(node);
                pruned_desc_dict[anc][anc_desc_ind] = desc_to_keep;
                # In the ancestral node, replace this node as a descendant with the node to keep                
            # print("-------")
        ## Construct subtrees for each internal node, adjusting when a node needs to be pruned
        #####

        if any(d in node_list for d in self.desc[self.root]):
        # If either of the nodes descending from the last (root) node are to be pruned, then the final
        # tree is just the subtree from the other node.
            if self.rooted:
                desc_to_keep = [d for d in pruned_desc_dict[self.root] if d not in node_list][0];
                pruned_tree = subtrees[desc_to_keep] + ";";
                # Find the node that is to be kept and lookup its subtree
                
                if self.has_bl:
                    pruned_tree = pruned_tree[:pruned_tree.rfind(":")];
                # Remove the trailing branch length which led to the root

            else:
                desc_to_keep = [d for d in pruned_desc_dict[self.root] if d not in node_list];
                pruned_tree = "(";
                for d in desc_to_keep:
                    pruned_tree += subtrees[d] + ",";
                pruned_tree = pruned_tree[:-1] + ");";


        else:
        # If both descendants from the last (root) node are to be kept, combine them
        # to the final tree
            if self.rooted:
                l = subtrees[pruned_desc_dict[self.root][0]];
                r = subtrees[pruned_desc_dict[self.root][1]];
                # Lookup both descendant subtrees

                pruned_tree = "(" + l + "," + r + ");";

            else:
                l = subtrees[pruned_desc_dict[self.root][0]];
                m = subtrees[pruned_desc_dict[self.root][1]];
                r = subtrees[pruned_desc_dict[self.root][2]];

                pruned_tree = "(" + l + "," + m + "," + r + ");";
            # if self.has_bl:
            #     pruned_tree = pruned_tree[:pruned_tree.rfind(":")] + ")";
            # Combine the subtrees and add the last branch length if necessary
        ## Make the final tree by combining the subtrees descending from the last (root) node
        #####

        pruned_tree = re.sub("<[\d]+>", "", pruned_tree);
        # Remove the added node labels

        return pruned_tree;

    ##########

    def Root(self, node_list):

        if len(node_list) > 1:
            clade_list = list(self.findClades(set(node_list)));

            if len(clade_list) == 2:
                if self.sis[clade_list[0]] != clade_list[1]:
                    return "cannot root";



            if len(clade_list) > 2:
                return "cannot root";     

        ## find clades... if 1, ok... if 2, check if split... if 3, or more can't root

    ##########

    def showType(self):
        for node in self.nodes:
            outline = node + "\t\t" + self.type[node];
            print(outline);

    ##########

    def showDesc(self):
        for node in self.nodes:
            outline = node + "\t\t" + ",".join(self.desc[node]);
            print(outline);

    ##########

    def showAnc(self):
        for node in self.nodes:
            outline = node + "\t\t" + self.anc[node];
            print(outline);

    ##########

    def showSis(self):
        for node in self.nodes:
            outline = node + "\t\t" + self.sis[node];
            print(outline);

    ##########

    def showClade(self):
        for node in self.nodes:
            outline = node + "\t\t" + str(self.clade[node]);
            print(outline);

    ##########

    def showSplit(self, n=False):
        for node in self.nodes:
            if not n or n == node:
                outline = node + "\t\t" + str(self.clade[node]) + "\t" + str(self.split[node]);
                print(outline);

    ##########

    def showQuartet(self, n=False):
        for node in self.internals:
            if node == self.root:
                continue;
            if not n or n == node:
                outline = node + "\t\t";
                for key in ['d1', 'd2', 's', 'q4']:
                    outline += "\t" + str(self.quartet[node][key]);
                print(outline);


#############################################################################

def remBranchLength(tree_str):
# Removes branch lengths from a tree.

    tree_str = re.sub('[)][\d\w<>/.eE_:-]+', ')', tree_str);
    tree_str = re.sub(':[\d.eE-]+', '', tree_str);

    return tree_str;

#############################################################################

def getSubtree(node, tree_str):
# Gets the subtree string at a given node from a labeled tree string from treeParse

    subtree = "";
    # Initialize the subtree string

    partree = tree_str[:tree_str.index(node)][::-1]
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
# Gets some info about branch lengths from the tree dict and retrieves support
# label

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

def mapNodes(query_tree, ref_tree):
# Species must have same tip labels, except for missing/pruned species
# Smaller tree must be query tree

    node_map = {};
    for query_node in query_tree.internals:
        max_node = "";
        max_node_set = {};

        # print(query_node);

        for ref_node in ref_tree.internals:
            check_node = False;

            # print("\t", ref_node);

            if not query_tree.rooted and query_node == query_tree.root:
                if query_tree.clade[query_node] <= ref_tree.clade[ref_node] or query_tree.clade[query_node] <= ref_tree.split[ref_node]:
                    check_node = True;

            elif query_tree.clade[query_node] <= ref_tree.clade[ref_node] and query_tree.split[query_node] <= ref_tree.split[ref_node]:
                check_node = True;

            # elif query_tree.split[query_node] <= ref_tree.clade[ref_node] and query_tree.clade[query_node] <= ref_tree.split[ref_node]:
            #     check_node = True;

            # if pruned_node == "<3>" and orig_node == "<9>":
            #     print(orig_node, check_node);
            #     print(cur_st.clade[pruned_node]);
            #     print(cur_st.split[pruned_node]);
            #     print(st.clade[orig_node]);
            #     print(st.split[orig_node]);                     

            if check_node:
                node_set = ref_tree.clade[ref_node] | ref_tree.split[ref_node];
                if len(node_set) > len(max_node_set):
                    max_node_set = node_set;
                    max_node = ref_node;

        node_map[query_node] = max_node;

    rev_map = { v : k for k,v in node_map.items() if v != "" };

    return node_map, rev_map;    

#############################################################################