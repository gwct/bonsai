#############################################################################
# Handles the reading of the input trees
# Gregg Thomas
#############################################################################

import lib.core as CORE
import lib.tree as TREE

#############################################################################

def readST(globs):

    if globs['st-input-type'] == "file":
        globs['orig-st-str'] = open(globs['st-input'], "r").read().strip();
    else:
        globs['orig-st-str'] = globs['st-input']
    # If the input type is a file, read the file here, otherwise take the string as the tree input

    try:
        globs['st-dict'], globs['labeled-st-str'], globs['root'] = TREE.treeParse(globs['orig-st-str']);
        #print(globs['labeled-tree-str']);
    except:
        print("\n\n");
        CORE.errorOut("1", "Error reading tree! Make sure it is formatted as a rooted, Newick tree.", globs);
    # Try and parse the tree with treeParse() and erorr out if fail

    globs['tips'] = [ n for n in globs['st-dict'] if globs['st-dict'][n][2] == 'tip' ];
    globs['num-tips'] = len(globs['tips']);
    globs['num-internals'] = len([ n for n in globs['st-dict'] if globs['st-dict'][n][2] == 'internal' ]);
    # Count the nodes

    return globs;

#############################################################################

def readGT(globs):

    line_num, tree_num = 0, 0;
    for line in open(globs['gt-input']):
        line = line.strip();
        line_num += 1;

        if not line:
            globs['gt-input-empty'] += 1;
            continue;

        try:
            gt_dict, gt_str, gt_root = TREE.treeParse(line);
            tree_num += 1;
            globs['gt-init'][tree_num] = { 'info' : gt_dict, 'str' : gt_str, 'root' : gt_root, 'orig-str' : line };
        except:
           globs['warnings'] += 1;
           globs['gt-input-skipped'] += 1;
           CORE.printWrite(globs['logfilename'], globs['log-v'], "# WARNING: Could not read line " + str(line_num) + " of gene tree file (-g) as a tree! Skipping.");
           # Warning statement if a line can't be read as a tree
           
    return globs;

#############################################################################

def readExempt(globs):

    for line in open(globs['exempt-file']):
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

                if not all(s in globs['tips'] for s in specs):
                    globs['warnings'] += 1;
                    CORE.printWrite(globs['logfilename'], globs['log-v'], "# WARNING: Label in exempt file not found in tree. Skipping line: " + line);
                # Throw a warning if both of the tips aren't found in the tree and skip

                else:
                    exempt_node, mono_flag = TREE.LCA(specs, globs['st-dict']);
                    # Get the least common ancestor of the given tips

                    CORE.printWrite(globs['logfilename'], globs['log-v'], "# INFO: " + line + " -> " + exempt_node);
                    # Info statement to show branch assignment

                    globs['exempt-branches'].append(exempt_node);
                    # Add the internal branch to the list of exempt branches

                    globs['exempt-clades'].append(set(TREE.getClade(exempt_node, globs['st-dict'])));
                    # Get and add the full clade descending from the internal branch to the list of exempt clades
                # If both tips are found, get the full clade
            ## End tip block

            else:
            # If there is no space, assume it is defining an internal branch by label

                if line not in globs['st-dict']:
                    globs['warnings'] += 1;
                    CORE.printWrite(globs['logfilename'], globs['log-v'], "# WARNING: Label in exempt file not found in tree. Skipping line: " + line);
                # Print a warning if the branch label isn't found in the tree

                else:
                    globs['exempt-branches'].append(line);
                    globs['exempt-clades'].append(set(TREE.getClade(line, globs['st-dict'])));
                # Otherwise, add the branch and its clade to their global lists
            ## End label block
        ## End line/branch loop   

    return globs;

#############################################################################