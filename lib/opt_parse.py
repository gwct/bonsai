#############################################################################
# Parsing and printing the options and meta-info for PhyloAcc.
# Much of the error checking is done here as well.
#############################################################################

import sys
import os
import math
import argparse
import multiprocessing as mp
import lib.core as CORE

#############################################################################

def optParse(globs):
# This function handles the command line options and prepares the output directory and files.
# Defaults are set in params.py
    try:
        import psutil
        globs['psutil'] = True;
    except:
        globs['psutil'] = False;
    # Check if psutil is installed for memory usage stats.

    ###############

    parser = argparse.ArgumentParser(description="Heuristic tree paring");

    parser.add_argument("-st", dest="st_input", help="A file or a string containing a single rooted, Newick formatted tree. REQUIRED", default=False);
    parser.add_argument("-gt", dest="gt_input", help="A file containing many rooted, Newick formatted gene trees, one per line. REQUIRED", default=False);
    # Input

    parser.add_argument("-o", dest="out_dest", help="Desired output directory. This will be created for you if it doesn't exist. Default: pare-[date]-[time]", default=False);
    # Output
    
    parser.add_argument("-b", dest="bl_percentile", help="The lower percentile of branch lengths to consider for paring at each iteration. Must be between 0 and 100. Default: 10.", type=int, default=0);
    parser.add_argument("-g", dest="gcf_threshold", help="The lower threshold of gene concordance factor to consider for paring at each iteration. Must be between 0 and 100. Default: 25.", type=float, default=0);
    parser.add_argument("-i", dest="max_iters", help="The maximum number of paring iterations to perform before stopping. Default: 3.", type=int, default=0);
    parser.add_argument("-s", dest="total_max_spec", help="The total maximum number of species that can be pruned while paring. Default: Unlimited.", type=int, default=0);
    parser.add_argument("-m", dest="branch_max_spec", help="The maximum number of species that can be pruned while paring any given single branch. Default: 10.", type=int, default=0);
    parser.add_argument("-p", dest="prune", help="A file containing branches to prune from the input tree(s), one branch per line defined as internal node labels with --labeltree or 2 tips that descend from that branch (e.g. 'spec1 spec2'). Lines in file starting with '#' are ignored.", default=False);
    parser.add_argument("-e", dest="exempt", help="A file containing branches to be exempt from pruning (note that they can still be PARED), one branch per line defined as internal node labels with --labeltree or 2 tips that descend from that branch (e.g. 'spec1 spec2'). Lines in file starting with '#' are ignored.", default=False);
    parser.add_argument("-d", dest="aln_dir", help="A directory of corresponding alignments to calculate site concordance factors on the species tree (optional).", default=False);
    parser.add_argument("-scf", dest="scf_samples", help="The number of quartets to sample around each branch for sCF calculations.", type=int, default=0);
    parser.add_argument("-n", dest="num_procs", help="The number of processes that this script should use. Default: 1.", type=int, default=1);
    # User params

    parser.add_argument("--cf", dest="cf_flag", help="Stop after calculating concordance factors.", action="store_true");
    parser.add_argument("--stats", dest="aln_stats_flag", help="When a directory of alignments is given with -scf, set this flag to write out some stats for the alignments.", action="store_true", default=False);
    #parser.add_argument("--count", dest="topo_count_flag", help="Exit after counting topologies in the gene tree file.", action="store_true");
    parser.add_argument("--labels", dest="labels_flag", help="Set this option to disable concordance factor calculations and use the labels on the input tree for pruning. Note that if these labels are NOT concordance factors, you may also need to adjust -g.", action="store_true", default=False);
    parser.add_argument("--labeltree", dest="labeltree", help="Just read the tree from the input, label the internal nodes, print, and exit.", action="store_true");
    parser.add_argument("--overwrite", dest="ow_flag", help="Set this to overwrite existing files.", action="store_true", default=False);
    parser.add_argument("--appendlog", dest="append_log_flag", help="Set this to keep the old log file even if --overwrite is specified. New log information will instead be appended to the previous log file.", action="store_true", default=False);
    # User options

    parser.add_argument("--info", dest="info_flag", help="Print some meta information about the program and exit. No other options required.", action="store_true", default=False);
    #parser.add_argument("--dryrun", dest="dryrun", help="With all options provided, set this to run through the whole pseudo-it pipeline without executing external commands.", action="store_true", default=False);
    parser.add_argument("--version", dest="version_flag", help="Simply print the version and exit. Can also be called as '-version', '-v', or '--v'", action="store_true", default=False);
    parser.add_argument("--quiet", dest="quiet_flag", help="Set this flag to prevent degenotate from reporting detailed information about each step.", action="store_true", default=False);
    # Run options

    parser.add_argument("--norun", dest="norun", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--debug", dest="debug_opt", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--debugtree", dest="debug_tree", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--nolog", dest="nolog_opt", help=argparse.SUPPRESS, action="store_true", default=False);
    # Performance tests
    args = parser.parse_args();
    
    ## The input options and help messages
    ###############

    globs['call'] = " ".join(sys.argv);
    # Save the program call for later

    ###############

    if args.info_flag:
        globs['info'] = True;
        globs['log-v'] = -1;
        startProg(globs);
        return globs;
    # Parse the --info option and call startProg early if set

    ###############

    if args.norun:
        globs['norun'] = True;
        globs['log-v'] = -1;
    # Check if norun is set

    globs['overwrite'] = args.ow_flag;
    # Check if overwrite is set
    
    if args.debug_tree:
        globs['debug-tree'] = True;
        globs['log-v'] = -1;
    # Check if tree debugging is turned on -- will just read the input species tree for testing in the tree library

    ###############

    if args.labeltree:
        globs['label-tree'] = True;
    # A flag to just print the labeled tree out and exit

    ###############

    if args.cf_flag:
        globs['cf-only'] = True;
    # A flag to exit after calculating concordance factors
    
    ###############

    if not args.st_input:
        CORE.errorOut("OP1", "A tree must be provided with -st.", globs);
    # Check that some input is provided to the tree option

    if os.path.isfile(args.st_input):
        globs['st-input-type'] = "file";
    else:
        globs['st-input-type'] = "string";
    # Check whether the tree input is a file or a string

    globs['st-input'] = args.st_input;
    # Set the tree input param

    ## Species tree inputs
    ###############

    if not globs['label-tree']:

        if args.labels_flag:
            globs['use-labels'] = True;
        # Check whether or not to use the labels on the input species tree for pruning

        ###############

        if not args.gt_input and not globs['use-labels'] and not globs['debug-tree'] and not args.prune:
            CORE.errorOut("OP2", "A set of gene trees must be provided with -gt", globs);
        globs['gt-input'] = args.gt_input;
        # Check if some input is provided to the gene tree option

        ###############

        globs['num-procs'] = CORE.isPosInt(args.num_procs, default=1);
        # Parse the number of processes to use

        ###############

        if args.aln_dir:
            if not os.path.isdir(args.aln_dir):
                CORE.errorOut("OP3", "Cannot find input alignment directory for (-d)", globs);
            globs['aln-dir'] = os.path.abspath(args.aln_dir);

            if args.scf_samples:
                if not CORE.isPosInt(args.scf_samples):
                    CORE.errorOut("OP4", "The number of sCF quartets sampled (-scf) must be a positive integer.", globs);
                globs['scf-quartets'] = args.scf_samples;

            globs['scf-pool'] = mp.Pool(processes=globs['num-procs']);
            if not globs['cf-only']:
                globs['scf-pool-prune'] = mp.Pool(processes=globs['num-procs']);
            # Create the pool of processes for sCF calculation here so we copy the memory profile of the parent process
            # before we've read any large data in
        # Read the sCF alignment directory

        else:
            if args.aln_stats_flag:
                CORE.errorOut("OP5", "An alignment directory must be given with -d to calculate alignment statistics (--stats)", globs);
            if args.scf_samples:
                CORE.errorOut("OP6", "An alignment directory must be given with -d to calculate sCF (-scf)", globs);
        # Check if a directory of alignments has been given to calculate alignment stats

        ## Check alignment directory input
        ###############

        # if args.topo_count_flag:
        #     if not globs['gt-input']:
        #         CORE.errorOut("OP5", "A file with gene trees must be given with -gt to count topologies with --count", globs);
        #     globs['count-only'] = True;
        # Option to count topologies in gene tree file and exit

        ###############

        if args.bl_percentile:
            if not CORE.isPosInt(args.bl_percentile, maxval=100):
                CORE.errorOut("OP7", "The lower branch length percentile (-b) must be an integer between 0 and 100.", globs);
            globs['bl-percentile'] = args.bl_percentile;
        # Check and set the branch length percentile option

        ###############

        if args.gcf_threshold:
            if not CORE.isPropFloat(args.gcf_threshold):
                CORE.errorOut("OP8", "The gCF threshold (-g) must be an integer between 0 and 100.", globs);
            globs['gcf-threshold'] = args.gcf_threshold;
        # Check and set the gcf threshold option

        ###############

        if args.max_iters:
            if not CORE.isPosInt(args.max_iters):
                CORE.errorOut("OP9", "The maximum number of paring iterations (-i) must be a positive integer.", globs);
            globs['max-iterations'] = args.max_iters;
        # Check and set the max iteration option

        ###############

        if args.total_max_spec:
            if not CORE.isPosInt(args.total_max_spec):
                CORE.errorOut("OP10", "The total maximum number of species to prune (-s) must be a positive integer.", globs);
            globs['total-max-spec'] = args.total_max_spec;
        # Check and set the total max number of species to remove

        ###############

        if args.branch_max_spec:
            if not CORE.isPosInt(args.branch_max_spec):
                CORE.errorOut("OP11", "The maximum number of species to prune per branch (-b) must be a positive integer.", globs);
            globs['branch-max-spec'] = args.branch_max_spec;
        # Check and set the max number of species to remove for a single branch

        ###############

        if args.prune:
            if not os.path.isfile(args.prune):
                CORE.errorOut("OP12", "Cannot find file with branches to prune (-p).", globs);
            globs['prune-file'] = args.prune;
        # Check the prune file

        ###############

        if args.exempt:
            if not os.path.isfile(args.exempt):
                CORE.errorOut("OP13", "Cannot find file with branches exempt from pruning (-e).", globs);
            globs['exempt-file'] = args.exempt;
        # Check the exempt file

        ###############

        if not args.out_dest:
            globs['outdir'] = "bonsai-out-" + globs['startdatetime'];
        else:
            globs['outdir'] = args.out_dest;
        # Get the output directory

        if not globs['overwrite'] and os.path.exists(globs['outdir']):
            CORE.errorOut("OP14", "Output directory already exists: " + globs['outdir'] + ". Specify new directory name OR set --overwrite to overwrite all files in that directory.", globs);
        if not os.path.isdir(globs['outdir']) and not globs['norun'] and not globs['info'] and not globs['debug-tree']:
            os.makedirs(globs['outdir']);
        # Main output dir
        
        ###############

        if globs['st-input-type'] == "file":
            globs['st-cf-tree-file'] = os.path.join(globs['outdir'], os.path.splitext(os.path.basename(globs['st-input']))[0] + "-cf.tre");
            globs['st-cf-stat-file'] = os.path.join(globs['outdir'], os.path.splitext(os.path.basename(globs['st-input']))[0] + "-cf-stats.tab");
            # CF files before pruning

            if not globs['use-labels'] or globs['scf-quartets']:
                globs['st-final-file'] = os.path.join(globs['outdir'], os.path.splitext(os.path.basename(globs['st-input']))[0] + "-cf-pruned.tre");
                globs['st-cf-stat-prune-file'] = os.path.join(globs['outdir'], os.path.splitext(os.path.basename(globs['st-input']))[0] + "-cf-pruned-stats.tab");
            else:
                globs['st-final-file'] = os.path.join(globs['outdir'], os.path.splitext(os.path.basename(globs['st-input']))[0] + "-pruned.tre");
            # CF files after pruning - tree file is also written as final tree

        else:
            globs['st-cf-tree-file'] = os.path.join(globs['outdir'], "species-tree-cf.tre");
            globs['st-cf-stat-file'] = os.path.join(globs['outdir'], "species-tree-cf-stats.tab");
            # CF files before pruning

            if not globs['use-labels'] or globs['scf-quartets']:
                globs['st-final-file'] = os.path.join(globs['outdir'], "pruned-species-tree-cf.tre");
                globs['st-cf-stat-prune-file'] = os.path.join(globs['outdir'], "pruned-species-tree-cf-stats.tab");
            else:
                globs['st-final-file'] = os.path.join(globs['outdir'], "pruned-species-tree.tre");
        
        ## Output files based on species tree input
        ###############
            
        if args.aln_stats_flag:
            globs['aln-stat-file'] = os.path.join(globs['outdir'], os.path.splitext(os.path.basename(globs['aln-dir']))[0] + "-aln-stats.csv");
            # Input align stats csv file

            globs['aln-stat-prune-file'] = os.path.join(globs['outdir'], os.path.splitext(os.path.basename(globs['aln-dir']))[0] + "-aln-stats-pruned.csv");
            # Pruned align stats csv file

            globs['aln-pool'] = mp.Pool(processes=globs['num-procs']);
            # Set the alignment stat pool here before we read anything large into memory

            if not globs['cf-only']:
                globs['aln-pool-prune'] = mp.Pool(processes=globs['num-procs']);
                # Set the alignment stat pool for pruned alns here before we read anything large into memory

        ## Output files for alignment stats
        ###############

        if globs['gt-input']:
            globs['prune-gt-file'] = os.path.join(globs['outdir'], os.path.splitext(os.path.basename(globs['gt-input']))[0] + "-pruned.tre");

            globs['topo-count-file'] = os.path.join(globs['outdir'], os.path.splitext(os.path.basename(globs['gt-input']))[0] + "-topo-counts.tab");
            globs['topo-count-pruned-file'] = os.path.join(globs['outdir'], os.path.splitext(os.path.basename(globs['gt-input']))[0] + "-pruned-topo-counts.tab");


        if globs['aln-dir']:
            globs['prune-aln-dir'] = os.path.join(globs['outdir'], os.path.splitext(os.path.basename(globs['aln-dir']))[0] + "-pruned");
            os.makedirs(globs['prune-aln-dir'], exist_ok=True);

        ## Output files for pruned gene trees and alignments (if given)
        ###############

        globs['run-name'] = os.path.basename(os.path.normpath(globs['outdir']));
        globs['logfilename'] = os.path.join(globs['outdir'], globs['run-name'] + ".log");
        # Log file

        if not args.append_log_flag and not globs['norun'] and not globs['debug-tree']:
            logfile = open(globs['logfilename'], "w");
            logfile.write("");
            logfile.close();
        # Prep the logfile to be overwritten

        ## Set the log file
        ###############
    ## End --labeltree check

    else:
        globs['log-v'] = -1;
    # If --labeltree is set, turn off log reporting

    ###############

    if args.quiet_flag:
        globs['quiet'] = True;
    # Check the quiet option

    ###############

    if globs['psutil']:
        globs['pids'] = [psutil.Process(os.getpid())];
    # Get the starting process ids to calculate memory usage throughout.

    ###############

    startProg(globs);
    # After all the essential options have been set, call the welcome function.

    return globs;

#############################################################################

def startProg(globs):
# A nice way to start the program.
    print("#");
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# Tree pruning with heuristic metrics.");
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# Version " + globs['version'] + " released on " + globs['releasedate']);
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# bonsai was developed by " + globs['authors']);
    #CORE.printWrite(globs['logfilename'], globs['log-v'], "# Citation:      " + globs['doi']);
    #CORE.printWrite(globs['logfilename'], globs['log-v'], "# Website:       " + globs['http']);
    #CORE.printWrite(globs['logfilename'], globs['log-v'], "# Report issues: " + globs['github']);
    CORE.printWrite(globs['logfilename'], globs['log-v'], "#");
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# The date and time at the start is: " + CORE.getDateTime());
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# Using Python version:              " + globs['pyver'] + "\n#");
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# The program was called as:         " + globs['call'] + "\n#");

    if globs['info']:
        return;
    # If --info is set, return after printing program info

    #######################

    pad = 40;
    opt_pad = 30;
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# INPUT INFO:");

    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Species tree input type:", pad) + globs['st-input-type']);
    if globs['st-input-type'] == "file":
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Species tree file", pad) + globs['st-input']);

    if globs['debug-tree']:
        print("\n--debugtree SET. READING TREE AND EXITING.\n");
        return;       

    if globs['label-tree']:
        print("\n--labeltree SET. READING TREE AND EXITING.\n");
        return;
    # If --labeltree is set, print a message and return.

    if globs['gt-input']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Gene tree input file", pad) + globs['gt-input']);
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Pruned gene tree file", pad) + globs['prune-gt-file']);
    if globs['aln-dir']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Alignment directory", pad) + globs['aln-dir']);

    # Input
    #######################

    CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# OUTPUT INFO:");

    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Output directory:", pad) + globs['outdir']);

    if globs['gt-input']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Pruned gene tree file", pad) + globs['prune-gt-file']);
    if globs['aln-dir'] and not globs['cf-only']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Pruned alignment directory", pad) + globs['prune-aln-dir']);

        if globs['aln-stat-file']:
            CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Alignment stats file", pad) + globs['aln-stat-file']);
            if not globs['cf-only']:
                CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Pruned alignment stats file", pad) + globs['aln-stat-file']);

    if globs['gt-input']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Concordance tree file", pad) + globs['st-cf-tree-file']);
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Concordance stat file", pad) + globs['st-cf-stat-file']);
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Topology count file", pad) + globs['topo-count-file']);
        
        if not globs['cf-only']:
            CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Pruned concordance stat file", pad) + globs['st-cf-stat-prune-file']);
            CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Pruned topology count file", pad) + globs['topo-count-pruned-file']);

    if not globs['cf-only']:
         CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Final pruned species tree file", pad) + globs['st-final-file']);

    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Log file:", pad) + os.path.basename(globs['logfilename']));

    if globs['cf-only']:
        print("\n--cf SET. CALCULATING CONCORDANCE FACTORS AND EXITING.\n");
        return;
    # If --cf is set, print a message and return.

    # Output
    #######################

    CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# OPTIONS INFO:");    
    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Option", pad) + CORE.spacedOut("Current setting", opt_pad) + "Current action");      

    if globs['use-labels']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --labels", pad) + 
                CORE.spacedOut(str(globs['use-labels']), opt_pad) + 
                "Labels on the input species tree will be used for paring as opposed to calculation gCF.");
    # Reporting labels option 

    if globs['scf-quartets']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# -scf", pad) +
                CORE.spacedOut(str(globs['scf-quartets']), opt_pad) + 
                "This many quartets will be sampled per branch to calculate sCF");        
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# -n", pad) +
                CORE.spacedOut(str(globs['num-procs']), opt_pad) + 
                "This many processes will be used to calculate sCF");
    elif globs['aln-dir']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# -scf", pad) +
                CORE.spacedOut(str(globs['scf-quartets']), opt_pad) + 
                "sCF will not be calculated"); 
    # Reporting sCF options          

    if not globs['prune-file']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Branch length percentile (-b)", pad) + 
                    CORE.spacedOut(str(globs['bl-percentile']), opt_pad) + 
                    "Branches will be considered for paring each iteration if they fall below this percentile of all branch lengths.");
        # Reporting the lower branch length percentile option

        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Lower gCF threshold (-g)", pad) + 
                CORE.spacedOut(str(globs['gcf-threshold']), opt_pad) + 
                "Branches will be considered for paring each iteration if they have a gCF lower than this.");
        # Reporting the gcf threshold option

        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Maximum iterations of paring (-i)", pad) + 
            CORE.spacedOut(str(globs['max-iterations']), opt_pad) + 
            "The program will stop paring after this many iterations, regardless of other settings.");
        # Reporting the max iters option

        total_max_spec_str = "Unlimited";
        if globs['total-max-spec'] != 999999999:
            total_max_spec_str = str(globs['total-max-spec']);
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Total maximum tips to prune (-s)", pad) + 
            CORE.spacedOut(total_max_spec_str, opt_pad) + 
            "The program will stop paring after this tips have been pruned, regardless of other settings.");
        # Reporting the total max species option

        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Maximum tips to prune per branch (-m)", pad) + 
            CORE.spacedOut(str(globs['branch-max-spec']), opt_pad) + 
            "Branches selected to prune that have more than this many species descending will be retained.");
        # Reporting the branch max species option

        if globs['exempt-file']:
            CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# -e", pad) +
                        CORE.spacedOut(globs['exempt-file'], opt_pad) + 
                        "bonsai will exempt the branches listed in this file from pruning.");
        # Reporting the exempt file option.
    ## Paring block, if no prune file is given

    else:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# -e", pad) +
                    CORE.spacedOut(globs['prune-file'], opt_pad) + 
                    "bonsai will only prune the branches listed in this file (no iterative pruning done).");
    # Reporting the prune file option. 

    if globs['overwrite']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --overwrite", pad) +
                    CORE.spacedOut("True", opt_pad) + 
                    "bonsai will OVERWRITE the existing files in the specified output directory.");
    # Reporting the overwrite option.

    if not globs['quiet']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --quiet", pad) + 
                    CORE.spacedOut("False", opt_pad) + 
                    "Time, memory, and status info will be printed to the screen while pare is running.");
    else:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --quiet", pad) + 
                    CORE.spacedOut("True", opt_pad) + 
                    "No further information will be printed to the screen while pare is running.");
        CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
        CORE.printWrite(globs['logfilename'], globs['log-v'], "# Running...");
    # Reporting the quiet option.

    # if globs['debug']:
    #     CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --debug", pad) + 
    #                 CORE.spacedOut("True", opt_pad) + 
    #                 "Printing out a bit of debug info.");
    # Reporting the debug option.

    if globs['norun']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --norun", pad) + 
                    CORE.spacedOut("True", opt_pad) + 
                    "ONLY PRINTING RUNTIME INFO.");
        CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    # Reporting the norun option.

    # Other options
    #######################

#############################################################################