#############################################################################
# This file holds some global variables for some of the input options.
# These global parameters should be read only -- they are not modified anywhere 
# else in the code except when reading the input options.
#
# This dictionary will also be used to track output.
#############################################################################

import sys
import timeit
import lib.core as CORE

#############################################################################

class StrictDict(dict):
# This prevents additional keys from being added to the global params dict in
# any other part of the code, just to help me limit it's scope
# https://stackoverflow.com/questions/32258706/how-to-prevent-key-creation-through-dkey-val
    def __setitem__(self, key, value):
        if key not in self:
            raise KeyError("{} is not a legal key of this StricDict".format(repr(key)));
        dict.__setitem__(self, key, value);

#############################################################################

def init():
    globs_init = {
        'version' : 'Beta 1.0',
        'releasedate' : "October 2021",
        'authors' : "Jeff Good, Gregg Thomas",
        'doi' : '',
        'http' : '',
        'github' : '',
        'starttime' : timeit.default_timer(),
        'startdatetime' : CORE.getOutTime(),
        # Meta info

        'pyver' :  ".".join(map(str, sys.version_info[:3])),
        # System info

        'call' : "",
        # Script call info

        'st-input' : False,
        'st-input-type' : False,
        'orig-st-str' : False,
        'st' : False,
        'use-labels' : False,
        'st-final-file' : False,
        # Species tree info

        'gt-input' : False,
        'gt-input-empty' : 0,
        'gt-input-skipped' : 0,
        'gts' : {},
        'gt-current' : False,
        # Gene tree params

        'aln-dir' : False,
        'alns' : False,
        'seq-compression' : False,
        'num-alns' : False,
        'aln-skip-chars' : ["-", "N", "n"],
        'aln-stats' : {},
        'aln-stats-pruned' : {},
        'aln-stat-file' : False,
        'aln-stat-prune-file' : False,
        'aln-pool' : False,
        'aln-pool-prune' : False,
        # Alignment params
        
        'cf-only' : False,
        'gcf-dict' : {},
        'scf-dict' : {},
        'scf-quartets' : False,
        'scf-pool' : False,
        'scf-pool-prune' : False,
        'seed' : False,
        'st-cf-tree-file' : False,
        'st-cf-stat-file' : False,
        'st-cf-stat-prune-file' : False,
        # CF options

        #'count-only' : False,
        'topo-count-file' : False,
        'topo-count-pruned-file' : False,
        # Topology counting options

        'prune-file' : False,
        'prune-branches' : [],
        'prune-clades' : [],
        'prune-st-file' : False,
        'prune-gt-file' : False,
        'prune-aln-dir' : False,
        # Specified branches to prune

        'exempt-file' : False,
        'exempt-branches' : [],
        'exempt-clades' : [],
        # Excluded branches

        'iter-trees' : [],
        'bl-thresholds' : [],
        'pared-branches' : [],
        'pruned-tips' : [],
        'total-pruned-tips' : 0,
        'total-pared-branches' : 0,
        # Tracking

        'branch-max-spec' : 10,
        'bl-percentile' : 10,
        'gcf-threshold' : 25,
        # Pruning thresholds

        'max-iterations' : 3,
        'total-max-spec' : 999999999,
        # Stopping conditions

        'outdir' : '',
        'run-name' : 'bonsai',
        'logfilename' : 'bonsai.errlog',
        'logdir' : '',
        'overwrite' : False,
        # I/O options

        'label-tree' : False,
        'info' : False,
        'dryrun' : False,
        'quiet' : False,
        # Other user options

        'num-procs' : 1,
        'pad' : 82,
        'endprog' : False,
        'warnings' : 0,
        'exit-code' : 0,
        'log-v' : 1,
        'stats' : True,
        'progstarttime' : 0,
        'stepstarttime' : 0,
        'pids' : "",
        'psutil' : False,
        'qstats' : False,
        'norun' : False,
        'debug' : False,
        'debug-tree' : False,
        'nolog' : False,
        # Internal stuff
    }

    globs_init['logfilename'] = "bonsai-" + globs_init['startdatetime'] + ".errlog";
    # Add the runtime to the error log file name.

    globs = StrictDict(globs_init);
    # Restrict the dict from having keys added to it after this

    return globs;

#############################################################################