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
        'authors' : "Jeffrey Good, Gregg Thomas",
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

        'tree-input' : False,
        'tree-input-type' : False,
        'orig-tree-str' : False,
        'tree-dict' : False,
        'labeled-tree-str' : False,
        'tips' : False,
        'root' : False,
        # Tree info

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
        'run-name' : 'pare',
        'logfilename' : 'pare.errlog',
        'logdir' : '',
        'overwrite' : False,
        # I/O options

        'label-tree' : False,
        'info' : False,
        'dryrun' : False,
        'quiet' : False,
        # Other user options

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
        'nolog' : False,
        # Internal stuff
    }

    globs_init['logfilename'] = "pare-" + globs_init['startdatetime'] + ".errlog";
    # Add the runtime to the error log file name.

    globs = StrictDict(globs_init);
    # Restrict the dict from having keys added to it after this

    return globs;

#############################################################################