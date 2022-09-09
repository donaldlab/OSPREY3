

'''
This module handles interactiing with [SLURM](https://slurm.schedmd.com/overview.html) clusters.
'''


import os
import sys
import subprocess
import math


def _env_or_none(key):
    try:
        return os.environ[key]
    except KeyError:
        return None


def _env_int_or_none(key):
    val = _env_or_none(key)
    if val is not None:
        val = int(val)
    return val


procid = _env_int_or_none('SLURM_PROCID')
'''
`int`: The current processor id of the SLURM job, or `None` if no SLURM job is currently running
'''

num_procs = _env_int_or_none('SLURM_NPROCS')
'''
`int`: The number of processors in the SLURM job, or `None` if no SLURM job is currently running
'''


is_slurm = procid is not None
'''
`bool`: `True` if a SLURM job is currently running, `False` if not.
'''


if is_slurm:
    print('Started Osprey on SLURM node %d/%d' % (procid + 1, num_procs))


def launch(num_nodes, parallelism, mem_mib, python=sys.executable, srun_args=[]):
    '''
    Launches the currently-running Python script as a SLURM job.
    The top-level python script, along with all its command-line arguments (ie `sys.argv`),
    will be launched on the SLURM cluster using `srun`.

    This function should be called on the SLURM login node.

    # Arguments
    num_nodes `int`: The number of SLURM nodes to use in the computation, ie the `--ntasks` flag.
    parallelism ${type_java(.parallelism.Parallelism)}: The compute resources to request for each node, ie the `--cpus-per-task` and `--gpus-per-task` flags.
    mem_mib `int`: The amount of memory to request for each node, in [Mebibytes](https://simple.wikipedia.org/wiki/Mebibyte), ie the `--mem` flag.
    python `str`: The path to the python interpreter, defaults to the currently-running python interpreter.
    srun_args `[str]`: Additional arguments to pass to `srun`.
    '''

    # convert MiB to MB for slurm
    # NOPE: empirically, it looks like SLURM uses MiB, GiB, etc after all
    #mem_mib = int(math.ceil(mem_mib*1024*1024/1000/1000))

    # python buffers stdout by default,
    # so flush it before spawning another process,
    # to keep the console in order
    sys.stdout.flush()

    # build the command to call the top-level script
    script_cmd = [python] + sys.argv

    # convert the resource requirements into arguments to srun
    args = [
        '--ntasks=%d' % num_nodes,
        '--cpus-per-task=%d' % parallelism.numThreads,
        '--mem=%dM' % mem_mib,
        '--export=ALL,OSPREY_PREAMBLE=false'
    ]
    if parallelism.numGpus > 0:
        args += ['--gres=gpu:%d' % parallelism.numGpus]

    # call that command with srun to start the cluster job
    srun_cmd = ['srun'] + args + srun_args + script_cmd
    subprocess.Popen(srun_cmd)
