
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
num_procs = _env_int_or_none('SLURM_NPROCS')

is_slurm = procid is not None

if is_slurm:
    print('Started Osprey on SLURM node %d/%d' % (procid + 1, num_procs))


def launch(num_nodes, parallelism, mem_mib, python=sys.executable, srun_args=[]):

    # convert MiB to MB for slurm
    # NOPE: empirically, it looks like SLURM uses MiB, GiB, etc after all
    #mem_mb = int(math.ceil(mem_mib*1024*1024/1000/1000))

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
        '--mem=%dM' % mem_mb,
        '--export=ALL,OSPREY_PREAMBLE=false'
    ]
    if parallelism.numGpus > 0:
        args += ['--gres=gpu:%d' % parallelism.numGpus]

    # call that command with srun to start the cluster job
    srun_cmd = ['srun'] + args + srun_args + script_cmd
    subprocess.run(srun_cmd, check=True)
