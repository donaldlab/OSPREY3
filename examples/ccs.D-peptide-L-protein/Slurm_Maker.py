# this script automates the process of creating slurm scripts for submitting jobs to the cluster
# intended for use in DL.py

def make_slurm(match_num: str, memory: int, cpus: int, partition: str, osprey_path: str, epsilon: float):

    # specify slurm params
    first_line = "#!/bin/csh -f\n"
    second_line = "#SBATCH --output=" + match_num + ".out\n"
    third_line = "#SBATCH --mem=" + str(memory) + "G\n"
    fourth_line = "#SBATCH --cpus-per-task=" + str(cpus) + "\n"
    fifth_line = "#SBATCH -A --partition=" + partition + "\n"

    # customize Kstar command to be run on cluster
    kstar_cmd = "time " + osprey_path + " bbkstar"
    complex = " --complex-confspace ./" + match_num + "-complex.ccsx"
    target = " --target-confspace ./" + match_num + "-target.ccsx"
    peptide = " --design-confspace ./" + match_num + "-peptide.ccsx"
    ensemble = " --ensemble-dir ./ensembles"
    eps = " -e " + str(epsilon)
    max_mut = " --max-simultaneous-mutations 100000\n"

    # combine + cleanup
    slurm_opener = first_line + second_line + third_line + fourth_line + fifth_line
    slurm_kstar = kstar_cmd + complex + target + peptide + ensemble + eps + max_mut
    cleanup = "rm *confdb\necho \"completed\"\n"

    slurm_str = slurm_opener + slurm_kstar + cleanup

    return slurm_str