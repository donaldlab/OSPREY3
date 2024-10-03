# this script serves to perform gradient descent to tune OSPREY's energy function on Thanatin binding data
# experimental binding data from https://www.sciencedirect.com/science/article/pii/S2590152423000077

# this script must be run within the OSPREY3 directory: https://github.com/donaldlab/OSPREY3 (must be nom branch)

# required bash scripts:
# ktune.sh: kstar specs for cluster
# submit_jobs.sh: copies ktune.sh into each run and submits to cluster
# get_KStar.sh: checks all runs and returns csv file with kstar score data for PM and variant (mutant or ortholog)

import math
import time
import subprocess
import os
import matplotlib.pyplot as plt
import scipy as sc
import glob
from collections import OrderedDict

import osprey
osprey.start()
import osprey.prep


# helper conversion functions
def kstar_to_ddg(score1, score2):
    # undo log10 for both
    raw_score1 = pow(10, score1)
    raw_score2 = pow(10, score2)

    # delta G = -RT ln K*
    dg1 = -1 * 8.3145 * 298 * math.log(raw_score1)
    dg2 = -1 * 8.3145 * 298 * math.log(raw_score2)

    # subtract: PM - ortholog/mutant
    return dg1 - dg2


def kd_to_ddg(score1, score2):
    # convert to Ka = 1 / Kd
    ka_score1 = (1 / score1)
    ka_score2 = (1 / score2)

    # delta G = - RT ln Ka
    dg1 = -1 * 8.3145 * 298 * math.log(ka_score1)
    dg2 = -1 * 8.3145 * 298 * math.log(ka_score2)

    # subtract: PM - ortholog/mutant
    return dg1 - dg2


# OSPREY confspace compiler
def compile(filename, ccsxname):

    compiler = osprey.prep.ConfSpaceCompiler(filename)

    # add the forcefields to use
    compiler.getForcefields().add(osprey.prep.Forcefield.Amber96)
    compiler.getForcefields().add(osprey.prep.Forcefield.EEF1)

    # run the compiler and wait for it to finish
    print('compiling %s ...' % filename)
    progress = compiler.compile()
    progress.printUntilFinish(5000)
    report = progress.getReport()

    # check for compiler errors
    if report.getError() is not None:
        raise Exception('Compilation failed', report.getError())

    # save the compiled conf space
    print('saving to %s ...' % ccsxname)
    open(ccsxname, 'wb').write(osprey.prep.saveCompiledConfSpace(report.getCompiled()))


# get kstar scores, or status every 60 seconds until completion
def get_kstars():

    kstar_done = False
    while not kstar_done:
        kstar_status = subprocess.run(["./get_Kstar.sh"], stdout=subprocess.PIPE)
        if "STILL_RUNNING" not in kstar_status.stdout.decode('utf-8'):
            return kstar_status.stdout.decode('utf-8')
        else:
            time.sleep(60)


# convert kstar csv values into python dictionary
def convert_kstars(score_str):

    kstar_dictionary = OrderedDict()

    for line in score_str.splitlines():
        data = line.split(",")
        name = data[0]
        pmscore = data[1]
        mutscore = data[2]
        scores = [pmscore, mutscore]
        kstar_dictionary[name] = scores

    return kstar_dictionary


# set up a new step
def setup_gd():

    # build a whl (default params)
    # NOTE: you cannot call gradlew from another directory. Must be in OSPREY3.
    build_whl = os.system("./gradlew pythonWheel")

    # install it (this can be done from any directory as long as filepath to whl is correct)
    # be sure to have jdk19 in your bash profile so gradle can find it
    install_whl = os.system("pip install build/python/osprey-3.3-py3-none-any.whl --force-reinstall")

    # compile each confspace (these hold all the forcefield info)
    for cf_file in glob.glob("PM-*/*/*confspace"):

        cf = osprey.prep.loadConfSpace(open(cf_file, 'r').read())
        cf_ccsx = cf_file[:-9] + ".ccsx"

        with osprey.prep.LocalService():
            compile(cf, cf_ccsx)


# get a round of scores for gradient descent
def run_gd():

    # # submit to cluster
    # subprocess.call(". ./submit_jobs.sh")

    # check scores every minute, save to csv when completed
    kstar_scores = get_kstars()

    # parse them
    kstar_dict = convert_kstars(kstar_scores)

    # convert each key to a single ddg value and add to list that matches binding data
    combined_scores = OrderedDict()
    for key, value in kstar_dict.items():
        # get kstar ddg
        new_ddg = kstar_to_ddg(float(value[0]), float(value[1]))

        # using key, combine ddg kstar with ddg binding

    # calculate Pearson
    # TODO loop over combined dictionary for this
    # corr, pvalue = sc.stats.pearsonr(kstar_data, binding_data)


run_gd()

# # visualize it
# plt.scatter(defaultKStar, bindingdata)
# plt.plot(defaultKStar, bindingdata)
# plt.xlabel("delta Kstar scores")
# plt.ylabel("delta Ka value")
#
# plt.show()

