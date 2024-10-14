# this files submits all IAS jobs to grisman and saves outputs for each residue to a csv
# also returns the consensus sequence


import glob
import os
import math
import subprocess
import time
import shutil
import re
import csv
from collections import OrderedDict
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


# give permissions to all files
os.system("chmod -R +x *")


# get the number of residues (same for each backbone)
def how_many_residues():

    num_residues = 0

    for d in glob.glob("match*"):

        res_folders = d + "/*"
        for f in glob.glob(res_folders):
            if "omol" in f:
                continue
            num_residues += 1

        return num_residues


# get the current number of matches (can change during runtime)
def how_many_matches():

    num_matches = 0

    for d in glob.glob("match*"):
        num_matches += 1

    return num_matches


# cancel job if run too long
def adjust_stat(job_id, time_entry, filename):

    # don't cancel this job
    if "match" not in filename:
        return

    num_fields = 0
    for s in time_entry:
        if s == ":":
            num_fields += 1

    # cancel if >= 8 hours
    if num_fields >= 2:
        num_hours = int(time_entry.split(":")[0])
        if num_hours >= 12:
            print("%s took >= 12 hours, so cancelling this job and marking for deletion" % filename)
            subprocess.run(['scancel', job_id], stdout=subprocess.PIPE)

            # add to log to be deleted
            short_file = filename[:-3]
            format_file = short_file.replace('-', '/')
            took_too_long.append(format_file)


# check squeue + cancel if needed
def check_queue():

    # get the id (%A), time (%M), and name (of the sh file, %j), -h for no header
    status = subprocess.run(['squeue', '--me', '-h', '--format=%A,%M,%j'], stdout=subprocess.PIPE)
    str_status = status.stdout.decode('utf-8')

    # check each running job
    for line in str_status.splitlines():
        data = line.split(",")
        adjust_stat(data[0], data[1], data[2])


# check how many jobs (including this one) are still running
def jobs_running():

    status = subprocess.run(['squeue',  '--me', '-h'], stdout=subprocess.PIPE)
    all_jobs = status.stdout.decode('utf-8')

    num_jobs = 0
    for line in all_jobs.splitlines():
        num_jobs += 1

    return num_jobs


# convert a list of dictionaries for a given residue to ddG
def convert_ddG(s_list, fname):

    wt_score = 1001

    for entry in s_list:
        if entry["RESID"].islower():
            wt_score = entry["DDG"]
            break

    if wt_score == 1001:
        print("ERROR! Was unable to find WT score for %s" % fname)
        exit()

    raw_wt = pow(10, wt_score)
    dg_wt = -1 * 8.3145 * 298 * math.log(raw_wt)
    dg_wt_kcal = dg_wt * (1 / 4184)

    for entry in s_list:
        curr_score = entry["DDG"]
        raw_curr = pow(10, curr_score)
        dg_curr = -1 * 8.3145 * 298 * math.log(raw_curr)
        dg_curr_kcal = dg_curr * (1 / 4184)
        entry["DDG"] = (dg_wt_kcal - dg_curr_kcal)

    return s_list


# parse a kstar output and return as a list w/ dictionaries for each residue
def parse_kstar(filename, res_num, match_num, res_type):

    score_list = []

    f = open(filename, "r")

    for line in f:

        tsv = line.split("   ")

        if tsv[0] == "sequence":

            curr_mut = "ERROR"
            curr_score = 1001

            info = re.findall("[^ ]+ [^ ]+", line)

            res_regex = "^" + res_num + " " + res_type
            for entry in info:
                res_info = re.search(res_regex, entry)
                k_info = re.search("log10", entry)
                if res_info:
                    curr_mut = entry[-3:]
                elif k_info:
                    curr_score = entry[11:]
                    if curr_score == 'none':
                        curr_score = 0
                    else:
                        curr_score = float(curr_score)
                    break

            if curr_mut == "ERROR":
                print("ERROR! Was unable to find a mutant id for %s" % filename)
                exit()
            if curr_score == 1001:
                print("ERROR! Was unable to find a kstar score for %s" % filename)
                exit()

            score_data = {}
            score_data["MATCH"] = match_num
            score_data["RESNUM"] = res_num
            score_data["RESID"] = curr_mut
            score_data["DDG"] = curr_score  # we'll convert K* -> ddG soon
            score_list.append(score_data)

    ddg_data = convert_ddG(score_list, filename)

    return ddg_data


# for each residue make a csv with ddG data
def get_ddG():

    all_scores = []

    for d in glob.glob("match*/*/"):

        matchname = d.split('/')[0]
        resnum = d.split('/')[1][3:]
        restype = d.split('/')[1][:3]

        filepath = d + 'submit.out'
        ddg_scores = parse_kstar(filepath, resnum, matchname, restype)
        all_scores += ddg_scores

    return all_scores


# make a graph from residue ddG data
def graph_ddG(csv_filename):

    # open the csv and load as df
    df = pd.read_csv(csv_filename, index_col=0)

    # uncomment if you want to include WT in graph
    # df['RESID'] = df["RESID"].str.isupper()

    # make the directory for saving graphs
    Path("ddg_graphs").mkdir(exist_ok=True)

    for i in range(1, (number_residues + 1)):

        # exclude the WT, so we're not ranking scores of 0 (ddG for WT is always 0)
        res_info = df.loc[(df['RESNUM'] == i) & (df['RESID'].str.isupper())]
        # res_info.plot.scatter(x='RESID', y='DDG', figsize=[9, 6])
        fig = res_info.boxplot(by='RESID', column=['DDG'], figsize=[11, 6], fontsize=10)

        # set title as residue number
        resname = "RESIDUE_%s" % i
        fig.get_figure().suptitle(resname)
        fig.set_ylabel("ΔΔG (kcal/mol)")
        fig.set_xlabel("Residue Identity")

        # save the figure as a png file
        filename = "ddg_graphs/" + resname + ".png"
        plt.savefig(filename, dpi=300)
        print("saved %s" % filename)


# check + delete if any errors occurred
def check_for_errors():

    failed_matches = []

    for d in glob.glob("match*/*/"):

        fpath = d + 'submit.out'

        try:
            f = open(fpath, "r")
        except:
            print("Error found: no log file for %s. Will delete this match." % d)
            del_dir = d.split('/')[0]
            if del_dir not in failed_matches:
                failed_matches.append(del_dir)
            continue

        for line in f:
            if "Error" in line or "exception" in line:
                print("The following error occured in %s. Will delete this match." % d)
                print(line + '\n')
                del_dir = d.split('/')[0]
                if del_dir not in failed_matches:
                    failed_matches.append(del_dir)
                break

    if len(failed_matches) == 0:
        print("No matches reported failures")

    else:
        print("Now deleting matches that threw errors")

        for f in failed_matches:
            shutil.rmtree(f)
            print("deleted: %s" % f)


# find the runtime of each residue
def find_runtime():

    runtime_dict = OrderedDict()

    for r in range(1, (number_residues + 1)):
        resname = "RES" + str(r)
        runtime_dict[resname] = 0

    for d in glob.glob("match*/*/submit.out"):

        resnum = "RES" + d.split('/')[1][3:]

        f = open(d, 'r')
        for line in f:
            if "KStar runtime" in line:
                new_runtime = float(line.split()[3])
                try:
                    old_runtime = runtime_dict[resnum]
                except:
                    old_runtime = 0
                runtime_dict[resnum] = (new_runtime + old_runtime)

    return runtime_dict


# find the # residues
number_residues = how_many_residues()

# submit all jobs
for d in glob.glob("match*/*/"):
    os.chdir(d)
    os.system("sbatch *sh")
    os.chdir("../..")


# check runtimes until all jobs finished
start_time = time.time()
took_too_long = []
while jobs_running() > 1:
    end_time = time.time()
    elapsed_time = (end_time - start_time) / 60
    print("Jobs running/queued: %s <elapsed time> %s minutes" % (jobs_running(), round(elapsed_time, 3)))
    check_queue()
    time.sleep(120)


# delete directories that took too long
print("\n---The following residues took too long, and will be deleted:")
for t in took_too_long:
    shutil.rmtree(t)
    print("deleted: %s" % t)


# delete any matches that had runtime errors
print("\n---All runs complete! Now checking logs for errors---")
check_for_errors()


# get total runtime across all residues
print("\n---Finding total runtime for each residue---")
res_runtimes = find_runtime()
total_runtime = 0
for key, value in res_runtimes.items():
    print(key, ":", value, "seconds")
    total_runtime += value
print("Total runtime was %s seconds" % total_runtime)


# get ddg data from submit.out files
print("\n---Fetching score data from logs---\n")
g_data = get_ddG()


# write all ddG to a csv
print("\n---Writing to ddG_data.csv---\n")
with open('ddG_data.csv', 'w', newline='') as csvfile:
    fieldnames = ['MATCH', 'RESNUM', 'RESID', 'DDG']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(g_data)


# get ddg for each residue + graph
print("\n---Graphing each residue from the csv---\n")
graph_ddG('ddG_data.csv')
