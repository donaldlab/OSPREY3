# this files submits all IAS jobs to grisman and saves outputs for each residue to a csv
# also returns the consensus sequence

import glob
import os
import csv
import time
import subprocess

# give permissions to all files
os.system("chmod -R +x *")

# find how many residues for each scan
num_residues = 0
for d in glob.glob("match*"):

    res_folders = d + "/*"
    for f in glob.glob(res_folders):
        if "omol" in f:
            continue
        num_residues += 1
    break

# find number of matches
num_matches = 0
for d in glob.glob("match*"):
    num_matches += 1


# create job status template
def make_status_template():

    all_status = []

    for d in glob.glob("match*"):
        stat_list = [False for i in range(num_residues)]
        match_data = {}
        match_data[d] = stat_list
        all_status.append(match_data)

    return all_status


# check and update the status of a match
def update_status(job_dict):

    matchname = next(iter(job_dict))
    filepath = matchname + "/*/submit.out"

    # if running, check the log file + update job status (if done)
    for s in glob.glob(filepath):
        f = open(s, "r")
        for line in f:
            if "completed" in line:
                holder = s.split("/")[1]
                resnum = int(holder[3:]) - 1
                job_dict[matchname][resnum] = True
                break


# get the best mutant (w/ score) for a given residue number in csv format
def make_res_csv(resnum):

    res_data = []

    # parse the scores using delimiters
    for d in glob.glob("match*/*/"):

        matchname = d.split('/')[0]

        mut = ""
        kstar = 0

        holder = d.split("/")[1]
        res_number = int(holder[3:])

        if res_number == resnum:
            filepath = d + "submit.out"
            f = open(filepath, "r")

            for line in f:
                tsv = line.split("   ")

                if tsv[0] == "sequence":

                    # get mutant id
                    ssv = tsv[2].split(" ")
                    for i in range(0, len(ssv)):
                        if ssv[i] == str(resnum):
                            curr_mut = ssv[(i + 1)][4:]

                    # get kstar score
                    ssv = tsv[3].split(" ")
                    curr_score = float(ssv[1])

                    if curr_score > kstar:
                        kstar = curr_score
                        mut = curr_mut

            # add the data
            best_data = [matchname, resnum, mut, kstar]
            res_data.append(best_data)

    # output as csv
    print("\nresidue " + str(resnum) + " complete! Saving CSV now.\n")
    fields = ["match_number", "residue_number", "residue_type", "kstar_score"]
    filename = str(resnum) + "_score.csv"
    with open(filename, 'w') as w:
        write = csv.writer(w)
        write.writerow(fields)
        write.writerows(res_data)


# check if any residues have finished, and update the overall log if so
def check_res_status(job_list):

    curr_res = 0

    while curr_res < num_residues:

        # don't recheck finished residues
        if overall_stats[curr_res]:
            curr_res += 1
            continue

        # check status
        num_finished = 0
        for entry in job_list:
            for key, value in entry.items():
                if value[curr_res]:
                    num_finished += 1
        if num_finished == num_matches:
            overall_stats[curr_res] = True
            make_res_csv(curr_res)

        curr_res += 1


# submit all jobs
for d in glob.glob("match*/*/"):
    os.chdir(d)
    os.system("sbatch *sh")
    os.chdir("../..")

# track + report progress
overall_stats = [False for i in range(num_residues)]
logs = make_status_template()

while False in overall_stats:

    for entry in logs:
        update_status(entry)

    check_res_status(logs)
