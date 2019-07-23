import sys
sys.path.insert(0, '..')

import argparse
import json
import subprocess 
import re
import ray
from pathlib import Path
from model.SBML import SBML
from experiment.ExperimentSet import ExperimentSet
from model.SBMLtoODES import sbml_to_odes
from marginal_likelihood.MarginalLikelihood import MarginalLikelihood
from model.PriorsReader import define_sbml_params_priors

def prepare_workers (workers, ray_path):
    current_path = str (Path ().absolute ())
    start_cluster_head_bin = current_path + "/start_cluster_head.sh"
    start_worker_bin = current_path + "/start_cluster_worker.sh"

    # Start the main node and the Redis server
    print ("Starting", workers[0], " as head node of the cluster...", \
            end="", flush=True)
    proc = subprocess.Popen ([start_cluster_head_bin, workers[0], \
            ray_path], stdout=subprocess.PIPE)
    (out, err) = proc.communicate ()

    redis_address_pattern = \
            re.compile (r".*redis-address (\d+\.\d+\.\d+\.\d+:\d+)")
    try:
        m = redis_address_pattern.match (str (out))
        redis_address = m.group (1)
        print ("[OK]")
    except:
        raise ConnectionError ("Could not find redis server address." \
                + " Are you sure you provided a valid worker machine?")

    # Now start ray on every other working machine and set the Redis 
    # server
    print ("Redis server address is:", redis_address)
    for worker in workers[1:]:
        print ("Starting worker " + worker + "...", end="", flush=True)
        proc = subprocess.Popen ([start_worker_bin, worker, ray_path, \
                redis_address], stdout=subprocess.PIPE)
        (out, err) = proc.communicate ()
        worker_up_pattern = re.compile (".*Started Ray on this node")
        if worker_up_pattern.search (str (out)):
            print ("[OK]")
        else:
            print ("Failed!")
            print ("Warning: Could not start ray on " + worker + ".")
            print (out.decode ("utf-8"))
            if err is not None:
                print (err.decode ("utf-8"))
    return redis_address


@ray.remote
def run_task (model_file, priors_file, experiments_file, \
        iterations_phase1, sigma_update_n, iterations_phase2, \
        iterations_phase3, nof_process):
    import os
    os.system ("export PYTHONPATH=\"${PYTHONPATH}\":/home/gestrela/cs/SigNetMS")
    print ("PAAATH!!" + sys.path)
    sbml = SBML ()
    sbml.load_file (model_file)
    odes = sbml_to_odes (sbml)
    experiments = ExperimentSet (experiments_file)
    theta_priors = define_sbml_params_priors (sbml, priors_file)
    ml = MarginalLikelihood (int (iterations_phase1), 
                             int (sigma_update_n), 
                             int (iterations_phase2), 
                             int (iterations_phase3), 20, 2, \
                             verbose=False, n_process=int (nof_process))
    log_l = ml.estimate_marginal_likelihood (experiments, odes, 
            theta_priors)
    return log_l


parser = argparse.ArgumentParser ()
parser.add_argument ("tasks_file", help="A JSON file that defines" \
        + " the tasks that should be performed.")
parser.add_argument ("cluster_definition_file", help="A JSON file" \
        + " that defines the cluster that should be used to perform" \
        + " tasks")
args = parser.parse_args ()
tasks_filename = args.tasks_file
cluster_filename = args.cluster_definition_file
tasks_file = open (tasks_filename, 'r')
cluster_file = open (cluster_filename, 'r')

tasks_json = json.load (tasks_file)
cluster_json = json.load (cluster_file)
workers = cluster_json["machines"]
ray_path = cluster_json["ray_path"]

# Prepares Ray on all machines
redis_address = prepare_workers (workers, ray_path)

# Initializes Ray for this program
ray.init (redis_address=redis_address)

sent_tasks = {}
for task in tasks_json:
    task_id = run_task.remote (*tasks_json[task])
    sent_tasks[task] = task_id

resulting_ml = {}
for task in tasks_json:
    resulting_ml[task] = ray.get (sent_tasks[task])
