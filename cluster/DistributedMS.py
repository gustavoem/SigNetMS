import sys
sys.path.insert(0, '..')

import argparse
import json
import subprocess 
import re
import ray
import warnings
from pathlib import Path

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

    print (out.decode ("utf-8"))
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
        
        print (out.decode ("utf-8"))

        if worker_up_pattern.search (str (out)):
            print ("[OK]")
        else:
            print ("Failed!")
            warnings.warn ("Could not start ray on " + worker + ".")
            if err is not None:
                print (err.decode ("utf-8"))
    return redis_address


def stop_clusters (workers, ray_path):
    current_path = str (Path ().absolute ())
    stop_worker_bin = current_path + "/stop_cluster_worker.sh"
    
    for worker in workers[::-1]:
        print ("Stoping worker ", worker, "...", end="", flush=True)
        proc = subprocess.Popen ([stop_worker_bin, worker, ray_path], \
                stdout=subprocess.PIPE)
        (out, err) = proc.communicate ()
        return_code = proc.returncode
        if return_code:
            print (out.decode ("utf-8"))
            warnings.warn ("Failed to stop ray!")
        else:
            print ("[OK]")


def abs_path (path, SIGNET_MS_PATH):
    return SIGNET_MS_PATH + "/" + path


@ray.remote
def run_task (model_file, priors_file, experiments_file, \
        iterations_phase1, sigma_update_n, iterations_phase2, \
        iterations_phase3, nof_process, SIGNET_MS_PATH):
    return 0
    # importing local modules...
    sys.path.insert (0, SIGNET_MS_PATH)
    from model.SBML import SBML
    from experiment.ExperimentSet import ExperimentSet
    from model.SBMLtoODES import sbml_to_odes
    from marginal_likelihood.MarginalLikelihood \
            import MarginalLikelihood
    from model.PriorsReader import define_sbml_params_priors
    from model.ODES import ODES

    # Now the actual code...
    model_abs_path = abs_path (model_file, SIGNET_MS_PATH)
    experiments_abs_path = abs_path (experiments_file, SIGNET_MS_PATH)
    priors_abs_path = abs_path (priors_file, SIGNET_MS_PATH)
    sbml = SBML ()
    sbml.load_file (model_abs_path)
    odes = sbml_to_odes (sbml)
    experiments = ExperimentSet (experiments_abs_path)
    theta_priors = define_sbml_params_priors (sbml, priors_abs_path)
    ml = MarginalLikelihood (iterations_phase1, 
                             sigma_update_n, 
                             iterations_phase2, 
                             iterations_phase3, 20, 2, \
                             verbose=False, n_process=nof_process)
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
process_by_task = cluster_json["n_process_by_task"]

# Prepares Ray on all machines
redis_address = prepare_workers (workers, ray_path)

# Initializes Ray for this program
ray.init (redis_address=redis_address)

current_path = str (Path ().absolute ())
SIGNET_MS_PATH = '/'.join (((current_path.split ('/'))[:-1]))

sent_tasks = {}
for task in tasks_json:
    task_id = run_task.remote (
            task["model_file"], 
            task["prior_file"],
            task["experiment_file"],
            int (task["phase1_it"]),
            int (task["sigma_update_n"]),
            int (task["phase2_it"]),
            int (task["phase3_it"]),
            int (process_by_task),
            SIGNET_MS_PATH)
    task_name = task["name"]
    sent_tasks[task_name] = task_id

resulting_ml = {}
for task in tasks_json:
    task_name = task["name"]
    resulting_ml[task_name] = ray.get (sent_tasks[task_name])
print (resulting_ml)
stop_clusters (workers, ray_path)
