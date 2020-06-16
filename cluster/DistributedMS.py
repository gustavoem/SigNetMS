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
    (out, _) = proc.communicate ()

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
        (out, _) = proc.communicate ()
        worker_up_pattern = re.compile (".*Started Ray on this node")
        
        print (out.decode ("utf-8"))
        if worker_up_pattern.search (str (out)):
            print ("[OK]")
        else:
            print ("Failed!")
            warnings.warn ("Could not start ray on " + worker + ".")
    return redis_address


def stop_clusters (workers, ray_path):
    current_path = str (Path ().absolute ())
    stop_worker_bin = current_path + "/stop_cluster_worker.sh"
    for worker in workers[::-1]:
        print ("Stoping worker ", worker + "...", end="", flush=True)
        proc = subprocess.Popen ([stop_worker_bin, worker, ray_path], \
                stdout=subprocess.PIPE)
        (out, _) = proc.communicate ()
        return_code = proc.returncode
        if return_code:
            print (out.decode ("utf-8"))
            warnings.warn ("Failed to stop ray!")
        else:
            print ("[OK]")


def abs_path (path, prefix_path):
    return prefix_path + "/" + path


def get_signetms_path (current_path):
    path_list = current_path.split ('/')
    signetms_idx = path_list.index ('SigNetMS')
    signetms_path_list = path_list[:signetms_idx + 1]
    signetms_path = '/'.join (signetms_path_list)
    return signetms_path


@ray.remote
def run_task (model_file, priors_file, experiment_file, \
        iterations_phase1, sigma_update_n, iterations_phase2, \
        iterations_phase3, nof_process, signetms_path, seed=42):
    import time
    start_time = time.time()
    # importing local modules...
    sys.path.insert (0, signetms_path)
    from model.SBML import SBML
    from experiment.ExperimentSet import ExperimentSet
    from model.SBMLtoODES import sbml_to_odes
    from marginal_likelihood.MarginalLikelihood \
            import MarginalLikelihood
    from model.PriorsReader import define_sbml_params_priors
    import seed_manager

    # Now the actual code...
    seed_manager.set_seed(seed)
    sbml = SBML ()
    sbml.load_file (model_file)
    odes = sbml_to_odes (sbml)
    experiments = ExperimentSet (experiment_file)
    theta_priors = define_sbml_params_priors (sbml, priors_file)
    ml = MarginalLikelihood (iterations_phase1, 
                             sigma_update_n, 
                             iterations_phase2, 
                             iterations_phase3, 20, 2, \
                             verbose=False, n_process=nof_process)
    log_l = ml.estimate_marginal_likelihood (experiments, odes, 
            theta_priors)
    elapsed_time = time.time () - start_time
    return log_l, elapsed_time


parser = argparse.ArgumentParser ()
parser.add_argument ("tasks_file", help="A JSON file that defines" \
        + " the tasks that should be performed.")
parser.add_argument ("cluster_definition_file", help="A JSON file" \
        + " that defines the cluster that should be used to perform" \
        + " tasks")
parser.add_argument ("--seed", type=int, nargs="?", default=0, \
        help="Random number generation seed")
args = parser.parse_args ()
tasks_filename = args.tasks_file
cluster_filename = args.cluster_definition_file
arg_seed = args.seed

tasks_file = open (tasks_filename, 'r')
cluster_file = open (cluster_filename, 'r')

tasks_json = json.load (tasks_file)
cluster_json = json.load (cluster_file)
workers_list = cluster_json["machines"]
ray_abs_path = cluster_json["ray_path"]
process_by_task = cluster_json["n_process_by_task"]

# Prepares Ray on all machines
redis_server_address = prepare_workers (workers_list, ray_abs_path)

# Initializes Ray for this program
ray.init (redis_address=redis_server_address)


# Find SigNetMS and input absolute path
current_abs_path = str (Path ().absolute ())
signetms_abs_path = get_signetms_path (current_abs_path)

# Send tasks
not_ready_tasks = []
id_to_name = {}
for task in tasks_json:
    model_abs_path = abs_path (task["model_file"], signetms_abs_path)
    prior_abs_path = abs_path (task["prior_file"], signetms_abs_path)
    experiments_abs_path = abs_path (task["experiment_file"], \
            signetms_abs_path)

    task_id = run_task.remote (
            model_abs_path, 
            prior_abs_path,
            experiments_abs_path,
            int (task["phase1_it"]),
            int (task["sigma_update_n"]),
            int (task["phase2_it"]),
            int (task["phase3_it"]),
            int (process_by_task),
            signetms_abs_path,
            seed=arg_seed)
    id_to_name[str (task_id)] = task["name"]
    not_ready_tasks.append (task_id)
    print ("Creating task", task["name"])

results = {}
while len (not_ready_tasks) > 0:
    ready, not_ready = ray.wait (not_ready_tasks)
    first_ready = ready[0]
    task_name = id_to_name[str (first_ready)]
    print ("Getting task", task_name)
    results[task_name] = ray.get (first_ready)
    not_ready_tasks = not_ready

print (results)
output_file = open("cluster_results.txt", "w")
for model in results:
    output_file.write(model + ": ")
    output_file.write(str(results[model]) + '\n')
output_file.close()

ray.shutdown ()
stop_clusters (workers_list, ray_abs_path)
