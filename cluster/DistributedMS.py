import argparse
import json
import subprocess 
import re
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

    return

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

prepare_workers (workers, ray_path)
