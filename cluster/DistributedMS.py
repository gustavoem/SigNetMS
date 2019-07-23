import argparse
import json
import subprocess 
import re
from pathlib import Path

def prepare_workers (workers, ray_path):
    # start the first worker as a host
    current_path = str (Path ().absolute ())
    start_cluster_head_bin = current_path + "/start_cluster_head.sh"
    proc = subprocess.Popen ([start_cluster_head_bin, workers[0], \
            ray_path], stdout=subprocess.PIPE)
    (out, err) = proc.communicate ()
    # print (out.decode ("utf-8"))
    redis_address_pattern = \
            re.compile (r".*redis-address (\d+\.\d+\.\d+\.\d+:\d+)")
    m = redis_address_pattern.match (str (out))
    redis_address = m.group (1)


    
    
    for worker in workers[1:]:
        # start other worker as a normal worker
        pass


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
