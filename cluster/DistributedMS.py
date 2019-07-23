import argparse
import json
import os

def prepare_workers (workers, ray_path):
    # start the first element as a host
    os.sys ("start_cluster_head.sh " + workers[0] + ray_path)


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
