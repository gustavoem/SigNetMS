#!/bin/sh

head=$1
ray_path=$2

echo "ssh $head $ray_path start --head 2>&1 > /dev/null"
ray_out=$(ssh $head $ray_path start --head 2>&1 >/dev/null)
echo "$ray_out"

