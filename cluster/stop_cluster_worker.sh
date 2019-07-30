#!/bin/sh

worker=$1
ray_path=$2
ray_out=$(ssh $worker $ray_path stop 2>&1 > /dev/null)
return_code=$?
echo "$ray_out"
return $return_code
