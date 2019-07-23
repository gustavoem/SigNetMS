#!/bin/sh

worker=$1
ray_path=$2
redis_address=$3
ray_out=$(ssh $worker $ray_path start --redis-address=$redis_address 2>&1 > /dev/null)
echo "$ray_out"
