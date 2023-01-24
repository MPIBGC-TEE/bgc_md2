#!/usr/bin/env bash

start_prob=0
end_prob=1

for ((prob=start_prob; prob <= end_prob; prob++))
do
	python CARDAMOM_pwc_mr_fd_parallel_script.py $prob
done

echo "done"
