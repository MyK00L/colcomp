#!/bin/sh

GEN=../../generator.out

n=30

for seed in 1 2 3 4 5
do
	for c in 3 6 9 12 15 18 21 24 27
	do
		for d in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
		do
			$GEN $seed $n $c $d > "s${seed}n${n}c${c}d${d}.graph"
		done
	done
done

