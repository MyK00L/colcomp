#!/bin/sh

GEN=../../generator.out

c=6
d=0.125

for seed in 1 2 3 4 5
do
	for n in 10 20 30 40 50 60 70 80 90 100
		do
		$GEN $seed $n $c $d > "s${seed}n${n}c${c}d${d}.graph"
	done
done

