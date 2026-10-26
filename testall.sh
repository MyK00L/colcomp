#!/bin/sh

COMMAND="timeout 180s ./main.out"

for var in "$@"
do
	echo "startgraph $var"
	# bnp bnp_nohot ilp
	for mode in stats karger bnp_hp
	do
		echo "start $mode"
		$COMMAND "$var" "$mode"
		echo "exit_code $?"
		echo "end $mode"
	done

	echo "endgraph $var"
done
