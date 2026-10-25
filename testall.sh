#!/bin/sh

COMMAND="timeout 180s ./main.out"

for var in "$@"
do
	echo "startgraph $var"

	for mode in karger ilp bnp bnp_nohot bnp_hp
	do
		echo "start $mode"
		$COMMAND "$var" "$mode"
		echo "exit_code $?"
		echo "end $mode"
	done

	echo "endgraph $var"
done
