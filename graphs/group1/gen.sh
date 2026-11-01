n=100
c=6
d=0.1

BIN="../../generator.out"

for seed in $(seq 1 1 20);
do
	$BIN $seed $n $c $d > "gen1seed${seed}n${n}c${c}d${d}.graph"
done

