n=70
c=6
pv=0.6
pe=0.7
#px=0.03

BIN="../../g2013.out"

for px in $(seq 0.03 0.01 0.5);
do
	$BIN $c $n $pv $pe $px > "g2013c${c}n${n}pv${pv}pe${pe}px${px}.graph"
done

