import statistics

def nsolved(data, mode):
	return len(list(filter(lambda item: (item[1][mode]["exit_code"] == "0"), data.items())))

def printstats(name, l):
	print(name, "mean:", statistics.mean(l))
	print(name, "median:", statistics.median(l))
	print(name, "min:", min(l))
	print(name, "max:", max(l))

def printstats_multi(names, stats, data):
	for name in names:
		for stat in stats:
			printstats(name+" "+stat, list(int(item[1][name][stat]) for item in data))

def process_data1(data):
	print("total: ", len(data))
	print("bnp solved: ", nsolved(data,"bnp"))
	print("bnp_hp solved: ", nsolved(data,"bnp_hp"))
	print("bnp_nohot solved: ", nsolved(data,"bnp_nohot"))
	byall = list(filter(lambda item: (item[1]["bnp"]["exit_code"] == "0" and item[1]["bnp_hp"]["exit_code"] == "0" and item[1]["bnp_nohot"]["exit_code"] == "0"), data.items() ))
	printstats_multi(["bnp", "bnp_hp", "bnp_nohot"], ["time_millis", "bnp_nodes", "pricing_bnb_nodes"], byall)


def main():
	file = open("data1.txt", "r")
	graph = ""
	mode = ""
	data = {}
	for line in file:
		sp = line.split()
		if sp[0]=="startgraph":
			graph=sp[1]
		elif sp[0]=="start":
			mode=sp[1]
		elif sp[0]=="endgraph":
			graph=""
		elif sp[0]=="end":
			mode=""
		else:
			if not graph in data:
				data[graph] = {}
			if not mode in data[graph]:
				data[graph][mode] = {}
			data[graph][mode][sp[0]]=sp[1]
	process_data1(data)

if __name__ == "__main__":
	main()