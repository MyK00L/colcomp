import statistics
import matplotlib.pyplot as plt
import numpy as np

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


def process_dataheur(data):
	print("total: ", len(data))
	printstats_multi(["karger", "merge"], ["time_millis","score"], data.items())

def process_data2012(data):
	print("total: ", len(data))
	print("bnp solved: ", nsolved(data,"bnp"))
	print("bnp_hp solved: ", nsolved(data,"bnp_hp"))
	print("bnp_nohot solved: ", nsolved(data,"bnp_nohot"))
	byall = list(filter(lambda item: (item[1]["bnp"]["exit_code"] == "0" and item[1]["bnp_hp"]["exit_code"] == "0" and item[1]["bnp_nohot"]["exit_code"] == "0"), data.items() ))
	printstats_multi(["bnp", "bnp_hp", "bnp_nohot"], ["time_millis", "bnp_nodes", "pricing_bnb_nodes"], byall)

def process_data2024(data):
	print("total: ", len(data))
	print("solved: ", nsolved(data,"bnp_hp"))
	printstats_multi(["bnp_hp"], ["time_millis", "bnp_nodes", "pricing_bnb_nodes"], list(filter(lambda item: (item[1]["bnp_hp"]["exit_code"] == "0"),data.items()) ))
	x = [0]+sorted(list(map(lambda item: int(item[1]["bnp_hp"]["time_millis"])/1000, filter(lambda item: (item[1]["bnp_hp"]["exit_code"] == "0"),data.items()))))
	y = range(len(x))
	plt.semilogx(x,y)
	plt.grid(axis='y')
	plt.xlabel("time in seconds")
	plt.ylabel("instances solved")
	plt.show()

	# TODO: charts

def process_datag1(data):
	print("total: ", len(data))
	for d in ["c\\d",0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
		print(d,end=" & ")
	print("\\\\")
	for c in [3,6,9,12,15,18,21,24,27]:
		print(c,end=" & ")
		for d in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
			pbnsum=0
			for seed in [1,2,3,4,5]:
				name = "./graphs/group1/s"+str(seed)+"n30c"+str(c)+"d"+str(d)+".graph"
				pbnsum += int(data[name]["bnp_hp"]["pricing_bnb_nodes"])
			print(int(pbnsum/5),end=" & ")
		print("\\\\")

def process_datag2(data):
	c = 6
	d = 0.125
	x=[]
	y=[]
	for n in [10,20,30,40,50,60,70,80,90]:
		pbnsum=0
		pbnnum=0
		for seed in [1,2,3,4,5]:
			name = "./graphs/group2/s"+str(seed)+"n"+str(n)+"c"+str(c)+"d"+str(d)+".graph"
			if data[name]["bnp_hp"]["exit_code"]=="0":
				pbnsum += int(data[name]["bnp_hp"]["pricing_bnb_nodes"])
				pbnnum += 1	
		print(int(pbnsum/pbnnum))
		x.append(n)
		y.append(pbnsum/pbnnum)

	plt.semilogy(x,y)
	plt.show()



def read_data(filename):
	file = open(filename, "r")
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
	return data

def main():
	data = read_data("dataheur.txt")
	process_dataheur(data)

if __name__ == "__main__":
	main()