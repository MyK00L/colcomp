import statistics
import matplotlib.pyplot as plt
import numpy as np

def nsolved(data, mode):
	return len(list(filter(lambda item: (item[1][mode]["exit_code"] == 0), data.items())))

def printstats(name, l):
	print(name, "mean:", statistics.mean(l))
	#print(name, "median:", statistics.median(l))
	#print(name, "min:", min(l))
	#print(name, "max:", max(l))

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
	byall = list(filter(lambda item: (item[1]["bnp"]["exit_code"] == 0), data.items() ))
	printstats_multi(["stats npp", "stats"], ["N", "M", "maxN", "maxM", "ccs"], data.items())
	printstats_multi(["bnp"], ["time_millis", "bnp_nodes", "pricing_bnb_nodes"], byall)

def process_data2024(data):
	smallish = list(filter(lambda item: (item[1]["stats npp"]["N"] >= 10 and item[1]["stats npp"]["N"] <= 210), data.items() ))
	print(len(smallish))


def process_datag(data):
	print("total: ", len(data))
	print("improved by preprocess: ", len(list(filter(lambda item: (item[1]["stats npp"]["M"] != item[1]["stats"]["M"]), data.items() ))))
	printstats_multi(["bnp", "bnp nhp", "bnp npp", "bnp nhs"], ["time_millis", "bnp_nodes", "pricing_steps", "pricing_bnb_runs", "pricing_bnb_nodes"], data.items())


def process_data_by_fnp(ax, data, fnp, stat, fun):
	print("total: ", len(data))
	print("bnp solved: ", nsolved(data,"bnp"))
	byall = list(filter(lambda item: (item[1]["bnp"]["exit_code"] == 0), data.items() ))
	printstats_multi(["stats npp", "stats"], ["N", "M", "maxN", "maxM", "ccs"], data.items())
	printstats_multi(["bnp"], ["time_millis", "bnp_nodes", "pricing_bnb_nodes"], byall)
	stats = list(item[1] for item in data.items())
	stats = sorted(stats, key=lambda x: float(x["fnp"][fnp]))
	x = list(float(item["fnp"][fnp]) for item in stats)
	y = list(fun(item["bnp"]) for item in stats)
	ax.grid()
	ax.set_xlabel(fnp)
	ax.set_ylabel("time in seconds")
	ax.scatter(x,y)
	ax.set_yscale("log")
	ax.set_xlim(left=0)

def process_data_bys():
	fig, ax = plt.subplots(2,3)
	py=0
	px=0
	for stat in ["n","c","px","pv","pe"]:
		data = read_data("data_by_"+stat+".txt")
		process_data_by_fnp(ax[py][px], data, stat, "time_millis", lambda x: float(x["time_millis"] if "time_millis" in x else 0)/1000)
		py+=1
		if py>=2:
			py=0
			px+=1
	fig.delaxes(ax[1][2])
	# plt.savefig('figure.pdf', format='pdf', bbox_inches='tight')
	plt.show()

def get_filename_params(name):
	res = {}
	key = ""
	i=0
	while(i<len(name)):
		i0=i
		while i<len(name) and name[i].isalpha():
			i+=1
		key=name[i0:i]
		i0=i
		while i<len(name) and not name[i].isalpha():
			i+=1
			try:
				float(name[i0:i])
				res[key]=name[i0:i]
			except:
				pass
	return res


def read_data(filename):
	file = open(filename, "r")
	graph = ""
	mode = ""
	data = {}
	for line in file:
		line = line.strip()
		if line.startswith("start ") or line.startswith("end "):
			sp = line.split(" ", 1)
		else:
			sp = line.split()
		if sp[0]=="startgraph":
			graph=sp[1]
			if not graph in data:
				data[graph] = {}
			data[graph]["fnp"] = get_filename_params(sp[1])
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
			data[graph][mode][sp[0]]=float(sp[1])
	return data

def main():
	data = read_data("data_group1.txt")
	process_datag(data)
	#process_data_bys()

if __name__ == "__main__":
	plt.rcParams['font.size'] = 12
	plt.rcParams['axes.labelsize'] = 14
	plt.rcParams['axes.titlesize'] = 16
	plt.rcParams['xtick.labelsize'] = 11
	plt.rcParams['ytick.labelsize'] = 11
	plt.rcParams['legend.fontsize'] = 11
	plt.rcParams['figure.dpi'] = 100
	plt.rcParams['savefig.dpi'] = 300
	#plt.rcParams['text.usetex'] = True
	plt.rcParams['font.family'] = 'serif'
	plt.rcParams['lines.linewidth'] = 2
	plt.rcParams['lines.markersize'] = 6
	main()