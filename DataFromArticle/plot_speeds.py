import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("darkgrid")
#plt.tight_layout()
ax = plt.subplot(1,2,1)

sns.set_style("darkgrid")




#acc_limited = [4.01,12.68,46.42,177.78,687.88,2724.41]
#jerk_limited = [3.37,12.49,49.94,201.04,799.08,3177.13]






#jerk_limited = [4.05, 12.51, 44.06, 162.63, 625.53, 2443.82, 9541.80]
#acc_limited = [ 2.83, 8.59, 28.89, 105.70, 391.64, 1514.82, 5925.46]
discs = [32, 64, 128, 256, 512, 1024]
jerk_limited = []
acc_limited = []
nrOfRuns = 100
for disc in discs:
	for ending in ["jerk", "acc"]:
		f = open(str(disc) + ending + ".txt", "r")
		text = f.read()
		lines = text.split("\n")
		t = float(lines[-2].split("= ")[-1].replace("milliseconds.", ""))
		print(t)
		if ending == "jerk":
			jerk_limited.append(t / nrOfRuns)
		else:
			acc_limited.append(t / nrOfRuns)
		f.close()

pol = [el*el/200 for el in [8]+discs + [2048, 2*2048]]

palette = sns.color_palette("rocket_r")

alpha = 0.8

plt.loglog(discs, jerk_limited, "-o", alpha=alpha, color=palette[5], label="Jerk Limited Dataset")
plt.loglog(discs, acc_limited, "--o", alpha=alpha, color=palette[4], label="Acceleration Limited Dataset")

plt.loglog([8]+discs + [2048, 2048*2], pol, "--", color="black", alpha=alpha, linewidth=0.8)
plt.yticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000])
#plt.xticks([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900,1000, 2000])
#plt.xticks([32, 64, 128, 256, 512, 1024])

#fig, ax = plt.subplots()
ax.set_xticks([32, 64, 128, 256, 512, 1024])
ax.set_xticklabels(["32", "64", "128", "256", "512", "1024"])

plt.xlim([29, 1100])
plt.ylim([1, 12500])
plt.legend()
plt.title("Running Time as a\nFunction of Velocity Resolution")
plt.xlabel("Number of Discrete Velocities")
plt.ylabel("Mean Running Time [ms]")

times = []
resolutions = []
import pandas as pd

#import seaborn as sns
sns.set_style("darkgrid")
plt.subplot(1,2,2)

from random import random, gauss
#import matplotlib.pyplot as plt
di = {}
datas = []
trues = {}

for ending in ["jerk", "acc"]:
	trues[ending] = []
	f = open(str(disc) + ending + ".txt", "r")
	text = f.read()
	lines = text.split("\n")
	lines.pop(-1)
	lines.pop(-1)
	f.close()
	for line in lines:
		parts = line.split(",")
		static = float(parts[0]) - float(parts[1])
		trues[ending].append(static)


for i in discs:
	print(i)

	for ending in ["jerk", "acc"]:
		print(str(i) + ending + ".txt")
		f = open(str(i) + ending + ".txt", "r")
		text = f.read()
		lines = text.split("\n")
		lines.pop(-1)
		lines.pop(-1)
		f.close()

		vec = []
		for line in lines:
			print(line)
			parts = line.split(",")
			static = float(parts[0]) - float(parts[1])
			vec.append(static)
		print(vec[0])

		for j in range(len(vec)):
			datas.append({"Lost Time [s]":vec[j], "Number of Discrete Velocities": i, "Dataset":"Acceleration Limited" if ending=="acc" else "Jerk Limited"})

df = pd.DataFrame.from_dict(datas)
print(df)
palette = sns.color_palette("rocket_r")
palette = [palette[5], palette[4]]
ax = sns.violinplot(data=df, cut=0, bw=.5, saturation=1.0, palette=palette, y="Lost Time [s]", x="Number of Discrete Velocities", hue="Dataset", split=True)
plt.setp(ax.collections, alpha=alpha)
plt.title("Lost Time as a\nFunction of Velocity Resolution")
plt.yticks([80 + i*20 for i in range(8)])
plt.show()
