
import numpy as np
import matplotlib.pyplot as plt
import sys
input_file=sys.argv[1]

def read_params(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

        i = 0
        for line in lines:
            if "Number of steps:" in line:
                steps = int(lines[i +1])

            if "Time evolution:" in line:
                file1 = lines[i + 1]
                file1_split = file1.split('"')[1]
                file2 = lines[i + 3]
                file2_split = file2.split('"')[1]
                break
            i = i + 1

    return steps,file1_split,file2_split

steps, file1_split, file2_split= read_params(input_file)

def read_values(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    time = []
    kin = []
    pot = []
    tot = []
    temp = []
    press = []

    for line in lines:
        if "#" in line:
            continue
        
        line_split = line.split()

        time.append(float(line_split[0]))
        kin.append(float(line_split[1]))
        pot.append(float(line_split[2]))
        tot.append(float(line_split[3]))
        temp.append(float(line_split[4]))
        press.append(float(line_split[5]))

    
    return time, kin, pot, tot, temp, press

time, kin, pot, tot, temp, press = read_values(file1_split)

def read_msd(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    time = []
    sqr_dist = []

    for line in lines:
        if "#" in line:
            continue
        
        line_split = line.split()

        time.append(float(line_split[0]))
        sqr_dist.append(float(line_split[1]))

    return time, sqr_dist

time_dist, sqr_dist = read_msd(file2_split)


plt.xlim(0,time[-1])
plt.xlabel("Time (ps)")
plt.ylabel("Energy (kJ/mol")
plt.plot(time,kin, label="Kinetic Energy")
plt.plot(time,pot, label="Potential Energy")
plt.plot(time,tot, label="Total Energy")
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12), ncol=3)
plt.savefig("energy_fig.png", bbox_inches="tight", dpi=300)
plt.clf()

plt.xlim(0,time[-1])
plt.plot(time,press, label="Pressure")
plt.ylabel("Pressure (Pa)")
plt.xlabel("Time (ps)")
plt.savefig("press_fig.png", bbox_inches="tight", dpi=300)
plt.clf()

plt.xlim(0,time[-1])
plt.plot(time,temp, label="Temperature")
plt.ylabel("Temperature (K)")
plt.xlabel("Time (ps)")
plt.savefig("temp_fig.png", bbox_inches="tight", dpi=300)
plt.clf()

plt.xlim(0,time_dist[-1])
plt.plot(time_dist,sqr_dist, label="Mean Square Distance")
plt.ylabel("Mean Square Distance (${\AA}^{2}$)")
plt.xlabel("Time (ps)")
plt.savefig("msd.png", bbox_inches="tight", dpi=300)
plt.clf()

