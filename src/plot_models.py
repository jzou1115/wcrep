import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm

def abline(slope, intercept, minS, maxS):
    x_vals = np.array([minS, maxS])
    y_vals = intercept + slope * x_vals
    return x_vals, y_vals

def plotWC(ax, slope, sd, minS, maxS):
    mean_x, mean_y = abline(slope, 0 , minS, maxS)
    p = ax.plot(mean_x, mean_y, color="#0571b0", label="WC")
    
    upper_x, upper_y = abline(slope, 2*sd, minS, maxS)
    ax.plot(upper_x, upper_y, linestyle="dashed", color="#0571b0")
    
    lower_x, lower_y = abline(slope, -2*sd, minS, maxS)
    ax.plot(lower_x, lower_y,linestyle="dashed", color="#0571b0")
    return p

def plotWCC(ax, slope, sd, minS, maxS):
    mean_x, mean_y = abline(slope, 0 , minS, maxS)
    p = ax.plot(mean_x, mean_y, color="#ca0020", label="WC+C")
    
    upper_x, upper_y = abline(slope, 2*sd, minS, maxS)
    ax.plot(upper_x, upper_y, linestyle="dashed", color="#ca0020")
    
    lower_x, lower_y = abline(slope, -2*sd, minS, maxS)
    ax.plot(lower_x, lower_y,linestyle="dashed", color="#ca0020")
    return p

def plotConditionalDistributions(ax, data, var_g, var_c1, var_c2, threshold, title, legend):
    
    n1 = data.iloc[0, 2]
    n2 = data.iloc[0,3]
    
    z = norm.ppf(threshold/2)
    sig = data[abs(data["s1"]) > abs(z)]
    
    z2 = norm.ppf(0.05/(2*sig.shape[0]))
    replicated = []
    not_replicated = []
    for i in sig.index:
        if abs(sig.loc[i, "s2"]) > abs(z2):
            #check if sign is the same:
            if np.sign(sig.loc[i, "s2"]) == np.sign(sig.loc[i, "s1"]):
                replicated.append(i)
            else:
                not_replicated.append(i)
        else:
            not_replicated.append(i)
        

    slope = (math.sqrt(n1*n2)*var_g)/(1+n1*var_g+n1*var_c1)
    sd = math.sqrt(n2*var_g+n2*var_c2+1-((n1*n2*var_g*var_g)/(n1*var_g+n1*var_c1+1)))
    slope_noConfounding = (math.sqrt(n1*n2)*var_g)/(1+n1*var_g)
    sd_noConfounding = math.sqrt(n2*var_g+1-((n1*n2*var_g*var_g)/(n1*var_g+1)))

    
    maxS = 1.2*max(max(abs(sig["s1"])), max(abs(sig["s2"])))
    minS = -1.2*maxS

    #plot summary statistics
    p = ax.scatter(sig.loc[replicated, "s1"], sig.loc[replicated, "s2"], s=10, c="k", label="Replicated")
    ax.scatter(sig.loc[not_replicated, "s1"], sig.loc[not_replicated, "s2"], s=10, c="grey", label="Not Replicated")
    
    #plot WC
    plotWC(ax, slope_noConfounding, sd_noConfounding, minS, maxS)
    
    #plot WC+C
    plotWCC(ax, slope, sd, minS, maxS)
    
    #axes
    ax.set_xlim([minS,maxS])
    ax.set_ylim([minS,maxS])
    ax.set_title(title)
    ax.set_xlabel("Discovery z-scores")
    ax.set_ylabel("Replication z-scores")
    

    return ax


def plotPheno(modelF, dataF, n1, n2, title, ax, legend=False):
    model = open(modelF, "r")
    lines = model.readlines()
    var_g = float(lines[-3])
    var_c1 = float(lines[-2])
    var_c2 = float(lines[-1])
    model.close()

    data = pd.read_table(dataF, sep="\t")

    data["n1"] = n1
    data["n2"] = n2

    threshold = 0.05/data.shape[0]

    data.columns = ["s1", "s2", "n1", "n2"]

    return plotConditionalDistributions(ax, data, var_g, var_c1, var_c2, threshold, title, legend)

def main():
    #parse args
    args = sys.argv
    modelF = args[1]
    dataF = args[2]
    n1 = int(args[3])
    n2 = int(args[4])
    outfile = args[5]
    title = args[6]

    fig, axes = plt.subplots(1,2)

    #main plot
    ax = axes[0]
    n=0
    plotPheno(modelF, dataF, n1, n2, title, ax)
    handles, labels = ax.get_legend_handles_labels()

    #Legend
    ax = axes[1]
    ax.axis('off')
    ax.legend(handles, labels, loc="center left")


    fig.set_size_inches(6, 4, forward = True)
    fig.tight_layout()
    plt.savefig(outfile)
    plt.clf()
    plt.close()

if __name__ == "__main__":
    main()
