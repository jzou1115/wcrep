# import libraries
import os
import math
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import minimize
from scipy.stats import norm
from scipy.stats import multivariate_normal


def get_near_psd(A):
    C = (A + A.T)/2
    eigval, eigvec = np.linalg.eig(C)
    eigval[eigval < 0] = 0

    return eigvec.dot(np.diag(eigval)).dot(eigvec.T)

def MLE(stats, n1, n2):
    
    def MLE_joint(params):
        sigma_g = params[0]
        sigma_c1 = params[1]
        sigma_c2 = params[2]
        
        mean = [0, 0]
        sigma_1 = n1*sigma_g*sigma_g + n1*sigma_c1*sigma_c1 + 1
        sigma_2 = n2*sigma_g*sigma_g + n2*sigma_c2*sigma_c2 + 1
        cov_s1_s2 = float(math.sqrt(n1)*math.sqrt(n2)*sigma_g*sigma_g)
        cov = np.array([[sigma_1, cov_s1_s2], [cov_s1_s2, sigma_2]])
        try:
            ll = -np.sum(multivariate_normal.logpdf(stats, mean, cov, allow_singular=True))
        except:
            cov = get_near_psd(cov) #estimation error, force covariance to be PSD matrix
            ll = -np.sum(multivariate_normal.logpdf(stats, mean, cov, allow_singular=True))
            
        #print(ll)
        return ll
    results = minimize(MLE_joint, [.1,.1,.1] , method = 'Nelder-Mead', options={'disp': False})
    
    sigma_g = results["x"][0]
    sigma_c1 = results["x"][1] 
    sigma_c2 = results["x"][2] 

    if sigma_g < 0:
        sigma_g = abs(sigma_g)
	
    if sigma_c1 < 0:
        sigma_c1 =abs(sigma_g)

    if sigma_c2 < 0:
        sigma_c2 = abs(sigma_g)
    
    return sigma_g, sigma_c1, sigma_c2





def estimate_parameters(s1, s2, n1, n2, z1):
    
    #estimate parameters
    stats = pd.concat([s1, s2], axis = 1)
    sigma_g, sigma_c1, sigma_c2 = MLE(stats, n1, n2)
    
    
    return sigma_g, sigma_c1, sigma_c2

def expected_rep_wc(d, parameters, n1, n2, z2):
    sigma_g = parameters[0]

    sigma_1 = n1*sigma_g*sigma_g  + 1
    sigma_2 = n2*sigma_g*sigma_g + 1
    cov_s1_s2 = float(math.sqrt(n1)*math.sqrt(n2)*sigma_g*sigma_g)
    
    d_pos = d[d["s1"]>0]
    mean = cov_s1_s2/sigma_1 *   d_pos["s1"]
    var =  sigma_2 - (cov_s1_s2*cov_s1_s2)/(sigma_1)
    prob_pos = 1 - norm.cdf(-z2, loc = mean, scale = np.sqrt(var) )

    d_neg = d[d["s1"]<=0]
    mean = cov_s1_s2/sigma_1 *   d_neg["s1"]
    prob_neg = norm.cdf(z2, loc = mean, scale = np.sqrt(var) )
    
    prob = np.sum(prob_pos)+  np.sum(prob_neg)
    prob = prob / d.shape[0]

    return prob  
    

def expected_rep_wcc(d, parameters, n1, n2, z2):
    sigma_g = parameters[0]
    sigma_c1 = parameters[1]
    sigma_c2 = parameters[2]
    
    sigma_1 = n1*sigma_g*sigma_g + n1*sigma_c1*sigma_c1 + 1
    sigma_2 = n2*sigma_g*sigma_g + n2*sigma_c2*sigma_c2 + 1
    cov_s1_s2 = float(math.sqrt(n1)*math.sqrt(n2)*sigma_g*sigma_g)
    
    d_pos = d[d["s1"]>0]
    mean = cov_s1_s2/sigma_1 *   d_pos["s1"]
    var =  sigma_2 - (cov_s1_s2*cov_s1_s2)/(sigma_1)
    prob_pos = 1 - norm.cdf(-z2, loc = mean, scale = np.sqrt(var) )
    
    d_neg = d[d["s1"]<=0]
    mean = cov_s1_s2/sigma_1 *   d_neg["s1"]
    prob_neg = norm.cdf(z2, loc = mean, scale = np.sqrt(var) )
    
    
    prob = np.sum(prob_pos)+  np.sum(prob_neg)
    prob = prob / d.shape[0]    
    return prob 

        
def fit_models(sig, n1, n2, threshold1=5e-8):

    z1 = norm.ppf(threshold1/2)

    if sig.shape[0]>0:
        threshold2 = 0.05/sig.shape[0]
    else:
        threshold2 = 0.05
    z2 = norm.ppf(threshold2/2)
    
    rep = sig[abs(sig["s2"])> abs(z2)]
    
    rep = rep[np.sign(rep["s1"])*np.sign(rep["s2"])==1]
   
    if sig.shape[0]>0: 
	    rep_rate_true = rep.shape[0]/sig.shape[0]
    else:
        rep_rate_true = 0
    
    parameters = estimate_parameters(sig["s1"], sig["s2"], n1, n2, z1)
    
    rep_rate_wc = expected_rep_wc(sig, parameters, n1, n2, z2)
    
    rep_rate_wcc = expected_rep_wcc(sig, parameters, n1, n2, z2)
    
    replication = [rep_rate_true, rep_rate_wc, rep_rate_wcc]
    
    return(replication, parameters, sig.shape[0])

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

def plotConditionalDistributions(ax, sig, n1, n2, var_g, var_c1, var_c2, threshold, title=""):
	 
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

def main():
	#parse command line input
	args = sys.argv
	dataF = args[1]
	prefix = args[2]
	print(dataF, prefix)

	data = pd.read_table(dataF)
	data.columns = ["s1", "s2", "n1", "n2", "t"]
	n1 = np.median(data["n1"])
	n2 = np.median(data["n2"])
	t = data["t"][0]

	replication, parameters, numSig = fit_models(data, n1, n2, t)

	#write output
	out_file = open(prefix+ ".txt", "w")
	out_file.write("file\tn1\tn2\tt\trep_true\trep_wc\trep_wcc\tsigma_g\tsigma_c1\tsigma_c2\tnumSig\n")
	out_file.write("%s\t%d\t%d\t%.10f\t%f\t%f\t%f\t%.10f\t%.10f\t%.10f\t%d\n" % (dataF, n1, n2, t, replication[0], replication[1], replication[2], parameters[0], parameters[1], parameters[2], numSig))
	out_file.close()

	#make plot
	out_plot = prefix+ ".png"
	fig, axes = plt.subplots(1,2)

    #main plot
	ax = axes[0]
	plotConditionalDistributions(ax, data, n1, n2, np.power(parameters[0], 2), np.power(parameters[1], 2), np.power(parameters[2], 2), t)
	handles, labels = ax.get_legend_handles_labels()

	#Legend
	ax = axes[1]
	ax.axis('off')
	ax.legend(handles, labels, loc="center left")


	fig.set_size_inches(6, 4, forward = True)
	fig.tight_layout()
	plt.savefig(out_plot)
	plt.clf()
	plt.close()

if __name__ == "__main__":
	main()
