import matplotlib.pyplot as plt
import numpy as np
from scipy.stats.distributions import chi2 as chi_squared
from scipy.stats import chisquare as chi2test
import scipy.stats

def model(x,k,l,i):
    return x - k*(l**i)

def chi2_old(x,k,l, xs, stdevs):
    return sum([((xs[i] - model(x,k,l,i))**2)/(stdevs[i]**2) for i in range(len(xs))])

def chi2(x,k,l, xs, stdevs):
    return sum([((xs[i] - model(x,k,l,i))**2)/(model(x,k,l,i)) for i in range(len(xs))])





# mt_best_fit = [0.6669888524756883, 0.04368502583636777, 0.538353712593029]
# xoshiro_best_fit = [0.666994594044332, 0.04173434349280304, 0.5460478925735742]

wagner_xs = [.6402, .6543, .6613, .6618, .6653, .6672, .6681, .6687]
wagner_stdvs = [.0024, .0015, .0006, .0003, .0002, .0001, .0001, .00004]
wagners_best_fit = [0.6700251595116361, 0.025248162255475347, 0.6544557120706123]
wagners_paper_fit = [0.6698, 0.0269, 0.63]
wagner_data = (wagner_xs, wagner_stdvs, wagners_best_fit[0],wagners_best_fit[1], wagners_best_fit[2], 5, "Wagner")
wagner_paper_data =  (wagner_xs, wagner_stdvs, wagners_paper_fit[0],wagners_paper_fit[1], wagners_paper_fit[2], 5, "Wagner Paper")

mt_xs = [0.60653, 0.63885, 0.65309, 0.66033, 0.66347, 0.66512, 0.66598, 0.666371, 0.666681]
mt_stdvs = [0.00916, 0.00265, 0.00144, 0.00074, 0.00025, 0.00013, 0.00005, 0.000039, 0.000032]
mt_best_fit = [0.6668783085496017, 0.054295045308503126, 0.5048466416984654]
mt_best_fit_no_25600 = [0.6667735182349045, 0.057655461658019, 0.48775143007859784]
mt_data = (mt_xs, mt_stdvs, mt_best_fit[0],mt_best_fit[1], mt_best_fit[2], 6, "MT")
mt_data_no_25600 = (mt_xs[:-1], mt_stdvs[:-1], mt_best_fit_no_25600[0],mt_best_fit_no_25600[1], mt_best_fit_no_25600[2], 5, "MT no 25600")


xoshiro_xs = [0.60193, 0.63981, 0.6536, 0.66031, 0.66349, 0.66502, 0.66596, 0.666354, 0.666686]
xoshiro_stdev = [0.01197,0.00332,0.00128,0.00033,0.00024,0.00013,0.00006, 0.000024, 0.000024]
xoshiro_best_fit = [0.6669349726043318, 0.04642388913080165, 0.5287992915131047]
xoshiro_best_fit_no_25600 = [0.6667557806714918, 0.05651638187681013, 0.49134860235577676]
xoshiro_data = (xoshiro_xs, xoshiro_stdev, xoshiro_best_fit[0],xoshiro_best_fit[1], xoshiro_best_fit[2], 6, "XOSHIRO")
xoshiro_data_no_25600 = (xoshiro_xs[:-1], xoshiro_xs[:-1], xoshiro_best_fit_no_25600[0],xoshiro_best_fit_no_25600[1], xoshiro_best_fit_no_25600[2], 5, "XOSHIRO no 25600")

xoshiro_xs_omit_100 = [0.63981, 0.6536, 0.66031, 0.66349, 0.66502, 0.66596, 0.666354, 0.666686]
xoshiro_stdev_omit_100 = [0.00332,0.00128,0.00033,0.00024,0.00013,0.00006, 0.000024, 0.000024]
xoshiro_best_fit_omit_100= [0.6669589383898806, 0.023570667997701633, 0.5368165492586843]
xoshiro_omit_100_data = (xoshiro_xs_omit_100, xoshiro_stdev_omit_100, xoshiro_best_fit_omit_100[0], xoshiro_best_fit_omit_100[1], xoshiro_best_fit_omit_100[2], 5, "XOSHIRO no 100")


def chi2stats(data):
    xs, stdevs, x, k, l, df, name = data

    sample_chi2 = chi2_old(x,k,l,xs,stdevs)
    print(name)
    print(sample_chi2)  
    print(chi_squared.isf(.05, df))
    print(1 - chi_squared.cdf(sample_chi2, df))

    plt.plot(xs, "o")
    plt.plot([model(x,k,l,i) for i in range(len(xs))])
    plt.title(name)
    plt.show()

# chi2stats(wagner_paper_data)
# chi2stats(mt_data)
chi2stats(xoshiro_data)
# chi2stats(mt_data_no_25600)
# chi2stats(xoshiro_data_no_25600)
chi2stats(xoshiro_omit_100_data)


# given sample_chi2 value, the following gives the approximate percent chance 
# of a new test_chi2 being larger than sample_chi2
# (i.e. being more different from expected distribution) 
# assuming the hyp is correct

# The p-value is the probability that a chi-square statistic having n degrees of freedom is more extreme than sample_chi2
# i.e. the integral of the chisquare dist from sample_chi2 to inf

# print(1 - chi_squared.cdf(sample_chi2, 6))