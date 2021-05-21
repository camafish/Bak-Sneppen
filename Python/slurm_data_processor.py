import matplotlib.pyplot as plt
import numpy as np
import math
import os

from sklearn.neighbors import KernelDensity
import pprint
pp = pprint.PrettyPrinter(indent=4)


def read_file(f):
    # given a reference to a slurm output file
    # return a dict with all relevant fields
    # which are 
    # 'rng', 'seed', 'n', 'k', 'time', 'relisted', 'rank_freqs', 'samples'
    data = {}
    fo = open(f.path, 'r')

    line = fo.readline()
    if (line[0] == 'i'):

        line = fo.readline()
        data['rng'] = line.split(": ")[1]
        line = fo.readline()
        data['seed'] = int(line.split(": ")[1][:-1])

        line = fo.readline()
        data['n'] = int(line.split(": ")[1][:-1])
        line = fo.readline()
        data['k'] = int(line.split(": ")[1][:-1])


        line = fo.readline()
        line = fo.readline()
        line = fo.readline()
        line = fo.readline()
        line = fo.readline()
        line = fo.readline()
        line = fo.readline()
        line = fo.readline()
        line = fo.readline()
        line = fo.readline()
        line = fo.readline()
            
        line = fo.readline()
        data['time'] = int(line.split(": ")[1][:-1])

        line = fo.readline()
        data['relisted'] = int(line.split(": ")[1][:-1])

        line = fo.readline()
        list_string = line.split(": ")[1][1:-2].split(',')
        data['rank_freqs'] = [int(list_string[i]) for i in range(len(list_string))]

        fo.readline()
        fo.readline()
        fo.readline()

        line = fo.readline()
        sample_strings = []
        for n in range(500):
            string = line[12 + 3*n + 21*data['n']*n:12 + 3*n+21*data['n']*(n+1)][1:-1]
            ls = string.split('; ')
            # ls[0] = ls[0][1:]
            # ls[-1] = ls[-1][:-1]
            ls = [float(ls[j]) for j in range(len(ls))]
            sample_strings.append(ls)

    else:
        data['rng'] = line.split(": ")[0].split(" ")[0]
        data['seed'] = int(line.split(": ")[1][:-1])

        line = fo.readline()
        data['n'] = int(line.split(": ")[1][:-1])
        line = fo.readline()
        data['k'] = int(line.split(": ")[1][:-1])
            
        line = fo.readline()
        data['time'] = int(line.split(": ")[1][:-1])

        line = fo.readline()
        data['relisted'] = int(line.split(": ")[1][:-1])

        line = fo.readline()
        list_string = line.split(": ")[1][1:-2].split(',')
        data['rank_freqs'] = [int(list_string[i]) for i in range(len(list_string))]

        fo.readline()
        fo.readline()
        fo.readline()

        line = fo.readline()
        sample_strings = []
        for n in range(500):
            string = line[12 + 3*n + 21*data['n']*n:12 + 3*n+21*data['n']*(n+1)][1:-1]
            ls = string.split('; ')
            # ls[0] = ls[0][1:]
            # ls[-1] = ls[-1][:-1]
            ls = [float(ls[j]) for j in range(len(ls))]
            sample_strings.append(ls)
            

    # samples_raw = line.split("= ")[1][:-1]
    # sample_strings = samples_raw.split(',')
    # sample_strings[0] = sample_strings[0][1:]
    # sample_strings[-1] = sample_strings[-1][:-1]

    # for i in range(len(sample_strings)):
    #     string = sample_strings[i]
    #     string = string[1:-1]
    #     ls = string.split('; ')
    #     ls[0] = ls[0][1:]
    #     ls[-1] = ls[-1][:-1]
    #     ls = [float(ls[j]) for j in range(len(ls))]
    #     sample_strings[i] = ls

    data['samples'] = sample_strings
    fo.close()
    return data

def make_histogram(samples, bins, normalized=True):
    # given list of samples (real numbers) and number of bins, create list of numbers
    # where the ith value counts how many items from all samples are in the interval 
    # [i/bins, (i+1)/bins). e.g. if bins = 1000 then histogram[4] counts how many values
    # are in [4/1000, 5/1000). So more bins means more precise segregation of real values

    # if normalized, scales values so that their integral is 1, normalized by default

    nsamples = len(samples)
    population = len(samples[0])

    histogram = [0 for i in range(bins)]

    for i in range(nsamples):
        for j in range(population):
            histogram[int(bins*samples[i][j])] += 1

    if normalized:
        for i in range(bins):
            histogram[i] = float(histogram[i]) * (bins/(population*nsamples))
    
    return histogram

def hist_to_dist(histogram, population, nsamples):
    # given the output of make_histogram, rescales values so that their integral is 1

    bins = len(histogram)
    dist = [0 for i in range(bins)]
    for i in range(bins):
        dist[i] = (bins*float(histogram[i]))/(population*nsamples)
    
    return dist

def threshold_index(dist):
    # estimates threshold by finding the first index with frequency crossing 1.5
    # reports that index/bins e.g. if bins=1000 then may return 665/1000 = .665
    t = 0
    bins = len(dist)
    while(dist[t] < 1.5 and t <= bins):
        t += 1
    return t

def est_threshold_2(dist):
    # estimates threshold by computing average height of 'flat' area of distribution
    # (eyeballed between .7 and .9), then uses that height to compute the left edge of
    # the flat area, assuming that the true (in the limit) flat area is a rectangle with
    # height as estimated and area equal to 1

    h = 0
    bins = len(dist)
    for i in range(int((7/10)*bins), int((9/10)*bins)):
        h += dist[i]
    h = h / ((2/10)*bins)
    return round((1 - (1/h)), int(math.log10(bins)))

def est_thresholds(trials, bins):
    # given a bunch of files of trials identical but in seed
    # create list of all threshold indices calculated
    out = []
    for data in trials:
        dist = make_histogram(data['samples'], bins)
        out.append(threshold_index(dist)/bins)
    return out

def simple_plot(cells):
    ___, axs = plt.subplots()
    axs.scatter(range(len(cells)), cells, c = 'red', s = 5, alpha = .5)
    axs.set_ylim(min(cells), max(cells))
    axs.set_yscale('log')
    return plt

def make_pretty_plot(data, bins):
    # given data from read_file, and a number of bins, determines thresholds and
    # returns a plot of the fitnesses histogram as well as threshold values

    dist = make_histogram(data['samples'], bins)

    threshold = threshold_index(dist)/bins
    threshold_string = '{number:.{digits}f}'.format(number=threshold, digits=int(math.log10(bins)))

    ___, axs = plt.subplots()
    
    xs = [x/bins for x in range(bins)]

    axs.bar(xs, dist, width = 1/bins, color = 'grey', edgecolor = 'black')
    axs.set_ylim(min(dist), max(dist))
    axs.set_xlim(0, 1)


    res = int(bins/100) # bigger number is larger tick chunks
    # plt.xticks([xs[res*i] for i in range(int(len(xs)/res))]+[1.0])

    plt.xlabel('Value')
    plt.ylabel('Frequency') 
    plt.title("Fitnesses histogram for n = " + str(data['n']) + " and k = "+str(data['k'])+" with RNG = "+ str(data['rng']) + "\n x* ~= " + threshold_string)#+ ", runtime ~= " +str(np.round(data['time']/60, 2)) + " minutes")
    plt.axvline(x=threshold)
    #splt.text(threshold - .2 , 3, "x* ~= " + threshold_string)

    return plt

def mean(ls):
    return sum(ls)/len(ls)

def stdev(ls):
    m = mean(ls)
    return (sum([(ls[i] - m)**2 for i in range(len(ls))])/len(ls))**(.5)

def rescale_alphas(alphas, n): # given rank_freqs from run with n population, make the alphas as if sum of all of them (including the missing ones) were 3
    k = int(len(alphas)/3) # alphas list is always of length 3k
    tail_avg = mean(alphas[2*k:]) # what is a principled choice for this cutoff?
    return [alphas[i]/alphas[0] for i in range(len(alphas))]

def estimate_all_alphas(alphas, n): # given 3k rank_freqs as raw counts, extend the tail so that the length is n and also rescale so that the sum is 3
    k = int(len(alphas)/3)
    tail_avg = mean(alphas[2*k:])
    return [alphas[i]/alphas[0] for i in range(len(alphas))] + [tail_avg/alphas[0] for i in range(n - 3*k)]

def estimate_p(data):
    n = data['n']
    k = data['k']
    alphas = data['rank_freqs']

    # first get average number of replacements over interval (k, 3k]
    tail_avg = mean(alphas[k:3*k - 1])   #off by one because indexing
    # next compute an estimate of the total replacements, assuming a flat tail with value tail_avg
    total_replacements = sum(alphas[0:3*k - 1]) + (n - 3*k)*tail_avg       #off by one because indexing
    # need to rescale the tail average 
    rescaled_tail = tail_avg*3/total_replacements 
    # we can now estimate p from [veerman prieto 2018 pg 25]
    p = 3 - rescaled_tail*(n - k)

    # print(tail_avg)
    # print(total_replacements)
    # print(rescaled_tail)
    # print(p)

    return p

def get_p(files):
    out = []
    for i in range(len(files)):
        data = read_file(files[i])
        out.append(estimate_p(data))
    return out

def get_thresholds(files, bins):
    out = []
    for i in range(len(files)):
        out.append(threshold_index(make_histogram(read_file(files[i])['samples'], bins))/bins)
        print(out)
    return out

def estimate_kp(alphas, p): # given alphas as true frequencies, etimate k_p
    kp = 0
    kp_sum = 0
    while kp_sum <= p and kp < len(alphas):
        kp_sum += alphas[kp]
        kp += 1
    if kp == len(alphas):
        kp = -666*3

    return kp

def dels(alphas, p, n): # given alphas as true frequencies and a guess for p, compute dels
    # given an p, compute differences between true alphas and estimate based on p. We have that the true p is the smallest
    # p such that lim n->infty of n* +/-del = 0
    kp = estimate_kp(alphas, p)
    #print(kp)
    try:
        diffs = [alphas[i] - (3 - p)/(n - kp) for i in range(kp, len(alphas))]
    except:
        print("kp too large -- kp = " + str(kp) + ", 3k = " + str(len(alphas)))
        return 0
    return (max(diffs), min(diffs)) 

def possible_ps(data,ps):
    # gets a vertical slice of fig 4.6 in wagner, run this on a variety of ns to observe the limit in order to estimate p
    n = data['n']
    alphas = rescale_alphas(data['rank_freqs'], n)
    
    out = []
    
    for p in ps:
        ds = dels(alphas, p, n)
        if ds != 0:
            out.append((n, p, n*ds[0], n*ds[1]))
    return out

def plot_dels(ps_data, ps): # given the output of 'possible_ps' or union of such outputs, plot series of dels appropriately

    colors=iter(plt.cm.rainbow(np.linspace(0,1, len(ps)+1)))
    fig, axs = plt.subplots(2)

    pos_series = {}
    neg_series = {}
    for p in ps:
        c = next(colors)
        pts = [pt for pt in ps_data if pt[1] == p]
        pos_series[str(p)] = [(pt[0], pt[2]) for pt in pts]
        neg_series[str(p)] = [(pt[0], pt[3]) for pt in pts]
        axs[0].plot(*zip(*pos_series[str(p)]), color = c, label= "p = "+str(p), linestyle = '--', marker = 'o')       
        axs[1].plot(*zip(*neg_series[str(p)]), color = c, label= "p = "+str(p), linestyle = '--', marker = 'o')

    axs[0].legend(loc="upper left")
    axs[1].legend(loc="lower right")

    axs[0].set_xlabel("Population (n)")
    axs[0].set_ylabel("Mean n*delta_plus(n)")

    axs[1].set_xlabel("Population (n)")
    axs[1].set_ylabel("Mean n*delta_minus(n)")

    fig.suptitle("Do the rank replacement frequencies form a (p,q) array?")

    plt.show()


path = "C:\\Users\\Cameron\\OneDrive\\Math\\Research\\Bak Sneppen\\Code\\slurm_output"

ns = [25600]
bins = 1000000

files = {}
mt_files = {}
xos_files = {}

for n in ns:
    files[str(n)] = [item for item in os.scandir(path + "\\" + str(n))]
    mt_files[str(n)] = [item for item in files[str(n)] if item.name[4] == 'm']
    xos_files[str(n)] = [item for item in files[str(n)] if item.name[4] == 'x']

def p_q_array_test(files, ps):
    ps_data = []
    for n in ns:
        runs = mt_files[str(n)]
        all_ps_data_lists = [possible_ps(read_file(f), ps) for f in runs]
        all_ps_data = [item for sublist in all_ps_data_lists for item in sublist]
        # all this mess gets the estimated dels for all runs and then averages the del values among all runs with the same n
        mean_ps_data = [(n, p, mean([item[2] for item in all_ps_data if item[1] == p]), mean([item[3] for item in all_ps_data if item[1] == p])) for p in ps]
        ps_data += mean_ps_data

    plot_dels(ps_data, ps)

ps = np.linspace(2.000,2.0016,8)
def get_all_ks(files, ps): # given some files, and some ps you want to test, tells you what k value to set on future runs in order to test those ps

    all_kps = {}
    for n in ns:
        runs = files[str(n)]
        all_kps[str(n)] = []
        for run in runs:
            data = read_file(run)
            alphas = rescale_alphas(data['rank_freqs'], data['n'])
            all_kps[str(n)] += [(p, estimate_kp(alphas, p)/3) for p in ps]
    pp.pprint(all_kps)

#p_q_array_test(mt_files, ps)
#get_all_ks(mt_files, ps)

#data = read_file(mt_files['1600'][1])

#print(data['rank_freqs'])
#print(estimate_all_alphas(data['rank_freqs'], data['n']))

def plot_kde_byhand(samples):
    data = [item for sublist in samples for item in sublist]
    xs = np.linspace(0,1,1000)
    density = sum(norm(x).pdf(xs) for x in data[:1000])
    plt.fill_between(xs, density, alpha=0.5)
    #plt.plot(data, np.full_like(data, -0.1), '|k', markeredgewidth=1)
    plt.axis([0, 1, -0.2, 500]);
    plt.show()

def plot_kde(data):
    resolution = 10000
    samples = np.array([item for sublist in data['samples'] for item in sublist])
    xs = np.linspace(0,1,resolution)
    kde = KernelDensity(bandwidth=0.01, kernel='epanechnikov')
    kde.fit(samples[:, None])

    
    logprob = kde.score_samples(xs[:, None])
    threshold = threshold_index(np.exp(logprob))/resolution
    threshold_string = '{number:.{digits}f}'.format(number=threshold, digits=int(math.log10(resolution)))
    print(threshold_string)

    plt.fill_between(xs, np.exp(logprob), alpha=0.5)
    plt.axvline(x=threshold)
    
    plt.ylim(0, 3)

    plt.xlabel('Value')
    plt.ylabel('Frequency') 
    plt.title("Smoothed fitness distribution for n = " + str(data['n']) + " and k = "+str(data['k'])+" with RNG = "+ str(data['rng']) + "\n x* ~= " + threshold_string+ ", runtime ~= " +str(np.round(data['time']/60, 2)) + " minutes")
    plt.axvline(x=threshold)
    plt.show()


#simple_plot(rescale_alphas(data['rank_freqs'],data['n'])).show()
#plot_kde(data)



#---------------------------

# plt = simple_plot(rescale_alphas(data['rank_freqs']))
# plt.axvline(x=kp)
# plt.show()


#print(rescale_alphas(read_file(mt_files[0])['rank_freqs']))

# ps = get_p(mt_files)

# print(str(n))
# print("ps: " + str(ps))
# print("mean of ps: " + str(mean(ps)) )
# print("stdev of ps: "+ str(stdev(ps)))


mt_thresholds = get_thresholds(mt_files[str(n)], bins)
xos_thresholds = get_thresholds(xos_files[str(n)], bins)

# mt_trials1 = [read_file(mt_files[i]) for i in [0,1,2,3,4]]

print(str(n))
print("mt_thresholds: " + str(mt_thresholds))
print("mean of mt: " + str(round(mean([mt_thresholds[i] for i in range(len(mt_thresholds))]), 1 + int(math.log10(bins)))))
print("stdev of mt: " + str(round(stdev([mt_thresholds[i] for i in range(len(mt_thresholds))]), 1 + int(math.log10(bins)))))

print("xos_thresholds: " + str(xos_thresholds))
print("mean of xos: " + str(round(mean([xos_thresholds[i] for i in range(len(xos_thresholds))]), 1 + int(math.log10(bins)))))
print("stdev of xos: " + str(round(stdev([xos_thresholds[i] for i in range(len(xos_thresholds))]), 1 + int(math.log10(bins)))))

a=.665
b=.667
___, axs = plt.subplots()
axs.scatter(range(len(mt_thresholds)), mt_thresholds, c = 'red', s = 5, alpha = .5)
axs.set_ylim(a,b)
plt.hlines(2/3, 0,len(mt_thresholds), color='blue', linestyles='dotted')
plt.title("Distribution of thresholds computed from 20 computations\n rng: mt, n: "+str(ns[0]))
plt.savefig("C:/Users/Cameron/OneDrive/Math/Research/Bak Sneppen/Images/threshold distribution/"+str(ns[0])+"_mt_thresh.png")

___, axs = plt.subplots()
axs.scatter(range(len(xos_thresholds)), xos_thresholds, c = 'red', s = 5, alpha = .5)
axs.set_ylim(a,b)
plt.hlines(2/3, 0,len(xos_thresholds), color='blue', linestyles='dotted')
plt.title("Distribution of thresholds computed from 20 computations\n rng: xoshiro, n: "+str(ns[0]))
plt.savefig("C:/Users/Cameron/OneDrive/Math/Research/Bak Sneppen/Images/threshold distribution/"+str(ns[0])+"_xos_thresh.png")

# plt = simple_plot(rescale_alphas(read_file(mt_files[0])['rank_freqs']))
# plt.axvline(x=24)
# plt.show()

# data = read_file(mt_files[0])['rank_freqs']

# print(sum([data[i]/data[0] for i in range(len(data))]))

#print(len(read_file(mt_files[0])['rank_freqs']))




# plt = make_pretty_plot(read_file(mt_files[str(ns[0])][0]), bins)

# plt.savefig("C:/Users/Cameron/OneDrive/Math/Research/Bak Sneppen/Images/fitness hist/"+str(ns[0])+"_mt_hist.png")

# plt = make_pretty_plot(read_file(xos_files[str(ns[0])][0]), bins)


# plt.savefig("C:/Users/Cameron/OneDrive/Math/Research/Bak Sneppen/Images/fitness hist/"+str(ns[0])+"_xos_hist.png")