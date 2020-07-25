import matplotlib.pyplot as plt
import numpy as np
import math
import os


def read_file(f):
    # given a reference to a slurm output file
    # return a dict with all relevant fields
    # which are 
    # 'rng', 'seed', 'n', 'k', 'time', 'relisted', 'rank_freqs', 'samples'
    data = {}
    fo = open(f.path, 'r')

    line = fo.readline()
    data['rng'] = line.split(": ")[0].split(' ')[0]
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
    samples_raw = line.split("= ")[1][:-1]
    sample_strings = samples_raw.split(',')
    sample_strings[0] = sample_strings[0][1:]
    sample_strings[-1] = sample_strings[-1][:-1]

    for i in range(len(sample_strings)):
        string = sample_strings[i]
        string = string[1:-1]
        ls = string.split('; ')
        ls[0] = ls[0][1:]
        ls[-1] = ls[-1][:-1]
        ls = [float(ls[j]) for j in range(len(ls))]
        sample_strings[i] = ls

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
    axs.scatter(range(len(cells)), cells, c = 'red', s = 2, alpha = .5)
    axs.set_ylim(min(cells), max(cells))
    return plt

def make_pretty_plot(data, bins):
    # given data from read_file, and a number of bins, determines thresholds and
    # returns a plot of the fitnesses histogram as well as threshold values

    dist = make_histogram(data['samples'], bins)

    threshold = threshold_index(dist)/bins
    threshold_string = '{number:.{digits}f}'.format(number=threshold, digits=int(math.log10(bins)))

    ___, axs = plt.subplots()
    
    xs = [x/bins for x in range(bins)]

    axs.scatter(xs, dist, c = 'red', s = 2, alpha = .5)
    axs.set_ylim(min(dist), max(dist))
    axs.set_xlim(0, 1)

    res = int(bins/10) # bigger number is larger tick chunks
    plt.xticks([xs[res*i] for i in range(int(len(xs)/res))]+[1.0])

    plt.xlabel('Value')
    plt.ylabel('Frequency') 
    plt.title("Fitnesses histogram for n = " + str(data['n']) + " and k = "+str(data['k'])+" with RNG = "+ str(data['rng']) + "\n x* ~= " + threshold_string+ ", runtime ~= " +str(np.round(data['time']/60, 2)) + " minutes")
    plt.axvline(x=threshold)
    plt.text(threshold + .01 , .5, "x* ~= " + threshold_string)

    return plt

def mean(ls):
    return sum(ls)/len(ls)

def stdev(ls):
    m = mean(ls)
    return (sum([(ls[i] - m)**2 for i in range(len(ls))])/len(ls))**(.5)

def rescale_alphas(alphas):
    n = len(alphas)
    
    k = int(n/3)

    tail_avg = mean(alphas[2*k:]) # what is a principled choice for this cutoff?
    total_replacements = sum(alphas) + (n - 3*k)*tail_avg 
    scaling_factor = 3/total_replacements
    return [scaling_factor*alphas[i] for i in range(n)]

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

def estimate_kp(alphas, p):
    kp = 0
    kp_sum = 0
    while kp_sum <= p and kp < len(alphas):
        kp_sum += alphas[kp]
        kp += 1
    return kp


path = "C:\\Users\\Cameron\\OneDrive\\Math\\Research\\Bak Sneppen\\Code\\slurm_output"

n = 800
bins = 10000

files = [item for item in os.scandir(path + "\\" + str(n))]
mt_files = [item for item in files if item.name[4] == 'm']
#xos_files = [item for item in files if item.name[4] == 'x']

def get_thresholds(files, bins):
    out = []
    for i in range(len(files)):
        out.append(threshold_index(make_histogram(read_file(files[i])['samples'], bins))/bins)
    return out

data = read_file(mt_files[0])


kp = estimate_kp(rescale_alphas(data['rank_freqs']), 2.8)

def dels(alphas, p, n): 
    # given an p, compute differences between true alphas and estimate based on p. We have that the true p is the smallest
    # p such that lim n->infty of n* +/-del = 0
    kp = estimate_kp(alphas, p)
    diffs = [alphas[i] - (3 - p)/(n - kp) for i in range(kp, len(alphas))]
    return (max(diffs), min(diffs)) 

def possible_ps(data):
    # gets a vertical slice of fig 4.6 in wagner, run this on a variety of ns to observe the limit in order to estimate p
    alphas = rescale_alphas(data['rank_freqs'])
    n = data['n']
    out = []
    ps = np.linspace(1.5,2.5,3)
    for p in ps:
        ds = dels(alphas, p, n)
        out.append((n, p, n*ds[0], n*ds[1]))
    return out

#print(possible_ps(data))

# plt = simple_plot(rescale_alphas(data['rank_freqs']))
# plt.axvline(x=kp)
# plt.show()


#print(rescale_alphas(read_file(mt_files[0])['rank_freqs']))

# ps = get_p(mt_files)

# print(str(n))
# print("ps: " + str(ps))
# print("mean of ps: " + str(mean(ps)) )
# print("stdev of ps: "+ str(stdev(ps)))


#mt_thresholds = get_thresholds(mt_files, bins)

# mt_trials1 = [read_file(mt_files[i]) for i in [0,1,2,3,4]]

# print(str(n))
# print("mt_thresholds: " + str(mt_thresholds))
# print("mean of mt: " + str(round(mean([mt_thresholds[i] for i in range(len(mt_thresholds))]), 1 + int(math.log10(bins)))))
# print("stdev of mt: " + str(round(stdev([mt_thresholds[i] for i in range(len(mt_thresholds))]), 1 + int(math.log10(bins)))))

# print("xos_thresholds: " + str(xos_thresholds))
# print("mean of xos: " + str(round(mean([xos_thresholds[i] for i in range(len(xos_thresholds))]), 1 + int(math.log10(bins)))))
# print("stdev of xos: " + str(round(stdev([xos_thresholds[i] for i in range(len(xos_thresholds))]), 1 + int(math.log10(bins)))))

# plt = simple_plot(rescale_alphas(read_file(mt_files[0])['rank_freqs']))
# plt.axvline(x=24)
# plt.show()

data = read_file(mt_files[0])['rank_freqs']

print(sum([data[i]/data[0] for i in range(len(data))]))

#print(len(read_file(mt_files[0])['rank_freqs']))

# plt = make_pretty_plot(read_file(mt_files[0]), bins)
# plt.show()