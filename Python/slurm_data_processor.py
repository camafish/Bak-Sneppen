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
    fo = open(f, 'r')

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

def make_histogram(samples, bins):
    # given list of samples (real numbers) and number of bins, create list of numbers
    # where the ith value counts how many items from all samples are in the interval 
    # [i/bins, (i+1)/bins). e.g. if bins = 1000 then histogram[4] counts how many values
    # are in [4/1000, 5/1000). So more bins means more precise segregation of real values
    nsamples = len(samples)
    pop = len(samples[0])

    histogram = [0 for i in range(bins)]

    for i in range(nsamples):
        for j in range(pop):
            histogram[int(bins*samples[i][j])] += 1
    
    return histogram

def hist_to_dist(histogram, population, nsamples):
    # given the output of make_histogram, rescales values so that their integral is 1

    bins = len(histogram)
    dist = [0 for i in range(bins)]
    for i in range(bins):
        dist[i] = (bins*float(histogram[i]))/(population*nsamples)
    
    return dist

def get_dist(data, bins):
    # given data dict, applies make_histogram and hist_to_dist to get a frequency distribution of fitness values
    # from the many samples found in data
    return hist_to_dist(make_histogram(data['samples'], bins), data['n'], len(data['samples']))

def est_threshold_1(dist):
    # estimates threshold by finding the first index with frequency crossing 1.5
    # reports that index/bins e.g. if bins=1000 then may return 665/1000 = .665
    t = 0
    bins = len(dist)
    while(dist[t] < 1.5 and t <= bins):
        t += 1
    return t/bins

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
    # create list of all pairs of thresholds
    out = []
    for data in trials:
        dist = get_dist(data, bins)
        out.append([est_threshold_1(dist), est_threshold_2(dist)])
    return out

def simple_plot(cells):
    fig, axs = plt.subplots()
    axs.scatter(range(len(cells)), cells, c = 'red', s = 2, alpha = .5)
    axs.set_ylim(0, 1)
    plt.show()

def mean(ls):
    return sum(ls)/len(ls)

def stdev(ls):
    m = mean(ls)
    return (sum([(ls[i] - m)**2 for i in range(len(ls))])/len(ls))**(.5)

path = "C:\\Users\\Cameron\\OneDrive\\Math\\Research\\Bak Sneppen\\Code\\C\\slurm_output"

files = [item for item in os.scandir(path + "\\6400")]
mt_files = [item for item in files if item.name[4] == 'm']
xos_files = [item for item in files if item.name[4] == 'x']

mt_trials = [read_file(mt_files[i]) for i in range(len(mt_files))]
xos_trials = [read_file(xos_files[i]) for i in range(len(xos_files))]

mt_thresholds = est_thresholds(mt_trials, 10000)
xos_thresholds = est_thresholds(xos_trials, 10000)

t = 1

print(mt_thresholds)
print(xos_thresholds)

print("mean of mt t2: " + str(mean([mt_thresholds[i][t] for i in range(len(mt_thresholds))])))
print("stdev of mt t2: " + str(stdev([mt_thresholds[i][t] for i in range(len(mt_thresholds))])))

print("mean of xos t2: " + str(mean([xos_thresholds[i][t] for i in range(len(xos_thresholds))])))
print("stdev of xos t2 " + str(stdev([xos_thresholds[i][t] for i in range(len(xos_thresholds))])))

