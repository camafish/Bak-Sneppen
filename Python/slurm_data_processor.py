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

path = "C:\\Users\\Cameron\\OneDrive\\Math\\Research\\Bak Sneppen\\Code\\slurm_output"

n = 200
bins = 1000

files = [item for item in os.scandir(path + "\\" + str(n))]
mt_files = [item for item in files if item.name[4] == 'm']
#xos_files = [item for item in files if item.name[4] == 'x']

mt_trials = [read_file(mt_files[i]) for i in range(len(mt_files))]
#xos_trials = [read_file(xos_files[i]) for i in range(len(xos_files))]

mt_thresholds = est_thresholds(mt_trials, bins)
#xos_thresholds = est_thresholds(xos_trials, bins)

print(str(n))
print("mt_thresholds: " + str(mt_thresholds))
print("mean of mt: " + str(round(mean([mt_thresholds[i] for i in range(len(mt_thresholds))]), 1 + int(math.log10(bins)))))
print("stdev of mt: " + str(round(stdev([mt_thresholds[i] for i in range(len(mt_thresholds))]), 1 + int(math.log10(bins)))))

# print("xos_thresholds: " + str(xos_thresholds))
# print("mean of xos: " + str(round(mean([xos_thresholds[i] for i in range(len(xos_thresholds))]), 1 + int(math.log10(bins)))))
# print("stdev of xos: " + str(round(stdev([xos_thresholds[i] for i in range(len(xos_thresholds))]), 1 + int(math.log10(bins)))))

#simple_plot(make_histogram(mt_trials[0]['samples'], 1000, True))

# plt = make_pretty_plot(read_file(mt_files[1]), 10000)
# plt.show()