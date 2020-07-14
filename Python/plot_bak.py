import matplotlib.pyplot as plt
import csv
import sys
import numpy as np
csv.field_size_limit(sys.maxsize)

def static_plot(cells):
    fig, axs = plt.subplots()
    axs.scatter(range(len(cells)), cells, c = 'red', s = 2, alpha = .5)
    axs.set_ylim(0, 1)
    plt.show()

def display_histogram(state, method, iterations, pop, bns, runtime):
    n, bins, patches = plt.hist(x=state, bins=bns, color='#0504aa', alpha=0.7, rwidth=1, density=True)
    threshold1 = np.round(bins[next(x for x, val in enumerate(n) if val > 1.5)], decimals=5)
    h = sum(n[7*int(np.round(bns/10)):9*int(np.round(bns/10))])/(2*int(np.round(bns/10)))
    threshold2 = np.round(1 - (1/h), decimals = 5)
    
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Value')
    plt.ylabel('Frequency') 
    plt.title("Fitnesses histogram for n = " + str(pop) + " and "+str(iterations)+" iterations using "+ str(method) + "\n x1* ~= " + str(threshold1)+ ", x2* ~= "+ str(threshold2)+ ", runtime ~= " +str(np.round(runtime/60, 2)) + " minutes")
    plt.show()


def get_rows(filename):
    with open(filename) as csvfile:
        reader = csv.reader(csvfile)
        rows = [r for r in reader]
        return rows


def estimate_threshold(state, bns):
    n, bins, patches = plt.hist(x=state, bins=bns, color='#0504aa', alpha=0.7, rwidth=1, density=True)
    threshold1 = np.round(bins[next(x for x, val in enumerate(n) if val > 1.5)], decimals=5)
    h = sum(n[7*int(np.round(bns/10)):9*int(np.round(bns/10))])/(2*int(np.round(bns/10)))
    threshold2 = np.round(1 - (1/h), decimals = 5)

    return (threshold1, threshold2)

# rows = get_states()

# i=7

# samples_string = rows[i-1][8].strip('][').split(', ')
# samples = [[float(samples_string[j][2:-2].split(';')[i]) for i in range(int(rows[i-1][2]))] for j in range(len(samples_string))]

# display_histogram([item for sublist in samples for item in sublist], rows[i-1][1], rows[i-1][3],rows[i-1][2], 1000, float(rows[i-1][5]))

# rows = get_rows("C:\\Users\\Cameron\\OneDrive\\Math\\Research\\Bak Sneppen\\Code\\C\\c_bak_data.txt")
# state_string = rows[0][0]
# list_of_strings = state_string.strip('][').split(",")
# state = [float(list_of_strings[i]) for i in range(len(list_of_strings))]
# #print(state)
# static_plot(state)

rows = get_rows("C:\\Users\\Cameron\\OneDrive\\Math\\Research\\Bak Sneppen\\Code\\C\\c_bak_data_6_3.txt")


i = 140
rank_freqs_string = rows[i-1][7].strip('][').split(',')
rank_freqs = [int(rank_freqs_string[i]) for i in range(len(rank_freqs_string))]
print(rank_freqs)

fig, axs = plt.subplots()
axs.scatter(range(len(rank_freqs)), rank_freqs, c = 'red', s = 4, alpha = .5)
axs.set_ylim(0, int(rows[i-1][3])*1.1)
plt.show()



# for i in [118,119,120,121,122]:
#     samples_string = rows[i-1][6].strip('][').split(',')
#     samples = [[float(samples_string[j][2:-2].split('; ')[i]) for i in range(int(rows[i-1][2]))] for j in range(len(samples_string))]
#     thresholds = estimate_threshold([item for sublist in samples for item in sublist], 1000)
#     #display_histogram([item for sublist in samples for item in sublist], rows[i-1][1], rows[i-1][3],rows[i-1][2], 1000, float(rows[i-1][5]))
#     print((rows[i-1][2], thresholds))