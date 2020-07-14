import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import statistics as stat
#matplotlib.use("Agg")
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

import numpy as np
np.random.seed(0)
import fibheap as fh
from timeit import default_timer as timer
from datetime import datetime
import copy
import heapq
import sys
import csv
import pprint
pp = pprint.PrettyPrinter(indent=4, width = 4)
from datetime import date

## HELPER STUFF
def create_ecosystem_heap(vals):
    cells = [[0, 0, 0] for i in range(N)]
    cell_heap = []
    firstnode = [vals[0], 0, 0]
    heapq.heappush(cell_heap, firstnode) 
    previous = firstnode
    for i in range(1,N):
        newnode = [vals[i], i, 0]
        heapq.heappush(cell_heap, newnode)
        cells[i][0] = newnode
        cells[i][1] = previous
        cells[i-1][2] = newnode
        previous = newnode
    cells[N-1][2] = firstnode
    cells[N-1][0] = previous
    cells[0][1] = previous
    cells[0][0] = firstnode
    #print(*cells, sep  = "\n")
    return [cell_heap, cells]

def terminal_print(cell_list, k): # prints crude graph of cells, k is number of height buckets
    n = len(cell_list)
    img = [[' ' for i in range(k)] for i in range(n)]
    minindex = cell_list.index(min(cell_list))
    for i in range(n):
        for j in range(k):
            if cell_list[i] > j/k:
                if i == minindex:
                    img[i][j] = '#'
                else:
                    img[i][j] = '.'
    imgstr = ''
    for j in range(k):     
        imgstr += '|'  
        for i in range(n):
            imgstr += img[i][k - 1 - j]
        imgstr += '|\n'
    imgstr = imgstr[:-1]
    #print([round(cell_list[i], 2) for i in range(n)])
    print(imgstr)

def unheap(h):
    nodes = []  
    #print(h)
    while len(h) > 0:
        nodes.append(pop(h))
    return nodes

def make_heap(nodes):
    h = []
    for node in nodes:
        heapq.heappush(h, node)
    return h

def print_heap(heap_cells):
    out = []
    for i in range(len(heap_cells)):
        out.append(heap_cells[i][0][0])
    #print(out)
    return copy.deepcopy(out)

def pop(heap):
    node = heapq.heappop(heap)
    while node[2] == -100 and len(heap) > 0:
        node = heapq.heappop(heap)
        #print("REMOVED A NODE:" + str(node))
        #print(len(heap))
    return node

def popalot(heap, n):
    nodes = []  
    i = 0
    while  i < n:
        nodes.append(pop(heap))
        i += 1
    for node in nodes:
        heapq.heappush(heap, node)
    return heap

def update_node(heap_cells, node, val): # replaces a node by marking old one as trash and updating neighbors connections in heap_cells
    index = node[1]
    node[2] = -100
    newnode = [val, index, 0]
    heap_cells[index][0] = newnode
    heap_cells[(index + 1) % N][1] = newnode
    heap_cells[(index - 1) % N][2] = newnode
    return newnode #still needs to be pushed onto heap

def update_nodes(cells, index, radius, rng, changed): # takes cell list, central node index, radius, random number function, and list of changed
    nodes = cells[index]
    vals = [rng() for i in range(0,2*radius+1)] # e.g. radius 2 => [1,2,  3  ,4,5]
    #input(vals)

    # always update me
    update_node(cells, nodes[0], vals[radius])
    changed.append(nodes[0][1])

    i = 1
    left, right = nodes[1], nodes[2] # will spread out in pairs updating
    while i <= radius:
        update_node(cells, left, vals[radius-i])
        update_node(cells, right, vals[radius+i])
        changed.append(left[1] % N)
        changed.append(right[1] % N)
        left, right = cells[left[1]][1], cells[right[1]][2] # left is now left's left, e.g.
        i += 1

    return vals # returns new vals

def update(cells, index, newvals):
    cells[(index - 1) % N] = newvals[0] 
    cells[index] = newvals[1]
    cells[(index + 1) % N] = newvals[2]

def which_smallest(ls, n):
    lst = ls.copy()
    tmplst = ls.copy()
    out = []
    for i in range(n):
        mn = min(tmplst)
        minind = lst.index(mn)
        out.append(minind)
        tmplst.remove(mn)
    return out

def is_minheap(arr):
    return all(arr[i] >= arr[(i-1)//2] for i in range(1, len(arr)))

def rng_maker(kind, a, b): # returns a rng function with certain parameters
    if kind == 'real':
        def rng():
            return ((b - a) * np.random.random_sample() + a)      
    elif kind == 'int':
        def rng():
            return np.random.randint(low=a+1, high=b)
    return rng
    
def make_random_stream(kind, a, b, n):
    if kind == 'real':
        return ((b - a) * np.random.random_sample((1, n)) + a).tolist()[0]
    if kind == 'int':
        return np.random.randint(low=a, high=b,size=(1, n)).tolist()[0]

def my_insert(ls, new):
    if new < ls[0]:
        return [new] + ls
    elif new > ls[-1]:
        return ls + [new]
    else:
        i = 0
        while new > ls[i]:
            i += 1 # at the end of this i is the index of the val just above new
        return ls[:i] + [new] + ls[i:]

## PLOTTING
def static_plot(cells, data, save, name):
    fig, axs = plt.subplots()
    axs.scatter(range(len(cells)), cells, c = 'red', s = 2, alpha = .5)

    axs.set_ylim((data['range'][0] -.05, data['range'][1] + .05))
    
    plt.title("n = " +str(len(cells)) +" with "+str(data['maxiter'])+" iterations\ntime = " + str(round(data['time'], 4))+" seconds, "+str(data['method'].__name__))
    
    if save:
        plt.savefig(name + '.png')

    plt.show()
   
def animate_plot(history, data, speed, save, name):
    n = len(history[0])
    m = len(history)

    fig, ax = plt.subplots()  
    sc, = ax.plot(range(n), history[0], 'rx', markersize = 3, ls="") 
    ax.set_ylim((data['range'][0] -.05, data['range'][1] + .05))
    ax.set_title("n = " +str(n) +" with "+str(data['maxiter'])+" iterations at "+str(data['framerate'])+" step intervals\ntime = " + str(round(data['time'], 4))+" seconds, "+str(data['method'].__name__))

    def plot(a):
        sc.set_data(range(n), history[a])

    anim = animation.FuncAnimation(fig, plot, frames=m, interval=speed, repeat_delay=1000) 

    if save:
        anim.save(name +'.mp4', writer='ffmpeg', fps=10)

    plt.show()

## ITERATING FUNCTIONS
def list_iterate(cells, rng):                                
    update(cells, cells.index(min(cells)), [rng() for i in range(3)])               

def heap_iterate(heap, heap_cells, rng):
    newvals = [rng() for i in range(3)]      
    if heap[0]:
        minnode = pop(heap)      
        minnode[0] = newvals[1]
        heapq.heappush(heap, minnode)
        
        nodes = heap_cells[minnode[1]]
        heapq.heappush(heap, update_node(heap_cells, nodes[1], newvals[0]))
        heapq.heappush(heap, update_node(heap_cells, nodes[2], newvals[2]))

def sub_avalanche(cells, index, rng, current_iteration, maxiter):
    sub_changed = set()
    length_of_sub_avalanche = 0

    old_min = cells[index]
    old_min_index = index
    
    newvals = [rng() for i in range(3)] # roll three new numbers
    smallest = min(newvals) # determine the new smallest one
    smallest_index = (old_min_index + newvals.index(smallest) - 1) % N # make note of its absolute index
    
    update(cells, old_min_index, newvals)

    sub_changed = sub_changed.union({old_min_index, (old_min_index + 1) % N, (old_min_index - 1) % N}) # update values and store which indices changed
    length_of_sub_avalanche += 1

    while smallest < old_min and current_iteration + length_of_sub_avalanche + 1 < maxiter:
        old_min = smallest
        old_min_index = smallest_index 
        
        newvals = [rng() for i in range(3)] # roll three new numbers
        smallest = min(newvals) # determine the new smallest one
        smallest_index = (old_min_index + newvals.index(smallest) - 1) % N # make note of its absolute index
        
        update(cells, old_min_index, newvals)

        sub_changed = sub_changed.union({old_min_index, (old_min_index + 1) % N, (old_min_index - 1) % N}) # update values and store which indices changed
        length_of_sub_avalanche += 1
 
    if AVALANCHE_REPORT:
        print("     sub avalanche of length " + str(length_of_sub_avalanche))

    return (sub_changed, length_of_sub_avalanche)

def avalanche(cells, rng, current_iteration, maxiter):
    changed = set()
    length_of_avalanche = 0

    if AVALANCHE_REPORT:
        input(' ')
        print("Composite avalanche START at time "+str(current_iteration))
    
    original_min = min(cells)
    original_min_index = cells.index(original_min)

    sub_result = sub_avalanche(cells, original_min_index, rng, current_iteration, maxiter)
    changed = changed.union(sub_result[0])
    length_of_avalanche += sub_result[1]

    smallest_changed = min([(cells[i], i) for i in changed]) # check the list of changed fitnessess to see if any are less than original_min, if so then the smallest is the new global min
    smallest_changed_index = smallest_changed[1] 

    while original_min > smallest_changed[0] and current_iteration + length_of_avalanche + 1 < maxiter:
        sub_result = sub_avalanche(cells, smallest_changed_index, rng, current_iteration + length_of_avalanche, maxiter)
        changed = changed.union(sub_result[0])
        length_of_avalanche += sub_result[1]

        smallest_changed = min([(cells[i], i) for i in changed]) # check the list of changed fitnessess to see if any are less than original_min, if so then the smallest is the new global min
        smallest_changed_index = smallest_changed[1] 

    if AVALANCHE_REPORT:
        print("Composite avalanche END with length "+str(length_of_avalanche))
        print(" cells changed: " + str(changed))

    return length_of_avalanche

def sub_heap_avalanche(cells, index, rng, current_iteration, maxiter):
    sub_changed = set()
    length_of_sub_avalanche = 0

    old_minnode = cells[index]
    old_min = old_minnode[0][0]
    old_min_index = index

    newvals = [rng() for i in range(3)] # roll three new numbers
    smallest = min(newvals) # determine the new smallest one
    smallest_index = (old_min_index + newvals.index(smallest) - 1) % N # make note of its absolute index
    
    update_node(cells, old_minnode[1], newvals[0])
    update_node(cells, old_minnode[0], newvals[1])
    update_node(cells, old_minnode[2], newvals[2])

    sub_changed = sub_changed.union({old_min_index, (old_min_index + 1) % N, (old_min_index - 1) % N})
    length_of_sub_avalanche += 1

    while smallest < old_min and current_iteration + length_of_sub_avalanche + 1 < maxiter:
        old_minnode = cells[smallest_index] #get a reference to the right min node
        old_min = smallest
        old_min_index = smallest_index 
        
        newvals = [rng() for i in range(3)] # roll three new numbers
        smallest = min(newvals) # determine the new smallest one
        smallest_index = (old_min_index + newvals.index(smallest) - 1) % N # make note of its absolute index

        update_node(cells, old_minnode[1], newvals[0])
        update_node(cells, old_minnode[0], newvals[1])
        update_node(cells, old_minnode[2], newvals[2])

        sub_changed = sub_changed.union({old_min_index, (old_min_index + 1) % N, (old_min_index - 1) % N})
        length_of_sub_avalanche += 1

    return (sub_changed, length_of_sub_avalanche)

def heap_avalanche(cells, minnode, rng, current_iteration, maxiter):
    changed = set()
    length_of_avalanche = 0

    original_min = minnode[0] # find the global minimum, save for later 
    original_min_index = minnode[1] # and its index

    sub_result = sub_heap_avalanche(cells, original_min_index, rng, current_iteration, maxiter)
    changed = changed.union(sub_result[0])
    length_of_avalanche += sub_result[1]

    smallest_changed = min([(cells[i][0][0], i) for i in changed]) # check the list of changed fitnessess to see if any are less than original_min, if so then the smallest is the new global min
    smallest_changed_index = smallest_changed[1] 

    while original_min > smallest_changed[0] and current_iteration + length_of_avalanche + 1 < maxiter:
        sub_result = sub_heap_avalanche(cells, smallest_changed_index, rng, current_iteration + length_of_avalanche, maxiter)
        changed = changed.union(sub_result[0])
        length_of_avalanche += sub_result[1]

        smallest_changed = min([(cells[i][0][0], i) for i in changed])  # check the list of changed fitnessess to see if any are less than original_min, if so then the smallest is the new global min
        smallest_changed_index = smallest_changed[1] 
        
    return (changed, length_of_avalanche)

def shortlist_sub_avalanche(cells, shortlist, rng, current_iteration, maxiter):
    sub_changed = set()
    length_of_sub_avalanche = 0
    threshold = shortlist[-1][0][0]
    old_minnode = shortlist[0] 

    old_min = old_minnode[0][0]
    old_min_index = old_minnode[0][1]
    
    newvals = [rng() for i in range(3)] # roll three new numbers
    smallest = min(newvals) # determine the new smallest one
    smallest_index = (old_min_index + newvals.index(smallest) - 1) % N # make note of its absolute index
    
    cells[(old_min_index - 1) % N][0][0] = newvals[0]
    cells[(old_min_index) % N][0][0] = newvals[1]
    cells[(old_min_index + 1) % N][0][0] = newvals[2]

    sub_changed = sub_changed.union({old_min_index, (old_min_index + 1) % N, (old_min_index - 1) % N}) # update values and store which indices changed
    length_of_sub_avalanche += 1

    while smallest < old_min and current_iteration + length_of_sub_avalanche + 1 < maxiter:
        old_min = smallest
        old_min_index = smallest_index 
        
        newvals = [rng() for i in range(3)] # roll three new numbers
        smallest = min(newvals) # determine the new smallest one
        smallest_index = (old_min_index + newvals.index(smallest) - 1) % N # make note of its absolute index
        
        cells[(old_min_index - 1) % N][0][0] = newvals[0]
        cells[(old_min_index) % N][0][0] = newvals[1]
        cells[(old_min_index + 1) % N][0][0] = newvals[2]

        sub_changed = sub_changed.union({old_min_index, (old_min_index + 1) % N, (old_min_index - 1) % N}) # update values and store which indices changed
        length_of_sub_avalanche += 1
 
    if AVALANCHE_REPORT:
        print("     sub avalanche of length " + str(length_of_sub_avalanche))

    for index in sub_changed:
        if cells[index] in shortlist:
            shortlist.remove(cells[index])

    for index in sub_changed:            
        if cells[index][0][0] < threshold:
            shortlist = my_insert(shortlist, cells[index])

    return (shortlist, sub_changed, length_of_sub_avalanche)

def shortlist_update(cells, shortlist, newvals): # shortlist is sorted smallest cells
    threshold = shortlist[-1][0][0]
    minnode = shortlist[0]
    min_index = minnode[0][1]

    # now for each of the three nodes, I need to compare with threshold to see about whether they stay in shortlist
    
    # first for the middle node
    minnode[0][0] = newvals[1] # the new value for this node

    shortlist.remove(minnode)       # it's in shortlist, so remove it
    if minnode[0][0] < threshold:   # if under threshold, insert
        shortlist = my_insert(shortlist, minnode) 

    # left node
    left = cells[(min_index - 1) % N]
    left[0][0] = newvals[0]

    if left in shortlist:           # if it's in shortlist, remove it
        shortlist.remove(left)
    if left[0][0] < threshold:      # if it should be in shortlist, insert
        shortlist = my_insert(shortlist, left)

    # right node
    right = cells[(min_index + 1) % N]
    right[0][0] = newvals[2]

    if right in shortlist:
        shortlist.remove(right)
    if right[0][0] < threshold:
        shortlist = my_insert(shortlist, right)

    #assert shortlist == sorted(shortlist)
    return shortlist
 
def shortlist_heap_iterate(cells, shortlist, rng): 
    net_length = 0 # want to keep track of number of real values in shortlist
    threshold = max(shortlist)[0] # just have to search for it
    minnode = pop(shortlist) # it's on the bottom so pop it
    
    min_index = minnode[1]

    newvals = [rng() for i in range(3)]   

    # now for each of the three nodes, I need to compare with threshold to see about whether they stay in shortlist
    
    # first for the middle node
    new_min = update_node(cells, minnode, newvals[1]) # note that update_node marks the old one for removal
    net_length -= 1
    if new_min[0] < threshold:       
        heapq.heappush(shortlist, new_min)
        net_length += 1

    # left node
    left = cells[(min_index - 1) % N][0]
    new_left = update_node(cells, left, newvals[0]) # note that update_node marks the old one for removal
    if left in shortlist:
        net_length -= 1
    if new_left[0] < threshold: # if it should be in shortlist, insert
        heapq.heappush(shortlist, new_left)
        net_length += 1

    # right node
    right = cells[(min_index + 1) % N][0]
    new_right = update_node(cells, right, newvals[2])
    if right in shortlist:
        net_length -= 1
    if new_right[0] < threshold: # if it should be in shortlist, insert
        heapq.heappush(shortlist, new_right)
        net_length += 1

    # assert is_minheap(shortlist)

    # heap = copy.deepcopy(shortlist)
    # heapq.heapify(heap)
    # assert heap == shortlist

    return (shortlist, net_length)

## METHODS

def list_method(data, i, history): # naive method
    #cells = data['cell_data']

    list_iterate(data['cell_data'], data['rng'])   

    if (data['framerate'] > 0 and i % data['framerate'] == 0):
        print(str(data['method'].__name__) + " Progress: " + str(100*i/data['maxiter'])+'%')
        if data['save']:
            history.append(data['cell_data'].copy())

    if not data['save'] and i + 1 == data['maxiter']:
        history.append(data['cell_data'].copy())

    return 1

def heap_method(data, i, history): # maintaining entire state as heap
    heap = data['cell_data'][0]
    heap_cells = data['cell_data'][1]

    if len(heap) > HEAP_LIMIT_FACTOR*data['popnum']:
        heap = popalot(heap, data['popnum']) 

    heap_iterate(heap, heap_cells, data['rng'])   

    # if data['reheap_timer'] > 0 and i % data['reheap_timer'] == 0:
    #     heap = popalot(heap, data['popnum']) 

    if (data['framerate'] > 0 and i % data['framerate'] == 0) or i > data['maxiter']:
        print(str(data['method'].__name__) + " Progress: " + str(100*i/data['maxiter'])+'%'+  " - heapsize: "+str(len(heap)))

        if data['save']:
            history.append(print_heap(heap_cells))

        #print(heap[0])
    if not data['save'] and i + 1 == data['maxiter']:
        history.append(print_heap(heap_cells))

    return 1

def avalanche_method(data, i, history): # following avalanches to find next minimum, storing state as list, iterative, and with compount avalanches
    cells = data['cell_data']

    if i + 1 < data['maxiter']:
        c = avalanche(cells, data['rng'], i, data['maxiter'])
        
    elif i < data['maxiter'] + 1:
        #print(i)
        list_iterate(data['cell_data'], data['rng'])   
        c = 1
    else:
        c = 0

    frame = i + c # tracks how many iterations the above function completed, i is which frame we started on in this call

    #data['avalengths'].append(c)

    if (data['framerate'] > 0 and frame % data['framerate'] in range(c)): # note that an avalanche can start when frame is less than maxiter, then it can continue, not a big deal, could pull i through the call?
        print(str(data['method'].__name__) + " Progress: " + str(100*frame/data['maxiter'])+'%')
        #print("   composite avalanche length " + str(c))
        if data['save']:
            history.append(cells.copy())

    if not data['save'] and frame == data['maxiter']:
        history.append(cells.copy())

    return c

def avalanche_single_method(data, i, history): # following avalanches to find next minimum, storing state as list, iterative, and with compount avalanches
    cells = data['cell_data']

    if i + 1 < data['maxiter']:
        __ , c = sub_avalanche(cells, cells.index(min(cells)), data['rng'], i, data['maxiter'])
    elif i < data['maxiter'] + 1:
        list_iterate(data['cell_data'], data['rng'])   
        c = 1
    else:
        c = 0

    frame = i + c # tracks how many iterations the above function completed, i is which frame we started on in this call

    data['avalengths'].append(c)

    if (data['framerate'] > 0 and frame % data['framerate'] in range(c)): # note that an avalanche can start when frame is less than maxiter, then it can continue, not a big deal, could pull i through the call?
        print(str(data['method'].__name__) + " Progress: " + str(100*frame/data['maxiter'])+'%')
        #print("   composite avalanche length " + str(c))
        if data['save']:
            history.append(cells.copy())

    if not data['save'] and frame == data['maxiter']:
        history.append(cells.copy())

    return c

def heap_avalanche_method(data, i, history):     
    heap = data['cell_data'][0]
    heap_cells = data['cell_data'][1]

    if len(heap) > HEAP_LIMIT_FACTOR*data['popnum']:
        heap = popalot(heap, data['popnum']) 

    if i  < data['maxiter']:
        changed, c = heap_avalanche(heap_cells, pop(heap), data['rng'], i, data['maxiter'])
        for n in changed:
            heapq.heappush(heap, heap_cells[n][0])
    elif i < data['maxiter'] + 1:
        list_iterate(print_heap(heap_cells), data['rng'])   
        c = 1
    else:
        c = 0

    frame = i + c # tracks how many iterations the above function completed, i is which frame we started on in this call

    # if data['reheap_timer'] > 0 and frame % data['reheap_timer'] in range(1,c):     
    #     heap = popalot(heap, data['popnum'])

    if (data['framerate'] > 0 and frame % data['framerate'] in range(c)): # note that an avalanche can start when frame is less than maxiter, then it can continue, not a big deal, could pull i through the call?
        print(str(data['method'].__name__) + " Progress: " + str(100*frame/data['maxiter'])+'%'+  " - heapsize: "+str(len(heap)))
        #print("   composite avalanche length " + str(c))

        if data['save']:
            history.append(print_heap(heap_cells))

    if not data['save'] and frame == data['maxiter']:
        history.append(print_heap(heap_cells))

    return c

def heap_avalanche_single_method(data, i, history):     
    heap = data['cell_data'][0]
    heap_cells = data['cell_data'][1]

    if len(heap) > HEAP_LIMIT_FACTOR*data['popnum']:
        heap = popalot(heap, data['popnum']) 

    if i < data['maxiter']:
        changed, c = sub_heap_avalanche(heap_cells, pop(heap)[1], data['rng'], i, data['maxiter'])
        for n in changed:
            heapq.heappush(heap, heap_cells[n][0])
    elif i < data['maxiter'] + 1:
        list_iterate(print_heap(heap_cells), data['rng'])   
        c = 1
    else:
        c = 0

    frame = i + c # tracks how many iterations the above function completed, i is which frame we started on in this call

    # if data['reheap_timer'] > 0 and frame % data['reheap_timer'] in range(1,c):     
    #     heap = popalot(heap, data['popnum'])

    if (data['framerate'] > 0 and frame % data['framerate'] in range(c)): # note that an avalanche can start when frame is less than maxiter, then it can continue, not a big deal, could pull i through the call?
        print(str(data['method'].__name__) + " Progress: " + str(100*frame/data['maxiter'])+'%'+  " - heapsize: "+str(len(heap)))
        #print("   composite avalanche length " + str(c))
        
        if data['save']:
            history.append(print_heap(heap_cells))

    if not data['save'] and frame == data['maxiter']:
        history.append(print_heap(heap_cells))
        
    return c

def shortlist_method(data, i, history):
    
    if len(data['shortlist']) < data['shortlength'][0]: #data['shortlength'] = (3*k_p, 5*k_p)
        data['shortlist'] = sorted(data['cell_data'])[:data['shortlength'][1]]
        data['relisted'] += 1
        print("Relisted " + str(data['relisted']) + " times")
        #print(i)

    data['shortlist'] = shortlist_update(data['cell_data'], data['shortlist'], [data['rng']() for i in range(3)]   ) # the 'rng' part is giving it three random numbers with the passed rng function
   # assert data['shortlist'] ==  sorted(data['cell_data'])[:len(data['shortlist'])]

    if (data['framerate'] > 0 and i % data['framerate'] == 0):
        print(str(data['method'].__name__) + " Progress: " + str(100*i/data['maxiter'])+'%')
        if data['save']:
            history.append(print_heap(data['cell_data']))

    if not data['save'] and i + 1 == data['maxiter']:
        history.append(print_heap(data['cell_data']))

    return 1

def shortlist_avalanche_single_method(data, i, history):
    
    if len(data['shortlist']) < data['shortlength'][0]: #data['shortlength'] = (3*k_p, 5*k_p)
        data['shortlist'] = sorted(data['cell_data'])[:data['shortlength'][1]]
        data['relisted'] += 1
        print("Relisted " + str(data['relisted']) + " times")
        #print(i)

    if i < data['maxiter']:
        data['shortlist'], __, c = shortlist_sub_avalanche(data['cell_data'],  data['shortlist'], data['rng'], i, data['maxiter'])
    elif i < data['maxiter'] + 1:
        list_iterate(print_heap(data['cell_data']), data['rng'])   
        c = 1
    else:
        c = 0

    frame = i + c

   # assert data['shortlist'] ==  sorted(data['cell_data'])[:len(data['shortlist'])]

    if (data['framerate'] > 0 and frame % data['framerate'] in range(c)):
        print(str(data['method'].__name__) + " Progress: " + str(100*i/data['maxiter'])+'%')
        if data['save']:
            history.append(print_heap(data['cell_data']))

    if not data['save'] and frame == data['maxiter']:
        history.append(print_heap(data['cell_data']))

    return c

def shortlist_heap_method(data, i, history): # maintaining shortlist as a heap

    if data['true_length'] < data['shortlength'][0]: #data['shortlength'] = (3*k_p, 5*k_p)
        #print(data['true_length'])
        data['shortlist'] = [(sorted(data['cell_data'])[:data['shortlength'][1]])[i][0] for i in range(data['shortlength'][1])]
        heapq.heapify(data['shortlist'])
        data['true_length'] = data['shortlength'][1]

        data['relisted'] += 1
        print("Relisted " + str(data['relisted']) + " times")
        
    if len(data['shortlist']) > HEAP_LIMIT_FACTOR*data['true_length']:
        data['shortlist'] = popalot(data['shortlist'], data['true_length']) 

    data['shortlist'], net_length = shortlist_heap_iterate(data['cell_data'], data['shortlist'], data['rng']) # the 'rng' part is giving it three random numbers with the passed rng function
    data['true_length'] += net_length

    # if data['reheap_timer'] > 0 and i % data['reheap_timer'] == 0:
    #     data['shortlist'] = popalot(data['shortlist'], data['true_length']) 

    if (data['framerate'] > 0 and i % data['framerate'] == 0):
        print(str(data['method'].__name__) + " Progress: " + str(100*i/data['maxiter'])+'%' +  " - heapsize: "+str(len(data['shortlist'])))

        if data['save']:
            history.append(print_heap(data['cell_data']))

    if not data['save'] and i + 1 == data['maxiter']:
        history.append(print_heap(data['cell_data']))

    return 1

## RUN
def run(data):
    data['history'] = [copy.deepcopy(data['cell_data'])]
    history = []
    if data['framerate'] == -1: #can use -1 for step by step running
        i = 0
        while input() == '':
            i += data['method'](data, i, history)   
    else:
        start = timer()
        i = 0 # i tracks how many iterations
        while i < data['maxiter']:                 
            i += data['method'](data, i, history)
            #print(i)

        end = timer()
        data['time'] = end - start
        print(str(data['method'].__name__)+" time: " + str(end - start))

    print()
    data['history'].extend(copy.deepcopy(history))
    return(copy.deepcopy(data))

N = 110
MAXITER = N**3
k_p = 5 #int(N**.5)#int(N / np.log(N))
#print(k_p)

FRAMERATE = MAXITER//100
HEAP_LIMIT_FACTOR = 3

NUM_TYPE = 'real'

RANGE = (0,1)
RADIUS = 1 # only chain method does this right now 3/4/20

RUN = 0
PLOT = 0
SAVE = 0
SAVEPLOT = 0
PLOTNAME = 'plt'
AVALANCHE_REPORT = 0

basic_data =  {
    'maxiter' : MAXITER, 
    'framerate' : FRAMERATE, 
    'reheap_timer' : 500,  ## around 500 seems to work # defunct now 5/6
    'popnum' : N,
    'save' : SAVE, 
    'range' : RANGE,
    'rng' : rng_maker(NUM_TYPE, RANGE[0], RANGE[1]),
    'radius': RADIUS
    }
    
def test_methods(basic_data, methods):
    np.random.seed(0) # set seed 0 for creating initial vals
    combo_data = {methods[i]:copy.deepcopy(basic_data) for i in range(len(methods))}
    init_vals = make_random_stream(NUM_TYPE, RANGE[0], RANGE[1], N) # List of N initial fitnesses
    init_heap = create_ecosystem_heap(init_vals) # each heap method will use this initial heap

    if 'list' in methods:
        this = 'list'
        combo_data[this]['cell_data'] = copy.deepcopy(init_vals)
        combo_data[this]['method'] = list_method

    if 'heap' in methods:
        this = 'heap'
        combo_data[this]['cell_data'] = copy.deepcopy(init_heap)
        combo_data[this]['method'] = heap_method

    if 'avalanche' in methods:
        this = 'avalanche'
        combo_data[this]['cell_data'] = copy.deepcopy(init_vals)
        combo_data[this]['method'] = avalanche_method   
        combo_data[this]['avalengths'] = []

    if 'avalanche_single' in methods:
        this = 'avalanche_single'
        combo_data[this]['cell_data'] = copy.deepcopy(init_vals)
        combo_data[this]['method'] = avalanche_single_method   
        combo_data[this]['avalengths'] = []

    if 'heap_avalanche' in methods:
        this = 'heap_avalanche'
        combo_data[this]['cell_data'] = copy.deepcopy(init_heap)
        combo_data[this]['method'] = heap_avalanche_method

    if 'heap_avalanche_single' in methods:
        this = 'heap_avalanche_single'
        combo_data[this]['cell_data'] = copy.deepcopy(init_heap)
        combo_data[this]['method'] = heap_avalanche_single_method

    if 'shortlist' in methods:
        this = 'shortlist'
        combo_data[this]['cell_data'] = copy.deepcopy(init_heap)[1]
        combo_data[this]['shortlength'] = (3*k_p, 5*k_p)
        combo_data[this]['shortlist'] = sorted(combo_data[this]['cell_data'])[:combo_data[this]['shortlength'][1]]
        combo_data[this]['method'] = shortlist_method
        combo_data[this]['relisted'] = 0
    
    if 'shortlist_avalanche_single' in methods:
        this = 'shortlist_avalanche_single'
        combo_data[this]['cell_data'] = copy.deepcopy(init_heap)[1]
        combo_data[this]['shortlength'] = (3*k_p, 5*k_p)
        combo_data[this]['shortlist'] = sorted(combo_data[this]['cell_data'])[:combo_data[this]['shortlength'][1]]
        combo_data[this]['method'] = shortlist_avalanche_single_method
        combo_data[this]['relisted'] = 0

    if 'shortlist_heap' in methods:
        this = 'shortlist_heap'
        combo_data[this]['cell_data'] = copy.deepcopy(init_heap)[1]
        combo_data[this]['shortlength'] = (3*k_p, 5*k_p)
        combo_data[this]['shortlist'] = [(sorted(combo_data[this]['cell_data'])[:combo_data[this]['shortlength'][1]])[i][0] for i in range(combo_data[this]['shortlength'][1])]
        heapq.heapify(combo_data[this]['shortlist'])
        combo_data[this]['true_length'] = combo_data[this]['shortlength'][1]
        combo_data[this]['method'] = shortlist_heap_method      
        combo_data[this]['relisted'] = 0

    for method in methods:
        np.random.seed(1) # reset seed to 1 to run each method
        combo_data[method] = run(combo_data[method])

    time = datetime.now()

    with open(r'bak_data.csv', 'a', newline='') as csvfile:
        logger = csv.writer(csvfile, delimiter=',')
        for method in methods:
            try:
                logger.writerow([time, method, combo_data[method]['popnum'], combo_data[method]['maxiter'], k_p, combo_data[method]['time'], init_vals, combo_data[method]['history'][-1]]) 
                #print(combo_data[method]['history'][-1])
                print(method + ' ' + str(combo_data[method]['time']))
            except:
                print('error in logging '+str(method))

    assert(len((set([tuple(combo_data[method]['history'][-1]) for method in methods])))==1) # This checks if all the endstates match for different methods
    print("Done!")
    return combo_data

methods = [
    #  'list',
    #  'heap',
    #  'avalanche',
    #  'avalanche_single',
    #  'heap_avalanche',
    # 'heap_avalanche_single',
    'shortlist',
    'shortlist_avalanche_single',
    'shortlist_heap'
    ]

results = test_methods(basic_data, methods)



def display_avalengths(results):
    lengths = results['avalanche']['avalengths']
    name = "avalanches"
    fig, axs = plt.subplots()
    axs.scatter(range(len(lengths)), lengths, c = 'red', s = 5, alpha = .5)
    plt.yscale('log')
    plt.title("Length of Composite Avalanches for n = " +str(N) + " and "+str(MAXITER)+" iterations")
    plt.savefig(name + '.png')
    plt.show()

#display_avalengths(results)


if RUN:
    basic_data =  {
        'maxiter' : MAXITER, 
        'framerate' : FRAMERATE, 
        'reheap_timer' : 1000,  ## around 500 seems to work
        'popnum' : N,
        'save' : ANIMATE, 
        'range' : RANGE,
        'rng' : rng_maker(NUM_TYPE, RANGE[0], RANGE[1]),
        'radius': RADIUS
        }

    #sys.setrecursionlimit(MAXITER) ## Stackoverflow when heapsize is about 573861 to 675758

    data = copy.deepcopy(basic_data)
    init_vals = make_random_stream(NUM_TYPE, RANGE[0], RANGE[1], N)

    if CHOICE == 'heap':
        heap = create_ecosystem_heap(init_vals)

        data['cell_data'] = copy.deepcopy(heap)
        data['method'] = heap_method
        history, data = run(data)

    elif CHOICE == 'chain':
        heap = create_ecosystem_heap(init_vals)

        data['cell_data'] = copy.deepcopy(heap)
        data['method'] = chain_method
        history, data = run(data)

    elif CHOICE == 'list':

        data['cell_data'] = copy.deepcopy(init_vals)
        data['method'] = list_method
        history, data = run(data)

    elif CHOICE == 'shortlist_heap':
    
        data['cell_data'] = copy.deepcopy(create_ecosystem_heap(init_vals))
        data['shortlength'] = SHORTLENGTH
        data['shortlist'] = heapq.nsmallest(data['shortlength'], data['cell_data'][0], key = lambda x: x[0])
        heapq.heapify(data['shortlist'])
        data['method'] = shortlist_heap_method      
        data['relisted'] = 0
    
        history, data = run(data)
        print("relisted " + str(data['relisted']) + " times after " + str(data['maxiter']) +" iterations")

    elif CHOICE == 'shortlist':
    
        data['cell_data'] = copy.deepcopy(create_ecosystem_heap(init_vals))
        data['shortlength'] = SHORTLENGTH
        data['shortlist'] = heapq.nsmallest(data['shortlength'], data['cell_data'][0], key = lambda x: x[0])
        data['method'] = shortlist_method      
        data['relisted'] = 0
    
        history, data = run(data)
        print("relisted " + str(data['relisted']) + " times after " + str(data['maxiter']) +" iterations")

    time = data['time']
    method = data['method']

    #print(str(list_method.__name__)+" time: " + str(list_time))
    print(str(method.__name__)+" time: " + str(time))
    #print(str(chain_method.__name__)+" time: " + str(chain_time))


    if PLOT:
        if ANIMATE:
            animate_plot(history, data, 100, SAVEPLOT, PLOTNAME)
        else:
            static_plot(history[-1], data, SAVEPLOT, PLOTNAME)

    #  elif data['save'] is 3:
    #     counts = [len([x for x in history[-1] if x == i]) for i in range(B)]
    #     fig, ax = plt.subplots()
    #     plt.bar(range(len(counts)), counts)
    #     ax.set_xscale('log')
    #     plt.show()
