#!/usr/bin/env python
"""
coverage -- explore the use of the 'coverage' measure for predicting which
            combination of interst operators should work best
"""
import os, sys, math, numpy, itertools
from mpi4py import MPI

# CLUSTER PREPARATION
# Before doing anything, get info about the cluster and your
# role in it.
comm  = MPI.COMM_WORLD
rank  = comm.Get_rank()
procs = comm.Get_size()

def load_locations (fn):
    fd = open (fn)
    lines = fd.readlines()
    fd.close()
    nf = int (lines[1])
    locs = numpy.zeros ((nf, 2))
    for f in range (0, nf):
        v = lines[f+2].split()
        # Each line contains x, y, and stuff we're not interested in.
        locs[f,1] = float (v[0])
        locs[f,0] = float (v[1])
    return locs

def coverage (locs):
    nlocs, loclen = locs.shape
    if loclen != 2: raise ValueError, "I don't have two values for the location!"
    dsum = 0.0
    npaths = 0
    ncoin = 0
    for i in range (0, nlocs):
        y1 = locs[i,0]
        x1 = locs[i,1]
        for j in range (i+1, nlocs):
            y2 = locs[j,0]
            x2 = locs[j,1]
            d = math.sqrt ((y2 - y1)**2 + (x2 - x1)**2)
            if d == 0.0:
                ncoin += 1
            else:
                dsum += 1.0 / d
                npaths += 1
    return npaths / dsum, npaths, ncoin

def distribute_coverage (locs):
    '''
    Function called by the root node (0), that distributes the calculation
    of coverage components (enumeration of combinations) and gathers the
    results, returning the coverage for all data.
    '''
    nlocs, loclen = locs.shape
    if loclen != 2: raise ValueError, "I don't have two values for the location!"
    dsum = 0.0
    npaths = 0
    ncoin = 0

    # DISTRIBUTE CALCULATION
    # broadcast the locs
    #comm.Bcast([locs, MPI.FLOAT], root=0)

    # calculate my part of the calculation
    partial_results = calc_my_partial_coverage(locs)

    comm.Barrier()
    # reduce the calculation
    results = numpy.array([0.0,0.0,0.0])
    comm.Reduce(partial_results, results, op.MPI.SUM, root=0)

    if rank == 0:
        # add my partial results to the reduction
        results = numpy.add(partial_results, results)


    # return the results
    dsum   = results[0]
    npaths = results[1]
    ncoin  = results[2]

    return npaths / dsum, npaths, ncoin

def calc_my_partial_coverage(locs):
    '''
    Calculates the partial coverage components (dsum, npaths, ncoin)
    based on the rank of the processor and the data that every other
    processor has access to.
    '''
    nlocs, loclen = locs.shape

    # get the length of sections (1 section of the data per processor)
    # therefore, we have the same amount of sections as processors
    sec_len= nlocs/procs

    # initialize partial coverage components placeholders
    dsum = 0.0
    npaths = 0
    ncoin = 0

    # COVERAGE OF MY SECTION
    lowerbound, upperbound = get_bounds(procs, rank, sec_len, nlocs)
    for i in range (lowerbound, upperbound):
        y1 = locs[i,0]
        x1 = locs[i,1]
        for j in range (i+1, upperbound):
            y2 = locs[j,0]
            x2 = locs[j,1]
            d = math.sqrt ((y2 - y1)**2 + (x2 - x1)**2)
            if d == 0.0:
                ncoin += 1
            else:
                dsum += 1.0 / d
                npaths += 1

    # COVERAGE OF COMBINATIONS OF SECTIONS
    # enumerate combinations of sections
    sec_combs = list(itertools.combinations(range(procs), 2))

    # get number of combinations to calculate
    nsec_comb = len(sec_combs)

    # get number of combinations per processor
    sec_comb_len = nsec_comb / procs
    sec_comb_len = 1 if sec_comb_len == 0 else sec_comb_len

    # calculate coverage for my combinations of sections
    comb_lower_b, comb_upper_b = get_bounds(sec_comb_len, nsec_comb)

    # make sure the master node does nothing if there are only two processors
    if procs == 2 and rank == 0:
        comb_lower_b = comb_upper_b
    elif procs == 2 and rank == 1:
        comb_lower_b = 0
        comb_upper_b = 2

    for ic in range(comb_lower_b, comb_upper_b):
        combination = sec_combs[ic]
        # get boundaries for sections to combine
        lowerbound1, upperbound1 = get_bounds(procs, combination[0], sec_len, nlocs)
        lowerbound2, upperbound2 = get_bounds(procs, combination[1], sec_len, nlocs)
        for i in range (lowerbound1, upperbound1):
            y1 = locs[i,0]
            x1 = locs[i,1]
            for j in range (lowerbound2, upperbound2):
                y2 = locs[j,0]
                x2 = locs[j,1]
                d = math.sqrt ((y2 - y1)**2 + (x2 - x1)**2)
                if d == 0.0:
                    ncoin += 1
                else:
                    dsum += 1.0 / d
                    npaths += 1

    return numpy.array([dsum, npaths, ncoin])

def get_bounds(section_length, data_length):
    '''
    Returns lower and upper bounds in data, like so:
     [lowerbound, upperbound)
    So upperbound is non-inclusive.

    Data was broken into 1 section per processor, and a
    section watch assigned to each processor contiguously, so:
     rank 0 has first section, rank 1 has second section, and so on.
    '''
    lowerbound = rank * section_length
    upperbound = data_length if (rank >= (procs-1)) else (rank+1)*section_length
    return lowerbound, upperbound

def mean_coverage1 (opname, imsets):
    csum = 0.0
    nfiles = 0
    for imset in imsets:
        for fno in range (1,7):
            fn = 'locations/%s_%s%d.txt' % (opname, imset, fno)
            locs = load_locations (fn)
            c, np, nc = distribute_coverage (locs)
            print '  ', fn, 'gave', np, 'paths and', nc, 'coincident points.'
            csum += c
            nfiles += 1
    return csum / nfiles

def mean_coverage2 (op1, op2, imsets):
    csum = 0.0
    nfiles = 0
    for imset in imsets:
        for fno in range (1,7):
            fn1 = 'locations/%s_%s%d.txt' % (op1, imset, fno)
            locs1 = load_locations (fn1)
            fn2 = 'locations/%s_%s%d.txt' % (op2, imset, fno)
            locs2 = load_locations (fn2)
            nlocs1 = len (locs1)
            nlocs = nlocs1 + len (locs2)
            locs = numpy.zeros ((nlocs, 2))
            locs[0:nlocs1,:] = locs1
            locs[nlocs1:,:] = locs2
            c, np, nc = distribute_coverage (locs)
            print '  ', fn1, fn2, 'gave', np, 'paths and', \
                  nc, 'coincident points.'
            csum += c
            nfiles += 1
    return csum / nfiles

def mean_coverage3 (op1, op2, op3, imsets):
    csum = 0.0
    nfiles = 0
    for imset in imsets:
        for fno in range (1,7):
            fn1 = 'locations/%s_%s%d.txt' % (op1, imset, fno)
            locs1 = load_locations (fn1)
            fn2 = 'locations/%s_%s%d.txt' % (op2, imset, fno)
            locs2 = load_locations (fn2)
            fn3 = 'locations/%s_%s%d.txt' % (op3, imset, fno)
            locs3 = load_locations (fn3)
            nlocs1 = len (locs1)
            nlocs = nlocs1 + len (locs2) + len(locs3)
            locs = numpy.zeros ((nlocs, 2))
            locs[0:nlocs1,:] = locs1
            locs[nlocs1:nlocs1+len(locs2),:] = locs2
            locs[nlocs1+len(locs2):,:] = locs3
            c, np, nc = distribute_coverage (locs)
            print '  ', fn1, fn2, fn3, 'gave', np, 'paths and', \
                  nc, 'coincident points.'
            csum += c
            nfiles += 1
    return csum / nfiles

def mean_coverage4 (op1, op2, op3, op4, imsets):
    csum = 0.0
    nfiles = 0
    for imset in imsets:
        for fno in range (1,7):
            fn1 = 'locations/%s_%s%d.txt' % (op1, imset, fno)
            locs1 = load_locations (fn1)
            fn2 = 'locations/%s_%s%d.txt' % (op2, imset, fno)
            locs2 = load_locations (fn2)
            fn3 = 'locations/%s_%s%d.txt' % (op3, imset, fno)
            locs3 = load_locations (fn3)
            fn4 = 'locations/%s_%s%d.txt' % (op4, imset, fno)
            locs4 = load_locations (fn4)
            nlocs1 = len (locs1)
            nlocs = nlocs1 + len (locs2) + len(locs3) + len(locs4)
            locs = numpy.zeros ((nlocs, 2))
            locs[0:nlocs1,:] = locs1
            locs[nlocs1:nlocs1+len(locs2),:] = locs2
            locs[nlocs1+len(locs2):nlocs1+len(locs2)+len(locs3),:] = locs3
            locs[nlocs1+len(locs2)+len(locs3):,:] = locs4
            c, np, nc = distribute_coverage (locs)
            print '  ', fn1, fn2, fn3, fn4, 'gave', np, 'paths and', \
                  nc, 'coincident points.'
            csum += c
            nfiles += 1
    return csum / nfiles

def sort_by_value (d):
    '''Return the keys of dictionary d sorted by their values'''
    items = d.items()
    backitems = [[v[1],v[0]] for v in items]
    backitems.sort()
    return [backitems[i][1] for i in range(0,len(backitems))]

#------------------------------------------------------------------------------
# Main program
#------------------------------------------------------------------------------

opnames = ['ebr', 'ibr', 'mser', 'sfop']
imsets = ['bark', 'bikes', 'boat', 'graf', 'leuv', 'trees', 'ubc', 'wall']
results = {}

# Calculate the coverage for each combination of quadruples of operators.
nops = len (opnames)
for i1 in range (0, nops):
    op1 = opnames[i1]
    for i2 in range (i1+1, nops):
        op2 = opnames[i2]
        for i3 in range (i2+1, nops):
            op3 = opnames[i3]
            for i4 in range (i3+1, nops):
                op4 = opnames[i4]
                mc = mean_coverage4 (op1, op2, op3, op4, imsets)
                print op1, op2, op3, op4, mc
                results[op1 + ' + ' + op2 + ' + ' + op3 + ' + ' + op4] = mc

# Calculate the coverage for each combination of triples of operators.
nops = len (opnames)
for i1 in range (0, nops):
    op1 = opnames[i1]
    for i2 in range (i1+1, nops):
        op2 = opnames[i2]
        for i3 in range (i2+1, nops):
            op3 = opnames[i3]
            mc = mean_coverage3 (op1, op2, op3, imsets)
            print op1, op2, op3, mc
            results[op1 + ' + ' + op2 + ' + ' + op3] = mc

# Calculate the coverage for each combination of pairs of operators.
nops = len (opnames)
for i1 in range (0, nops):
    op1 = opnames[i1]
    for i2 in range (i1+1, nops):
        op2 = opnames[i2]
        mc = mean_coverage2 (op1, op2, imsets)
        print op1, op2, mc
        results[op1 + ' + ' + op2] = mc

# Calculate the coverage for each individual operator.
for op in opnames:
    mc = mean_coverage1 (op, imsets)
    print op, mc
    results[op] = mc

# Output what we've calculated in descending order of coverage.
ds = sort_by_value (results)
for k in reversed (ds):
    print k, results[k]

#------------------------------------------------------------------------------
# End of coverage
#------------------------------------------------------------------------------
