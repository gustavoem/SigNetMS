from pathos.multiprocessing import ProcessPool

# Formerly, we used a solution that was inspired in this thread:
# https://stackoverflow.com/questions/3288595/multiprocessing-how-to-use-pool-map-on-a-function-defined-in-a-class

def parallel_map (f, X, nof_process):
    """ Runs a map of X to f in parallel. 
    
    Parameters
        f: a callable object to which X should be applied.
        X: a list of objects to which apply f.
        nof_process: a integer representing the number of process that 
            should be spawned to calculate f(X)

    Returns
        results: a list containing the resulting application of f to 
            every element in X.
    """
    pool = ProcessPool (nof_process)
    pool.clear () # why would you clear something you just created?
    # anser: pathos.multiprocessing caches worker pools. (yeahhh... 
    # maybe this should be a class) see issue #49
    results = pool.map (f, X)
    return results
