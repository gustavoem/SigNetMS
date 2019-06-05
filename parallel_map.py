import multiprocessing

# This solution is inspired on the solution of klaus se on the Stack 
# Overflow thread:
# https://stackoverflow.com/questions/3288595/multiprocessing-how-to-use-pool-map-on-a-function-defined-in-a-class

def __fun_runner (f, q_in, q_out):
    while True:
        i, x = q_in.get ()
        if i is None:
            break
        run_result = (i, f (x))
        q_out.put (run_result)


def parallel_map (f, X, nof_process):
    """ Runs a map of X to f in parallel, using nof_process process. """
    q_in = multiprocessing.Queue (1)
    q_out = multiprocessing.Queue ()

    proc = [multiprocessing.Process (target=__fun_runner, 
        args=(f, q_in, q_out)) for _ in range (nof_process)]

    for p in proc:
        p.daemon = True
        p.start ()
    
    sent = 0
    for i, x in enumerate (X):
        q_in.put ((i, x))
        sent += 1

    # sends a final signal to every worker
    for _ in range (nof_process):
        q_in.put ((None, None), block=True)
    q_in.close ()

    res = [q_out.get () for _ in range (sent)]

    for p in proc:
        p.join ()
    return [x for i, x in sorted (res)]
