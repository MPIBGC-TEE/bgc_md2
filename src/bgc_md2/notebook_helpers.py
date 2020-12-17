import sys
import time
import multiprocessing

import shutil
import zarr

import numpy as np


def write_to_logfile(logfilename, *args):
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    with open(logfilename, 'a') as f:
        t = (current_time,) + args
        f.write(" ".join([str(s) for s in t]) + '\n')


def write_header_to_logfile(logfile_name, var_da, time_limit_in_min):
#    c = var_da.chunks
#    nr_chunks = np.prod([len(val) for val in c])
#    print('nr_chunks:', nr_chunks)
    nr_singles = np.prod(var_da.shape[:3])
    write_to_logfile(
        logfile_name,
        'starting',
#        nr_chunks, "chunks, ",
        nr_singles, "singles,",
        "timeout (min) =",
        time_limit_in_min
    )
    
    s = "nr_singles: " + str(nr_singles) + "\n"
    s += "timeout (min): " + str(time_limit_in_min) + "\n"
    return s


def _custom_timeout_target(queue, function, *args, **kwargs):
    try:
        queue.put((True, function(*args, **kwargs)))
    except:
        queue.put((False, sys.exc_info()[1]))

def custom_timeout(seconds, function, *args, **kwargs):
                q = multiprocessing.Queue(1)
                args = (q, function) + args
                
                p = multiprocessing.Process(
                    target=_custom_timeout_target,
                    args=args,
                    kwargs=kwargs
                )
                p.daemon = True

                timeout = time.time() + seconds

                def cancel():
                    if p.is_alive():
                        p.terminate()

                def ready():
                    if timeout < time.time():
                        cancel()
                        raise(TimeoutError)

                    return q.full() and not q.empty()
                
                p.start()
                while not ready():
                    time.sleep(0.01)

                flag, result = q.get()
                if not flag:
                    raise result

                return result


def create_zarr_archive(zarr_path, result_shape, result_chunks, overwrite=False):
    if overwrite & zarr_path.exists():
        shutil.rmtree(zarr_path)

    z = zarr.create(
        result_shape,
        chunks=result_chunks,
        dtype=np.float64,
        fill_value=-np.inf, # -np.inf indicates incomplete computation
        store=str(zarr_path)
    )
    return z


###############################################################################


if __name__ == "__main__":
    pass

