import sys
import time
import multiprocessing


def write_to_logfile(logfilename, *args):
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    with open(logfilename, 'a') as f:
        t = (current_time,) + args
        f.write(" ".join([str(s) for s in t]) + '\n')


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


###############################################################################


if __name__ == "__main__":
    pass

