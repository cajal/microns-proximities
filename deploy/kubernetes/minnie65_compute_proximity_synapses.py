if __name__ == '__main__':
    import time
    import numpy as np
    time.sleep(np.random.randint(5000)) # for starting 100 jobs simultaneously (to avoid database overload)
    from microns_proximities.minnie_proximities import minnie_proximities as mp
    mp.ProximitySynapse.populate(reserve_jobs=True, order='random', suppress_errors=True)

