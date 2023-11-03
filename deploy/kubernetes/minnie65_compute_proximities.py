if __name__ == '__main__':
    from microns_proximities.minnie_proximities import minnie_proximities as mp
    mp.ProximityMaker.SkeletonProcessedSetChunk.populate(reserve_jobs=True, order='random', suppress_errors=True)

