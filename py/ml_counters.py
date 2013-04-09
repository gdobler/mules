import numpy as np

class counters():

    """
      global counters and cell lists to keep track of the number of
      cells and individual stars used
    """

    cellcnt = 0
    starcnt = 0
    clist   = []
    mmrange = np.arange(100)
    glones  = np.ones(100000)
