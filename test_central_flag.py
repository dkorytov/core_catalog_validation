#!/usr/bin/env python

from __future__ import print_function

import dtk
import sys
import numpy as np
import matplotlib.pyplot as plt

def test(fname):
    dtk.gio_inspect(fname)
    centrals = dtk.gio_read(fname, 'central')
    infall_step = dtk.gio_read(fname, 'infall_step')
    print(centrals)
    print(np.unique(centrals))
    print(np.sum(centrals), len(centrals), np.float(np.sum(centrals))/len(centrals))


if __name__ == "__main__":
    test(sys.argv[1])

