#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
#import matplotlib
import numpy as np
#import scipy
#import math


x = np.linspace(0.0, 5.0, 1000)
y = [
        np.exp(0.9 * x),
        np.exp(1.0 * x),
        np.exp(1.1 * x),
        np.exp(1.2 * x)
    ]


fig = plt.figure(figsize=(6, 3))
for i in xrange(len(y)):
    plt.plot(x, y[i], '-r', lw=3)

plt.show()
