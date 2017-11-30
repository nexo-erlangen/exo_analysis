#!/usr/bin/env python
import ternary
import matplotlib.pyplot as plt
import numpy as np

fig, tax = ternary.figure(scale=scale)
tax.heatmapf(:boundary=True, style='triangular')
tax.boundary(linewidth=2.)

tax.show()

