#!/usr/bin/env python

import numpy as np

disc = 0.5
block_x = np.arange(10, 70+disc, disc)
block_y = np.arange(-60, 45+disc, disc)
block_z = np.arange(6370, 6300, -1)

np.savetxt("block_x", block_x, fmt="%.2f", header="1\n"+str(len(block_x)), comments='')
np.savetxt("block_y", block_y, fmt="%.2f", header="1\n"+str(len(block_y)), comments='')
np.savetxt("block_z", block_z, fmt="%.2f", header="1\n"+str(len(block_z)), comments='')
