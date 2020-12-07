#!/usr/bin/env python3

import sys
import random

NUM_VARS=int(sys.argv[1])
NUM_SAMPLES=int(sys.argv[2])


for i in range(0, NUM_SAMPLES):
    sample = [str(i)]
    sample += [str(random.randint(0,1))]
    for i in range(NUM_VARS):
        sample += [str(random.randint(0,1))]
    print(",".join(sample))
