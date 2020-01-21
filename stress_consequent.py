import os
import random
import numpy as np
import time
from tqdm import tqdm
from random import choice
import sys

nodeNum = 8

nodeList = [
    "cd /data/ost0;", "cd /data/ost1;", "cd /data/ost2;", "cd /data/ost3;"#,
    #"cd /data/ost4;", "cd /data/ost5;", "cd /data/ost6;", "cd /data/ost7;"
]

print("stress generator start...\n")


#for it in tqdm(range(numIterations)):
def stress_generator(n, s, t):
    print("stress generator with s=" + str(s) + " t=" + str(t) +
          " concurrent straggler =" + str(n))
    i = 0
    while True:
        sleepT = 0
        for each in range(n):
            #directory = choice(nodeList)
            directory = nodeList[(each+i)%len(nodeList)]
            hdd = str(s)
            stressTime = t
            stressCmd = "nohup stress --hdd " + hdd + " --timeout " + str(
                stressTime) + "&"
            cmd = directory + " " + stressCmd
            #print("cmd is " + cmd)
            print(directory + " is straggler, degree is " + hdd + ", time is " +
                str(stressTime))
            sleepT = max(sleepT, stressTime)
            os.system(cmd)
        time.sleep(sleepT)
        i = (i + n) % len(nodeList)


###stress_node.py num s t
###stragglers at the same time :
###s       : 0-6
###t       : 5 10 20 40 80
stress_generator(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
