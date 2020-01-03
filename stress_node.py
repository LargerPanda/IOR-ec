import os
import random
import numpy as np
import time
from tqdm import tqdm
import sys

nodeNum = 8

S = ["2", "3", "4", "5", "6", "7"]

nodeList = [
    "cd /data/ost0;", "cd /data/ost1;", "cd /data/ost2;", "cd /data/ost3;",
    "cd /data/ost4;", "cd /data/ost5;", "cd /data/ost6;", "cd /data/ost7;"
]

print("stress generator start...\n")


#for it in tqdm(range(numIterations)):
def stress_generator(n, s, t):
    print("stress generator with p=" + str(p) + " t=" + str(t)+" n=" + str(n))
    while True:
        sleepT = 0
        i = n      
        directory = nodeList[i]
        hdd = S[s]
        stressTime = t
        stressCmd = "nohup stress --hdd " + hdd + " --timeout " + str(
            stressTime) + "&"
        cmd = directory + " " + stressCmd
        #print("cmd is " + cmd)
        print("ost" + str(i) + " is straggler, degree is " + hdd +
                ", time is " + str(stressTime))
        sleepT = max(sleepT, stressTime)
        os.system(cmd)
        time.sleep(sleepT)

###stress_node.py ost_num s t
###ost_num : 0-7
###s       : 0-5
###t       : 5 10 20 40 80
stress_generator(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
