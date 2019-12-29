import os
import random
import numpy as np
import time
from tqdm import tqdm
import sys

nodeNum = 8


S =     ["2",
         "3",
         "4",
         "5",
         "6",
         "7"]

nodeList = ["cd /data/ost0;",
            "cd /data/ost1;",
            "cd /data/ost2;",
            "cd /data/ost3;",
            "cd /data/ost4;",
            "cd /data/ost5;",
            "cd /data/ost6;",
            "cd /data/ost7;"]

print("stress generator start...\n")

#for it in tqdm(range(numIterations)):
def stress_generator(p, t):
    print("stress generator with p="+str(p)+" t="+str(t))
    while True:
        sleepT = 0
        for i in range(nodeNum):
            is_straggler = random.random()
            if is_straggler >= 0 and is_straggler < p:
                directory = nodeList[i]
                hdd = S[random.randint(0, len(S) - 1)]
                stressTime = random.randint(0, t)
                stressCmd = "nohup stress --hdd " + hdd + " --timeout " + str(
                    stressTime) + "&"
                cmd = directory + " " + stressCmd
                #print("cmd is " + cmd)
                print("ost" + str(i) + " is straggler, degree is " + hdd +
                      ", time is " + str(stressTime))
                sleepT = max(sleepT, stressTime)
                os.system(cmd)
        time.sleep(sleepT)


stress_generator(float(sys.argv[1]), int(sys.argv[2]))
