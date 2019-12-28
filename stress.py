import os
import random
import numpy as np
import time
from tqdm import tqdm

nodeNum = 8

P = 0.2
np.random.seed(0)
p = np.array([P, 1-P])
p = p.ravel()

T = 60


S =     ["1",
         "2",
         "3",
         "4",
         "5",
         "6"]

nodeList = ["cd /data/ost0;",
            "cd /data/ost1;",
            "cd /data/ost2;",
            "cd /data/ost3;",
            "cd /data/ost4;",
            "cd /data/ost5;",
            "cd /data/ost6;",
            "cd /data/ost7;"]

numIterations = 100;

print("stress generator start...\n")

for it in tqdm(range(numIterations)):
    sleepT = 0
    for i in range(nodeNum):
        is_straggler = random.random()
        if is_straggler>=0 and is_straggler<P:
            directory = nodeList[i]
            hdd = S[random.randint(0, len(S)-1)]
            stressTime = random.randint(0, T)
            stressCmd = "nohup stress --hdd " + hdd +" --timeout "+ stressTime + "&"
            cmd = directory + " " + stressCmd
            print("cmd is " + cmd)
            print("ost" + str(i) + " is straggler, degree is " + hdd +
                  ", time is " + str(stressTime))
            sleepT = max(sleepT, stressTime)
            os.system(cmd)
    time.sleep(sleepT)