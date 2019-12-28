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
    for i in range(nodeNum):
        is_straggler = np.random.choice([1.0, 0.0], p)
        if is_straggler:

            directory = nodeList[i]
            hdd = hddList[random.randint(0, len(hddList)-1)]
            time = str(random.randint(0, T))
            stressCmd = "nohup stress --hdd " + hdd +" --timeout "+ time + "&"
            cmd = directory + " " + stressCmd
            print("cmd is " + cmd)
            print("ost" + i + " is straggler, degree is " + hdd +
                  ", time is " + time)
            os.system(cmd)
    time.sleep(T)