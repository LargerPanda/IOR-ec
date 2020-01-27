import os
import random
import numpy as np
import time
from tqdm import tqdm
from random import choice
import sys

nodeNum = 8

S = ["2", "3", "4", "5", "6", "7", "8"]

nodeList = [
    "cd /data/ost0;", "cd /data/ost1;", "cd /data/ost2;", "cd /data/ost3;",
    "cd /data/ost4;", "cd /data/ost5;", "cd /data/ost6;", "cd /data/ost7;"
]

print("stress generator start...\n")


#for it in tqdm(range(numIterations)):
def stress_generator():
    while True:
        sleepT = 0
        directory = nodeList[0]
        hdd = str(choice([4,5,6,7,8]))
        stressTime = choice([3,4,5,6,7,8])
        stressCmd = "nohup stress --hdd " + hdd + " --timeout " + str(
            stressTime) + "&"
        cmd = directory + " " + stressCmd
        #print("cmd is " + cmd)
        print(directory + " is straggler, degree is " + hdd +
              ", time is " + str(stressTime))
        
        sleepT = max(sleepT, stressTime) + choice([2,3,4,5])
        os.system(cmd)
        time.sleep(sleepT)


###stress_node.py ost_num s t
###ost_num : 0-7
###s       : 0-5
###t       : 5 10 20 40 80
stress_generator()
