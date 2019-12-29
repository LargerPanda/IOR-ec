import os
import random
import numpy as np
import time
from tqdm import tqdm
import subprocess

nodeNum = 8

P = 0.2

T = 60

S = ["1", "2", "3", "4", "5", "6"]

nodeList = [
    "cd /data/ost0;", "cd /data/ost1;", "cd /data/ost2;", "cd /data/ost3;",
    "cd /data/ost4;", "cd /data/ost5;", "cd /data/ost6;", "cd /data/ost7;"
]

numIterations = 100

result = []

for it in tqdm(range(numIterations)):
    print("stress generator start...\n")
    sleepT = 0
    for i in range(nodeNum):
        is_straggler = random.random()
        if is_straggler >= 0 and is_straggler < P:
            directory = nodeList[i]
            hdd = S[random.randint(0, len(S) - 1)]
            stressTime = random.randint(0, T) + 60
            stressCmd = "nohup stress --hdd " + hdd + " --timeout " + str(
                stressTime) + "&"
            cmd = directory + " " + stressCmd
            #print("cmd is " + cmd)
            print("ost" + str(i) + " is straggler, degree is " + hdd +
                  ", time is " + str(stressTime))
            sleepT = max(sleepT, stressTime)
            os.system(cmd)
    #run task
    print("task start...\n")
    f = subprocess.Popen("mpiexec -n 1 src/C/IOR -f read_8n_3g", shell=True, stdout=subprocess.PIPE)
    f.wait()
    print("process result...\n")
    lines = f.stdout.readlines()
    temp = []
    for line in lines:
        if line[0:4] == "read":
            strList = line.split()
            temp.append(float(strList[6]))

    mid = np.median(temp)
    largest = max(temp);
    variance = (largest-mid)/mid
    result.append(variance)
    print(result)
    time.sleep(sleepT)
    print("------------------------.\n")
    print("------------------------.\n")