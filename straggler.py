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



def test_straggler(numIterations, script):
    print("-----------------------------------------------------.")
    result = []
    cmd = "mpiexec -n 1 src/C/IOR -f " + script
    for it in tqdm(range(numIterations)):
        #run task
        os.system("sudo ./refresh.sh")
        print("refreshing...")
        print("task start...")
        f = subprocess.Popen(cmd,
                         shell=True,
                         stdout=subprocess.PIPE)
        f.wait()
        print("process result...\n")
        lines = f.stdout.readlines()
        temp = []
        for line in lines:
            line = line.decode('utf-8')
            if line[0:5] == "#read":
                strList = line.split()
                #print(float(strList[6]))
                temp.append(float(strList[6]))

        mid = np.median(temp)
        #print("mid=%f", mid)
        largest = max(temp)
        #print("largest=%f", largest)
        variance = (largest - mid) / mid
        result.append(variance)
        #print(result)
        #time.sleep(sleepT)
    print(result)
    print("-----------------------------------------------------.")
