import os
import random
import numpy as np
import time
from tqdm import tqdm
import multiprocessing
import subprocess


P = [0.1, 0.15, 0.2, 0.25, 0.3]
T = [10, 20, 40, 80, 100, 120]
S = ["2", "3", "4", "5", "6", "7"]

nodeNum = 8
nodeList = [
    "cd /data/ost0;", "cd /data/ost1;", "cd /data/ost2;", "cd /data/ost3;",
    "cd /data/ost4;", "cd /data/ost5;", "cd /data/ost6;", "cd /data/ost7;"
]

fileSizeList = ["64MB", "256MB", "1GB", "4GB", "16GB", "64GB"]
scalesizeList = [2, 4, 8]

fd = open("data.txt", 'w')

def stress_generator(p, t):
    print("stress generator with p=%f, t=%d start...", p, t)
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

print("test start...")

def test_straggler(numIterations, fileSize, scaleSize):
    #print("-----------------------------------------------------.")
    result = []
    generate_script = "generate_"+fileSize+"_"+str(scaleSize)
    generate_cmd = "mpiexec -n 1 src/C/IOR -f " + generate_script
    print(generate_cmd)
    read_script = "read_" + fileSize + "_" + str(scaleSize)
    read_cmd = "mpiexec -n 1 src/C/IOR -f " + read_script
    print(read_cmd)
    f = subprocess.Popen(generate_cmd, shell=True, stdout=subprocess.PIPE)
    f.wait()
    for it in tqdm(range(numIterations)):
        #run task
        os.system("sudo ./refresh.sh")
        print("refreshing...")
        print("task start...")
        f = subprocess.Popen(read_cmd,
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
    fd.write(str(result)+'\n')
    fd.flush()
    print(result)
    #print("-----------------------------------------------------.")

for p in P:
    for t in T:
        stress_process = multiprocessing.Process(target=stress_generator, args=(p, t))
        stress_process.start()
        for fileSize in fileSizeList:
            for scaleSize in scalesizeList:
                test_straggler(1, fileSize, scaleSize)
        stress_process.terminate()
        stress_process.join()

fd.close()