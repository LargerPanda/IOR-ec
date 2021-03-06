import os
import random
import numpy as np
import time
from tqdm import tqdm
import multiprocessing
import subprocess


P = [0.1, 0.2, 0.3]
T = [5, 10, 20, 40, 80]
S = ["2", "3", "4", "5", "6", "7"]

nodeNum = 8
nodeList = [
    "cd /data/ost0;", "cd /data/ost1;", "cd /data/ost2;", "cd /data/ost3;",
    "cd /data/ost4;", "cd /data/ost5;", "cd /data/ost6;", "cd /data/ost7;"
]

fileSizeList = ["64MB", "256MB", "1GB", "4GB", "16GB"]
scalesizeList = [2, 4, 8]

fd = open("data.txt", 'a+')

def stress_generator(p, t):
    print("stress generator with p=" + str(p) + " t=" + str(t))
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
        variance = largest / mid
        result.append(variance)
        #print(result)
        #time.sleep(sleepT)
    fd.write(str(result)+'\n')
    fd.flush()
    print(result)
    #print("-----------------------------------------------------.")

countP = 0
countT = 0
countF = 0
countS = 0

startP = 2
startT = 0
startF = 4
startS = 1

for p in P:
    if countP>=startP:
        for t in T:
            if countT >= startT:
                stress_process = multiprocessing.Process(target=stress_generator, args=(p, t))
                stress_process.start()
                for fileSize in fileSizeList:
                    if countF>=startF:
                        for scaleSize in scalesizeList:
                            if countS >= startS:
                                print(
                                    str(countP) + ' ' + str(countT) + ' ' +
                                    str(countF) + ' ' + str(countS))
                                fd.write(
                                    str(countP) + ' '+str(countT) +' '+ str(countF) + ' ' +
                                    str(countS) + '\n')
                                fd.flush()
                                test_straggler(5, fileSize, scaleSize)
                            countS = countS + 1
                        countS = 0
                        startS = 0
                    countF = countF+1
                countF=0
                startF = 0
            if countT >= startT:
                time.sleep(10)
                stress_process.terminate()
                stress_process.join()
            countT = countT+1
        countT = 0
        startT = 0
    countP = countP+1
startP = 0
countP = 0

fd.close()