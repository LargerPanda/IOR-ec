import os
import random
import numpy as np
import time
from tqdm import tqdm

os.system("sudo ./refresh.sh")
os.system("mpiexec -n 1 src/C/IOR -f read_8n_3g")
os.system("sudo ./refresh.sh")
os.system("mpiexec -n 1 src/C/IOR -f read_8n_3g_ec")