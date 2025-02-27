# SNSPD Input Counter Read
# Python modules needed

import os
import time
import socket
import zmq
import itertools
import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd

def check_host(address, port):
    s = socket.socket()
    s.settimeout(5)
    try:
        s.connect((address, port))
        s.settimeout(None)
        return True
    except socket.error as e:
        return False

def error(message):
    print('Error: ' + message)
    exit(1)

def zmq_exec(zmq, cmd):
    print(cmd)
    zmq.send_string(cmd)
    ans = zmq.recv().decode('utf-8')
    if ans: print('>', ans)
    return ans

#################################################################
###################   INITIALIZE ID-900   #######################
#################################################################

def init_ID900():

    # IP address of the ID900 
    ID900_IP = "169.254.102.131"
    ID900_PORT = 5555
    ID900_ADDR = 'tcp://' + ID900_IP + ':' + str(ID900_PORT)

    # Check if ID900 is listening
    if not check_host(ID900_IP, ID900_PORT):
        error('Unable to connect to ID900 "' + ID900_IP + '" on port ' + str(ID900_PORT) + '.')

    # Create zmq socket and connect to the ScpiClient
    context = zmq.Context()
    tc = context.socket(zmq.REQ)
    tc.connect(ID900_ADDR)

    zmq_exec(tc, "INPUt4:COUP DC")
    zmq_exec(tc, "INPUt4:EDGE FALLING")
    zmq_exec(tc, "INPUt4:THRE 1V")
    zmq_exec(tc, "INPUt4:SELE UNSHAPED")
    zmq_exec(tc, "INPUt4:ENAB ON")
    zmq_exec(tc, "INPUt4:Counter:MODE CYCLE")
    zmq_exec(tc, "INPUt4:Counter:INTEgrationtime 100")
    
    return tc

#################################################################
################   QUERY PHOTON COUNTER   #######################
#################################################################

def query_photon_counter(tc):

    countersValue = 0

    print("STARTING COLLECTION")

    n_repeats = 3; 

    for i in range(n_repeats):

        time.sleep(0.5)
        countersValue = countersValue + np.float64(zmq_exec(tc, "INPUt2:COUNter?"))
        #print(countersValue)
        
    
    countersValue = countersValue/n_repeats

    return countersValue



