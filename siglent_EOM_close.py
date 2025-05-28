# import socket for communication with SIGLENT waveform generator
import socket
import sys
import time

#******************************************************************#
#******************************************************************#
# Interface with Siglent 
SIGLENT_IP = "10.2.0.11" # should match the instrument’s IP address
SIGLENT_PORT = 5024 # the port number of the instrument service
 
def SocketConnect():
    try:
        #create an AF_INET, STREAM socket (TCP)
        s_SIGLENT = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    except socket.error:
        print ('Failed to create socket.')
        sys.exit();
    try:
        #Connect to remote server
        s_SIGLENT.connect((SIGLENT_IP , SIGLENT_PORT))
    except socket.error:
        print ('Failed to connect to ip ' + SIGLENT_IP)
    return s_SIGLENT
	
def SocketSend(Sock, cmd):
    try :
        #Send cmd string
        Sock.sendall(cmd)
        Sock.sendall(b'\n')
        time.sleep(1)
    except socket.error:
        #Send failed
        print ('Send failed')
        sys.exit()
    #reply = Sock.recv(4096)
    #return reply
 
def SocketClose(Sock):
    #close the socket
    Sock.close()
    time.sleep(1)

#******************************************************************#
#******************************************************************#

# Initializes DC voltage setting 
def initialize_siglent_waveform_generator():

    s_SIGLENT = SocketConnect()
    qStr = SocketSend(s_SIGLENT, b'C1:BSWV WVTP,DC') #Set CH1 Wavetype to DC
    qStr = SocketSend(s_SIGLENT, b'C1:BSWV OFST,'+str.encode(str(0.5))) #Set CH1 DC Voltage
    qStr = SocketSend(s_SIGLENT, b'C1:OUTP LOAD,50') #Switch ON 
    qStr = SocketSend(s_SIGLENT, b'C1:OUTP ON') #Switch ON 

    return s_SIGLENT

# print("Initializing SIGLENT")
SocketClose(s_SIGLENT) #Close socket