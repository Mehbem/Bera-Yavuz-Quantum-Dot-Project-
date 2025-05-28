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

def siglent_DC_offset(s_SIGLENT, DC_voltage_val):
    print('Siglent DC Amp. Set')
    qStr = SocketSend(s_SIGLENT, b'C1:BSWV OFST,'+str.encode(str(DC_voltage_val))) #Set CH1 DC Voltage
    time.sleep(0.5)

set_DC_val = siglent_DC_offset(s_SIGLENT, DC_voltage_val)

# print("Setting DC voltage value")
# DC_voltage_val = 0.24
# siglent_DC_offset(s_SIGLENT, DC_voltage_val)