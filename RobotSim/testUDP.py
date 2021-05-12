# import socket

# UDP_IP = "127.0.0.1"
# UDP_PORT = 3030
# MESSAGE = b"Hello, World!"

# print("UDP target IP: %s" % UDP_IP)
# print("UDP target port: %s" % UDP_PORT)
# print("message: %s" % MESSAGE)

# sock = socket.socket(socket.AF_INET, # Internet
#                      socket.SOCK_DGRAM) # UDP
# sock.sendto(MESSAGE, (UDP_IP, UDP_PORT))


import socket
import numpy as np
import struct

sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
# server_address = ('localhost', 5006)
# print('starting up on {} port {}'.format(*server_address))
# # sock.bind(server_address)


message = struct.pack("B", 1)
sock.connect(('127.0.0.1', 43209))
# message = 'TEST'
print(message)

sent = sock.send(message)
# sent = sock.sendto(message, server_address)
print("Sent")

sock.close()