import socket

my_socket= socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
my_socket.connect(('127.0.0.1', 8822))

MESSAGE='test1'
for i in range(1,10):
    my_socket.send(MESSAGE)
    print i

my_socket.close
