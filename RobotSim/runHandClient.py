import socket
import interfaces
import numpy as np
from HandEnv import HandEnv

sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
server_address = ('localhost', 5006)
print('starting up on {} port {}'.format(*server_address))
sock.bind(server_address)


env = HandEnv()


while True:
	data, address = sock.recvfrom(4096)

	command = data[0]
	val1 = data[1]
	val2 = data[2]

	if command == 0:
		if val1 == 0:		# Open Robot 
			if robot_open == 0:
				interface.open()
				robot_open = 1
		if val1 == 1:		# Reset 
			env.reset()
	
	if command == 1:	# Set Target
		target_pos[0] = ((val1-1) *val2 + val3/100)/ 1250
		target_pos[1] = -((val4-1) *val5 + val6/100)/ 1250
		target_pos[2] = ((val7-1) *val8 + val9/100)/ 1250

		# interface.create_target3D(target_pos,1 )
		print(val1, val2, val3, val4, val5, val6)
		print(target_pos)
	if command == 2:	# Set Dirr
		key = val1
		interface.update_joystick(key)
		interface.render()
	if command == 3:
		key = val1
		interface.update_joystick(key)
		interface.render()
	if command == 4:	# Set Target
		key = val10
		robot_pos[0] = ((val1-1) *val2 + val3/100)/ 1250
		robot_pos[1] = -((val4-1) *val5 + val6/100)/ 1250
		robot_pos[2] = ((val7-1) *val8 + val9/100)/ 1250

		interface.updateRobotPos(robot_pos,key )
		interface.render()

sock.shutdown()
sock.close() 
