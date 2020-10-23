import socket
import interfaces
import numpy as np

sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
server_address = ('localhost', 5006)
print('starting up on {} port {}'.format(*server_address))
sock.bind(server_address)

interface = interfaces.DiscreteActionsRobot()
interface.open()
robot_open = 1
target_pos = np.array([0.0, 300.0])

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
			interface.reset()
		if val1 == 2:		# Change hold time on debug lines
			updateRate = 1.0 /val2;
			interface.updateRefresh(updateRate);
	if command == 1:	# Set Target
		# interface.reset()
		target_pos[0] = (val1 - 128) / 100
		target_pos[1] = (val2 - 128) / 100
		interface.create_target(target_pos)
		print(val1, val2)
		print(target_pos)
	if command == 2:	# Set Dirr
		key = val1
		interface.update_joystick(key)
		interface.render()

sock.shutdown()
sock.close() 