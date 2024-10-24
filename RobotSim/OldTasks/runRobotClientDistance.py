import socket
import interfacesDist as interfaces
import numpy as np

sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
server_address = ('localhost', 5006)
print('starting up on {} port {}'.format(*server_address))
sock.bind(server_address)

interface = interfaces.DiscreteActionsRobot()
interface.mode = 4
interface.angle = 0
interface.debugLines = 0
interface.open()

robot_open = 1
target_pos = np.array([0.0, 300.0, 0.0])

while True:
	data, address = sock.recvfrom(4096)
	command = data[0]
	val1 = data[1]
	val2 = data[2]

	if len(data) > 3:
		val3 = data[3]
	else:
		val3 = 128

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
		if val1 == 3:		# Change hold time on debug lines
			interface.updateMode(val2);
		if val1 == 4:		# Change hold time on debug lines
			interface.updateDebugLines(val2);
		if val1 == 5:
			interface.updateDistanceDec(key, 0)
	if command == 1:	# Set Target
		target_pos[0] = (val1 - 128) / 100
		target_pos[1] = (val2 - 128) / 100
		target_pos[2] = (val3 - 128) / 100
		print(target_pos)
		interface.create_target3D(target_pos,1 )
		print(val1, val2, val3)
		print(target_pos)
	if command == 2:	# Set Dirr
		key = val1
		interface.updateDistanceDec(key, 1)
		interface.render()

sock.shutdown()
sock.close() 
