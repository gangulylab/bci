import socket
import interfaces
import numpy as np
import math

sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
server_address = ('localhost', 5006)
# server_address = ('0.0.0.0', 5006)
print('starting up on {} port {}'.format(*server_address))
sock.bind(server_address)

interface = interfaces.DiscreteActionsRobot()
interface.mode = 3
interface.angle = 0
interface.debugLines = 0
interface.open()
interface.letterMode = 0

robot_open = 1
target_pos = np.array([0.0, 300.0, 0.0])
robot_pos = np.array([0.0, 0.0, 0.0])
robot_orn = np.array([0.0, 0.0, 0.0])
p =[0.0, 0.0, 0.0]
 

while True:
	data, address = sock.recvfrom(4096)

	command = data[0]
	val1 = data[1]
	val2 = data[2]

	if len(data) > 3:
		val3 = data[3]
	else:
		val3 = 128
	if len(data) > 4:
		val4 = data[4]
	if len(data) > 5:
		val5 = data[5]
		val6 = data[6]
		val7 = data[7]
		val8 = data[8]
		val9 = data[9]
		val10 = data[10]

	else:
		val4 = 0

	if command == 0:
		if val1 == 0:		# Open Robot 
			if robot_open == 0:
				interface.open()
				robot_open = 1
		if val1 == 1:		# Reset 
			interface.reset()
		if val1 == 2:		# Change hold time on debug lines
			updateRate = 1.0 /val2
			interface.updateRefresh(updateRate)
			interface.updateRefresh(.2)
		if val1 == 3:		# Mode
			interface.updateMode(val2)
			print("MODE:", val2)
		if val1 == 4:		# Use debug lines
			interface.updateDebugLines(val2)
			print("DEBUG LINES", val2)
		if val1 == 5:		# Target type
			if interface.letterMode == 1:
				interface.create_targetLetter(target_pos, 0)
			elif interface.mode == 9:
				interface.create_targetRing(target_pos,0)
			else:
				interface.create_target3D(target_pos, 0, target)

		if val1 == 6:
			if interface.letterMode == 1:
				interface.create_targetLetter(target_pos, 2)
			elif interface.mode == 9:
				interface.create_targetRing(target_pos,2)
			else:
				interface.create_target3D(target_pos, 2, target)
		if val1 == 7:
			if interface.letterMode == 1:
				interface.create_targetLetter(target_pos, 3)
			elif interface.mode == 9:
				interface.create_targetRing(target_pos,3)
			else:
				interface.create_target3D(target_pos,3, target)
		if val1 == 8:
			interface.setMode(val2)
		if val1 == 9:
			interface.grabCube()
		if val1 == 15:
			interface.update_color(target_pos, val2)
		if val1 == 16:
			interface.targetID = val2
		if val1 == 17:
			interface.letterMode = val2
		if val1 == 18:
			interface.setTargetRad(val2/1000.0)
		if val1 == 19:
			interface.setPath(val2)
			interface.drawPath()
		if val1 == 20:
			interface.setGoCue(val2)
		if val1 == 21:
			interface.drawPath()
		if val1 == 22:
			interface.robotenv.removeDebug()
		if val1 == 23:
			betaVal = ((val2-1) *(val3 + val4/100))
			print(betaVal)
			interface.robotenv.drawBetaLine(betaVal)
		if val1 == 24:
			startOrn = val2
			if startOrn == 1:
				interface.robotenv.set_robotOrn([math.pi, 0, math.pi])
			if startOrn == 2:
				interface.robotenv.set_robotOrn([0, -math.pi/2, 0])
			interface.render()


	if command == 1:	# Set Target
		target_pos[0] = ((val1-1) *(val2 + val3/100))/ 1000
		target_pos[1] = -((val4-1) *(val5 + val6/100))/ 1000
		target_pos[2] = ((val7-1) *(val8 + val9/100))/ 1000
		target = val10
		if interface.letterMode == 1:
			interface.create_targetLetter(target_pos, 1)
		elif interface.mode == 9:
			interface.create_targetRing(target_pos,1 )
		else:
			interface.create_target3D(target_pos, 1, target)


	if command == 11:	# Set Target
		target_pos[0] = ((val1-1) *(val2 + val3/100))/ 1000
		target_pos[1] = -((val4-1) *(val5 + val6/100))/ 1000
		target_pos[2] = ((val7-1) *(val8 + val9/100))/ 1000
		print(target_pos)
		interface.create_target(target_pos)

	if command == 4:	
		if interface.mode == 10 or interface.mode == 12:
			interface.robotenv.drawMode(1)
		key = val10
		robot_pos[0] = ((val1-1) *(val2 + val3/100))/ 1000
		robot_pos[1] = -((val4-1) *(val5 + val6/100))/ 1000
		robot_pos[2] = ((val7-1) *(val8 + val9/100))/ 1000
		interface.robotenv.opMode = 0
		interface.updateRobotPos(robot_pos,key )
		print(key)
		interface.render()

	
	if command == 7:	
		interface.robotenv.drawMode(2)
		key = val10
		
		robot_ornz = ((val1-1) *(val2 + val3/100))/10

		# interface.robotenv.opMode = 1
		# interface.robotenv.set_robotOrn([math.pi, 0, robot_ornz])

		interface.robotenv.set_robotOrn([robot_ornz,-math.pi/2,0])

		print(robot_ornz)

		grasp =  val5

		if grasp == 1:
			interface.robotenv.setFing(0)
		elif grasp == 2:
			interface.robotenv.setFing(1.3)

		interface.render()

	if command == 5:
		interface.displayCue(val1,val2)

	if command == 6:	# Set Target
		fing = ((val1-1) *val2 + val3/100)/80
		interface.setFing(fing)
		# print(fing)
		interface.render()

	if command == 9:	
		
		robot_orn[0] = (val1-1) *(val2 + val3/100)/10
		robot_orn[1] = -(val4-1) *(val5 + val6/100)/10
		robot_orn[2] = (val7-1) *(val8 + val9/100)/10

		robot_pos

		print(robot_orn)
		interface.robotenv.set_robotOrn(robot_orn)
		interface.render()

	if command == 14:	
		key = val10
		robot_pos[0] = ((val1-1) *(val2 + val3/100))/ 1000
		robot_pos[1] = -((val4-1) *(val5 + val6/100))/ 1000
		robot_pos[2] = ((val7-1) *(val8 + val9/100))/ 1000
		interface.updateRobotPos(robot_pos,key )
		print(robot_pos)
		interface.render()

	if command == 15:		#set path coordinates
		ind = val10
		p[0] = ((val1-1) *(val2 + val3/100))/ 1000
		p[1] = -((val4-1) *(val5 + val6/100))/ 1000
		p[2] = ((val7-1) *(val8 + val9/100))/ 1000
		interface.setPath_ES(p, ind)


	if command == 16:		#set robot orn
		ind = val10
		robot_theta_delta = ((val1-1) *(val2 + val3/100))/100
		interface.robotenv.set_robotRotation(1, robot_theta_delta)
		print(robot_theta_delta)
		interface.render()

sock.shutdown()
sock.close() 
