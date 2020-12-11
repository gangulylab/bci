import interfaces
import numpy as np
import time
from sixDirLabels_UI import UI
import pygame as pg

interface = interfaces.DiscreteActionsRobot()
interface.mode = 3
interface.angle = 0
interface.debugLines = 1
interface.open()
robot_open = 1
target_pos = np.array([0.0, 300.0])

logCursor = 0
updateRate = 0.125
pg.init()
pg.display.set_mode((600,500))
pg.display.set_caption("Control Interface")
runUI = UI(logCursor, updateRate)

target_pos = np.array([[0., -0.2, 0.0], [0.0, 0.2, 0.0],  [0.2, 0.0, 0.0], [-0.2, 0.0, 0.0], [ 0.0, 0.0, 0.2], [ 0.0, 0.0, -0.2]])
interface.updateRefresh(0.125);

runUI.mode = interface.mode

while True:
	i = np.random.randint(6)
	current_target = target_pos[i]
	print(current_target)
	interface.reset()
	interface.create_target3D(current_target,1)


	pos = np.array([0., 0., 0.])
	pos[0] = interface.robotenv.pos[0] - interface.robotenv.center[0]
	pos[1] = interface.robotenv.pos[1] - interface.robotenv.center[1]
	pos[2] = interface.robotenv.pos[2] 
	dist = np.linalg.norm(pos - current_target)

	if interface.mode == 1:
		notdone = 1
		while notdone == 1:
			time.sleep(0.125)

			runUI.update()
			key = runUI.state
			interface.update_joystick(key)
			interface.render()

			if interface.target == 1:
				
				if interface.robotenv.fing < .01:
					notdone = 0
			if interface.target == 3:
				if interface.robotenv.fing > 1.34:
					notdone = 0
			if interface.target == 2:
				if interface.robotenv.pos[2] > 0.4:
					notdone = 0
			if interface.target == 4:
				if interface.robotenv.pos[2] < 0.05:
					notdone = 0
	else:
		while dist > 0.04:
			pos[0] = interface.robotenv.pos[0] - interface.robotenv.center[0]
			pos[1] = interface.robotenv.pos[1] - interface.robotenv.center[1]
			pos[2] = interface.robotenv.pos[2] - interface.robotenv.center[2]

			dist = np.linalg.norm(pos - current_target)
			# print(dist)
			time.sleep(0.125)

			runUI.update()
			key = runUI.state
			# print(key)
			interface.update_joystick(key)
			interface.render()
	time.sleep(1.0)


	# if command == 0:
	# 	if val1 == 0:		# Open Robot 
	# 		if robot_open == 0:
	# 			interface.open()
	# 			robot_open = 1
	# 	if val1 == 1:		# Reset 
	# 		interface.reset()
	# 	if val1 == 2:		# Change hold time on debug lines
	# 		updateRate = 1.0 /val2;
	# 		interface.updateRefresh(updateRate);
	# if command == 1:	# Set Target
	# 	# interface.reset()
	# 	target_pos[0] = (val1 - 128) / 100
	# 	target_pos[1] = (val2 - 128) / 100
	# 	interface.create_target(target_pos)
	# 	print(val1, val2)
	# 	print(target_pos)
	# if command == 2:	# Set Dirr
	# 	key = val1
