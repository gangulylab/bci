import pybullet as p
import pybullet_data

class HandEnv(object):
	def __init__(self):
		p.connect(p.GUI)
		p.setAdditionalSearchPath(pybullet_data.getDataPath())
		p.resetSimulation()
		p.configureDebugVisualizer(p.COV_ENABLE_GUI,0)
		p.resetDebugVisualizerCamera(cameraDistance=0.22, cameraYaw=0, cameraPitch=-80, cameraTargetPosition=[-0.5,0.13,0.2])


		p.loadURDF("table/table.urdf", basePosition=[-0.5,-0.1,-0.65])
		self.handId = p.loadURDF("URDFs/asr_hand/model_hand_right2.urdf", [-0.5,0,-1.0],  useFixedBase=True)
		
		p.configureDebugVisualizer(p.COV_ENABLE_SHADOWS, 0)

		self.handStates = [1,2,3]
		self.jointState = [0,0,-0.4,0]
		self.hideHand()
		p.resetJointState(self.handId, 3, -0.4)
		self.showDecodes = 0

	def setJointPosition(self, joint, angle):
		self.jointState[joint] = angle

	def setGrasp(self, action, s):

		print(action)
		print(s)


		l = 0.6

		if action == 1:		# thumb
			p.resetJointState(self.handId, 25, l*s+0.1)
			p.resetJointState(self.handId, 26, l*s+ 0.5)
			p.resetJointState(self.handId, 27, l*s+ 0.2)

		if action == 2:		# index
			p.resetJointState(self.handId, 21, s)
			p.resetJointState(self.handId, 22, s)
			p.resetJointState(self.handId, 23, s)

		if action == 3:	# middle
			p.resetJointState(self.handId, 16, s)
			p.resetJointState(self.handId, 17, s)
			p.resetJointState(self.handId, 18, s)

		elif action == 4:  # ring 
			p.resetJointState(self.handId, 11, s)
			p.resetJointState(self.handId, 12, s)
			p.resetJointState(self.handId, 13, s)

		if action == 5:  # pinky
			p.resetJointState(self.handId, 6, s)
			p.resetJointState(self.handId, 7, s)
			p.resetJointState(self.handId, 8, s)

		elif action == 6:  # Power Grasp	
			p.resetJointState(self.handId, 3, -.4)

			# p.resetJointState(self.handId, 5, -0.5-0.2*s )
			p.resetJointState(self.handId, 6, s)
			p.resetJointState(self.handId, 7, s)
			p.resetJointState(self.handId, 8, s)

			# p.resetJointState(self.handId, 10, -0.2-0.1*s )
			p.resetJointState(self.handId, 11, s)
			p.resetJointState(self.handId, 12, s)
			p.resetJointState(self.handId, 13, s)

			# p.resetJointState(self.handId, 15, 0.1+0.1*s )
			p.resetJointState(self.handId, 16, s)
			p.resetJointState(self.handId, 17, s)
			p.resetJointState(self.handId, 18, s)

			# p.resetJointState(self.handId, 20, 0.4+0.2*s )
			p.resetJointState(self.handId, 21, s)
			p.resetJointState(self.handId, 22, s)
			p.resetJointState(self.handId, 23, s)

			p.resetJointState(self.handId, 25, s+0.1)
			p.resetJointState(self.handId, 26, s+ 0.5)
			p.resetJointState(self.handId, 27, s+ 0.2)

		elif action == 7:  #PINCH GRASP

			p.resetJointState(self.handId, 3, 0.2)
			p.resetJointState(self.handId, 24, 0.5)
			p.resetJointState(self.handId, 25, s*0.1)
			p.resetJointState(self.handId, 26, 0)
			p.resetJointState(self.handId, 27, s+ 0.2)


			p.resetJointState(self.handId, 20, 0.4+0.2*s )
			p.resetJointState(self.handId, 21, s)
			p.resetJointState(self.handId, 22, s)
			p.resetJointState(self.handId, 23, s)

		elif action == 8:  #TRIPOD GRASP

			p.resetJointState(self.handId, 3, 0.)
			p.resetJointState(self.handId, 24, 0.5)
			p.resetJointState(self.handId, 25, s*0.1)
			p.resetJointState(self.handId, 26, 0)
			p.resetJointState(self.handId, 27, s+ 0.2)

			p.resetJointState(self.handId, 20, 0.4+0.2*s )
			p.resetJointState(self.handId, 21, s)
			p.resetJointState(self.handId, 22, s)
			p.resetJointState(self.handId, 23, s)

			p.resetJointState(self.handId, 15, 0.1+0.1*s )
			p.resetJointState(self.handId, 16, s)
			p.resetJointState(self.handId, 17, s)
			p.resetJointState(self.handId, 18, s)
		
		elif action == 9:  #WRIST ADDUCTION (OUT)
			p.resetJointState(self.handId, 1, -s*0.5)
		
		elif action == 10:  #WRIST ABDUCTION (IN)
			p.resetJointState(self.handId, 1, s*0.5)

		elif action == 11:  #WRIST FLEX
			p.resetJointState(self.handId, 2, -s*0.5)
		
		elif action == 12:  #WRIST EXTEND
			p.resetJointState(self.handId, 2, s*0.5)

		elif action == 13:  # WRIST SUPINATE
			p.resetJointState(self.handId, 25, 0)
			p.resetJointState(self.handId, 26,  0)
			p.resetJointState(self.handId, 27,  0)
			p.resetJointState(self.handId, 3, 2*s+1)


		elif action == 14:  #WRIST PRONATE
			p.resetJointState(self.handId, 25, 0)
			p.resetJointState(self.handId, 26, 0)
			p.resetJointState(self.handId, 27, 0)
			p.resetJointState(self.handId, 3, -2*s + 1)

	def reset(self):
		# p.resetJointState(self.handId, 3, -0.4)
		print('RESET')
		self.hideHand()
		p.resetJointState(self.handId, 3, -0.4)

		p.resetJointState(self.handId, 1, 0)
		p.resetJointState(self.handId, 2, 0)
		
		p.resetJointState(self.handId, 4, 0)
		p.resetJointState(self.handId, 5, 0)
		p.resetJointState(self.handId, 6, 0)
		p.resetJointState(self.handId, 7, 0)
		p.resetJointState(self.handId, 8, 0)

		p.resetJointState(self.handId, 10, 0) 
		p.resetJointState(self.handId, 11, 0)
		p.resetJointState(self.handId, 12, 0)
		p.resetJointState(self.handId, 13, 0)

		p.resetJointState(self.handId, 15, 0)
		p.resetJointState(self.handId, 16, 0)
		p.resetJointState(self.handId, 17, 0)
		p.resetJointState(self.handId, 18, 0)

		p.resetJointState(self.handId, 20, 0)
		p.resetJointState(self.handId, 21, 0)
		p.resetJointState(self.handId, 22, 0)
		p.resetJointState(self.handId, 23, 0)
		p.resetJointState(self.handId, 24, 0.0)
		p.resetJointState(self.handId, 25, 0.1)
		p.resetJointState(self.handId, 26, 0.5)
		p.resetJointState(self.handId, 27, 0.2)
		p.removeAllUserDebugItems()
		
		if self.showDecodes:
			self.displayDecodes()


	def displayDecodes(self):
		# self.reset()
		c = [1, 1, 1]
		t = 0
		d = 10

		p1 = [-0.58,0.11,0.21]
		p2 = [-0.58,0.15, 0.21]
		p.addUserDebugLine(p1, p2, c, d, t)

		p1 = [-0.515,0.19,0.21]
		p2 = [-0.515,0.23, 0.21]
		p.addUserDebugLine(p1, p2, c, d, t)
	
		p1 = [-0.49,0.20,0.21]
		p2 = [-0.49,0.24, 0.21]
		p.addUserDebugLine(p1, p2, c, d, t)
	
		p1 = [-0.47,0.19,0.21]
		p2 = [-0.47,0.23, 0.21]
		p.addUserDebugLine(p1, p2, c, d, t)
	
		p1 = [-0.45,0.16,0.21]
		p2 = [-0.45,0.20, 0.21]
		p.addUserDebugLine(p1, p2, c, d, t)

		p1 = [-0.6,0.0,0.21]
		p2 = [-0.7,0.04, 0.21]
		p.addUserDebugLine(p1, p2, c, d, t)
	
		p1 = [-0.4,0.0,0.21]
		p2 = [-0.3,0.04, 0.21]
		p.addUserDebugLine(p1, p2, c, d, t)


		p1 = [-0.43,0.06,0.21]
		p2 = [-0.43,0.01, 0.21]
		p.addUserDebugLine(p1, p2, c, d, t)
	
		p1 = [-0.43,0.0,0.21]
		p2 = [-0.43,-0.05, 0.21]
		p.addUserDebugLine(p1, p2, c, d, t)


		p1 = [-0.6,-0.02,0.21]
		p2 = [-0.7,-0.02, 0.21]
		p.addUserDebugLine(p1, p2, c, d, t)
	
		p1 = [-0.4,-0.02,0.21]
		p2 = [-0.3,-0.02, 0.21]
		p.addUserDebugLine(p1, p2, c, d, t)


	def displayCurrentDecode(self, action):
		c = [1, 0, 0]
		t = 0.2
		d = 15

		if action == 1:
			p1 = [-0.58,0.11,0.2101]
			p2 = [-0.58,0.15, 0.2101]
			p.addUserDebugLine(p1, p2, c, d, t)
		if action == 2:
			p1 = [-0.515,0.19,0.2101]
			p2 = [-0.515,0.23, 0.2101]
			p.addUserDebugLine(p1, p2, c, d, t)
		if action == 3:
			p1 = [-0.49,0.20,0.2101]
			p2 = [-0.49,0.24, 0.2101]
			p.addUserDebugLine(p1, p2, c, d, t)
		if action == 4:
			p1 = [-0.47,0.19,0.2101]
			p2 = [-0.47,0.23, 0.2101]
			p.addUserDebugLine(p1, p2, c, d, t)
		if action == 5:
			p1 = [-0.45,0.16,0.2101]
			p2 = [-0.45,0.20, 0.2101]
			p.addUserDebugLine(p1, p2, c, d, t)
		if action == 6:
			p1 = [-0.58,0.11,0.2101]
			p2 = [-0.58,0.15, 0.2101]
			p.addUserDebugLine(p1, p2, c, d, t)
			p1 = [-0.515,0.19,0.2101]
			p2 = [-0.515,0.23, 0.2101]
			p.addUserDebugLine(p1, p2, c, d, t)
			p1 = [-0.49,0.20,0.2101]
			p2 = [-0.49,0.24, 0.2101]
			p.addUserDebugLine(p1, p2, c, d, t)
			p1 = [-0.47,0.19,0.2101]
			p2 = [-0.47,0.23, 0.2101]
			p.addUserDebugLine(p1, p2, c, d, t)
			p1 = [-0.45,0.16,0.2101]
			p2 = [-0.45,0.20, 0.2101]
			p.addUserDebugLine(p1, p2, c, d, t)
		if action == 7:  # pinch
			p1 = [-0.58,0.11,0.2101]
			p2 = [-0.58,0.15, 0.2101]
			p.addUserDebugLine(p1, p2, c, d, t)
			p1 = [-0.515,0.19,0.2101]
			p2 = [-0.515,0.23, 0.2101]
			p.addUserDebugLine(p1, p2, c, d, t)	
		if action == 8:
			p1 = [-0.58,0.11,0.2101]
			p2 = [-0.58,0.15, 0.2101]
			p.addUserDebugLine(p1, p2, c, d, t)
			p1 = [-0.515,0.19,0.2101]
			p2 = [-0.515,0.23, 0.2101]
			p.addUserDebugLine(p1, p2, c, d, t)
			p1 = [-0.49,0.20,0.2101]
			p2 = [-0.49,0.24, 0.2101]
			p.addUserDebugLine(p1, p2, c, d, t)
		if action == 10:	
			p1 = [-0.6,0.0,0.2101]
			p2 = [-0.7,0.04, 0.2101]
			p.addUserDebugLine(p1, p2, c, d, t)
		if action == 9:
			p1 = [-0.4,0.0,0.2101]
			p2 = [-0.3,0.04, 0.2101]
			p.addUserDebugLine(p1, p2, c, d, t)

		if action == 11:	
			p1 = [-0.43,0.06,0.2101]
			p2 = [-0.43,0.01, 0.2101]
			p.addUserDebugLine(p1, p2, c, d, t)
		if action == 12:	
			p1 = [-0.43,0.0,0.2101]
			p2 = [-0.43,-0.05, 0.2101]
			p.addUserDebugLine(p1, p2, c, d, t)

		if action == 13:
			p1 = [-0.6,-0.02,0.2101]
			p2 = [-0.7,-0.02, 0.2101]
			p.addUserDebugLine(p1, p2, c, d, t)
		if action == 14:
			p1 = [-0.4,-0.02,0.2101]
			p2 = [-0.3,-0.02, 0.2101]
			p.addUserDebugLine(p1, p2, c, d, t)
			


	def hideHand(self):
		p.resetBasePositionAndOrientation(self.handId,[0,0,-1.0],[0,0,0,1])

	def showHand(self):
		a = p.getBasePositionAndOrientation(self.handId)
		p.resetBasePositionAndOrientation(self.handId,[-0.5,0,0.2],a[1])

	def displayCue(self, cue,c):
		p.removeAllUserDebugItems()
		if cue == 1:
			p.addUserDebugText('THUMB', [-0.6,0.25,0.2],  [0,0,c], 12)
		elif cue == 2:
			p.addUserDebugText('INDEX', [-0.61,0.25,0.2],  [0,0,c], 12)
		elif cue == 3:
			p.addUserDebugText('MIDDLE', [-0.65,0.25,0.2],  [0,0,c], 12)
		elif cue == 4:
			p.addUserDebugText('RING', [-0.6,0.25,0.2],  [0,0,c], 12)
		elif cue == 5:
			p.addUserDebugText('PINKY', [-0.65,0.25,0.2],  [0,0,c], 12)
		elif cue == 6:
			p.addUserDebugText('POWER GRASP', [-0.75,0.25,0.2],  [0,0,c], 12)
		elif cue == 7:
			p.addUserDebugText('PINCH GRASP', [-0.75,0.25,0.2],  [0,0,c], 12)
		elif cue == 8:
			p.addUserDebugText('TRIPOD GRASP', [-0.75,0.25,0.2],  [0,0,c], 12)
		elif cue == 9:
			p.addUserDebugText('WRIST OUT', [-0.7,0.25,0.2],  [0,0,c], 12)
		elif cue == 10:
			p.addUserDebugText('WRIST IN', [-0.7,0.25,0.2],  [0,0,c], 12)
		elif cue == 11:
			p.addUserDebugText('WRIST FLEX', [-0.7,0.25,0.2],  [0,0,c], 12)
		elif cue == 12:
			p.addUserDebugText('WRIST EXTEND', [-0.75,0.25,0.2],  [0,0,c], 12)
		elif cue == 13:
			p.addUserDebugText('WRIST PRONATE', [-0.75,0.25,0.2],  [0,0,c], 12)
		elif cue == 14:
			p.addUserDebugText('WRIST SUPINATE', [-0.75,0.25,0.2],  [0,0,c], 12)

		self.setGrasp(cue, 0.1)
		if self.showDecodes:
			self.displayDecodes()

	def step(self):
		
		for i in self.handStates:
			p.resetJointState(self.handId, i, self.jointState[i])
			# print(self.jointState[i])

