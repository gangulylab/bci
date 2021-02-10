import pybullet as p
import pybullet_data

class HandEnv(object):
	def __init__(self):
		p.connect(p.GUI)
		p.setAdditionalSearchPath(pybullet_data.getDataPath())
		p.resetSimulation()
		p.configureDebugVisualizer(p.COV_ENABLE_GUI,0)
		p.resetDebugVisualizerCamera(cameraDistance=0.3, cameraYaw=0, cameraPitch=-80, cameraTargetPosition=[-0.5,0.1,0.2])


		p.loadURDF("table/table.urdf", basePosition=[-0.6,0.0,-0.65])
		self.handId = p.loadURDF("URDFs/asr_hand/model_hand_right2.urdf", [-0.5,0,0.2],  useFixedBase=True)
		p.resetJointState(self.handId, 3, -0.4)

		self.handStates = [1,2,3]
		self.jointState = [0,0,0,0]

	def setJointPosition(self, joint, angle):
		self.jointState[joint] = angle

	def reset(self):
		p.resetBasePositionAndOrientation(self.handId, [-1., -1., -1.], [0,0,0,1])
		p.resetJointState(handId, 3, -0.4)

	def step(self):
		
		for i in self.handStates:
			p.resetJointState(self.handId, i, self.jointState[i])
			# print(self.jointState[i])

