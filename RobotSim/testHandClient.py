from HandEnv import HandEnv

env = HandEnv()


s = 0
while True:

	s = s+ .00001
	joint = 3
	env.setJointPosition(2,s)
	env.setJointPosition(3,s)
	env.setJointPosition(1,s)

	env.step()




