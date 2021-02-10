import pybullet as p
import time
import math
from datetime import datetime
import pybullet_data
# import pygame as pg
import random
import numpy as np
import os
import csv
from scipy.spatial.transform import Rotation as R

use2D   = 0
logData = 1

clid = p.connect(p.SHARED_MEMORY)
if (clid<0):
  p.connect(p.GUI)

p.setAdditionalSearchPath(pybullet_data.getDataPath())
p.resetSimulation()
p.configureDebugVisualizer(p.COV_ENABLE_GUI,0)

# p.loadURDF("plane.urdf",[0,0,-.65])
p.loadURDF("table/table.urdf", basePosition=[-0.6,0.0,-0.65])
handId = p.loadURDF("URDFs/asr_hand/model_hand_right2.urdf", [-0.5,0,0.2],  useFixedBase=True)

basePos = [0,0,0]
# p.resetBasePositionAndOrientation(handId,basePos,[0,0,0,1])
p.resetDebugVisualizerCamera(cameraDistance=0.3, cameraYaw=0, cameraPitch=-80, cameraTargetPosition=[-0.5,0.1,0.2])


# p.resetJointState(handId, 2, -0.0)
p.resetJointState(handId, 3, -0.4)


l1 = [-0.5, 0.5]
l2 = [-0.5, 1.2]
l3 = [-0.5, 3.14]
l2 = [-1.2, 0.5]

l4 = [-1.0, 0.1]
s = 0
dir  = 1
while True:

#   if dir == 1: 
#     if s <= l1[1]:
#       s = s+.00001
#     else:
#       dir = 0
#   else:
#     if s >= l1[0]:
#       s = s-.00001
#     else:
#       dir = 1
# # 
#   # print(s)
#   p.resetJointState(handId, 1, s)


  # if dir == 1: 
  #   if s <= l2[1]:
  #     s = s+.00001
  #   else:
  #     dir = 0
  # else:
  #   if s >= l2[0]:
  #     s = s-.00001
  #   else:
  #     dir = 1

  # # print(s)
  # p.resetJointState(handId, 1, s)

  if dir == 1: 
    if s <= l4[1]:
      s = s+.0001
    else:
      dir = 0
      time.sleep(1.0)
  else:
    if s >= l4[0]:
      s = s-.0001
    else:
      dir = 1
      time.sleep(1.0)


  p.resetJointState(handId, 5, -0.5-0.2*s )
  p.resetJointState(handId, 6, s)
  p.resetJointState(handId, 7, s)
  p.resetJointState(handId, 8, s)


  p.resetJointState(handId, 10, -0.2-0.1*s )
  p.resetJointState(handId, 11, s)
  p.resetJointState(handId, 12, s)
  p.resetJointState(handId, 13, s)

  p.resetJointState(handId, 15, 0.1+0.1*s )
  p.resetJointState(handId, 16, s)
  p.resetJointState(handId, 17, s)
  p.resetJointState(handId, 18, s)

  p.resetJointState(handId, 20, 0.4+0.2*s )
  p.resetJointState(handId, 21, s)
  p.resetJointState(handId, 22, s)
  p.resetJointState(handId, 23, s)


  p.resetJointState(handId, 25, s+0.1)
  p.resetJointState(handId, 26, s+ 0.5)
  p.resetJointState(handId, 27, s+ 0.2)



  # if dir == 1: 
  #   if s <= l3[1]:
  #     s = s+.00001
  #   else:
  #     dir = 0
  # else:
  #   if s >= l3[0]:
  #     s = s-.00001
  #   else:
  #     dir = 1

  # # print(s)
  # p.resetJointState(handId, 3, s)

  # p.resetJointState(handId, 3, s)
  # p.resetJointState(handId, 2, s)



  # p.resetJointState(handId, , 1.0)
# numJoints = 10
# jacoArmJoints = [2, 3, 4, 5, 6, 7]
# jacoFingerJoints = [9, 11, 13]

# #joint damping coefficents
# jd = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
# rp = [0,math.pi/4,math.pi,1.0*math.pi, 1.8*math.pi, 0*math.pi, 1.75*math.pi, 0.5*math.pi]

# wu = [0.1, 0.5, 0.5]
# wl = [-.66, -.5, 0.00]

# for i in range(8):
#   p.resetJointState(jacoId,i, rp[i])

# ls = p.getLinkState(jacoId, jacoEndEffectorIndex)
# p.setGravity(0,0,-10)

# t = 0.
# prevPose = [0, 0, 0]
# prevPose1 = [0, 0, 0]
# hasPrevPose = 0
# useNullSpace = 0

# useOrientation = 1
# useSimulation = 1
# useRealTimeSimulation = 1
# ikSolver = 0
# p.setRealTimeSimulation(useRealTimeSimulation)
# p.configureDebugVisualizer(p.COV_ENABLE_RENDERING,1) 
# trailDuration = 5

# pos = list(ls[4])
# orn = list(ls[5])

# # print(ls)
# i=0

# JP = list(rp[2:9])
# fing = 0
# wri = 0

# newPosInput = 1
# keyT = time.time()
# kTotal = np.zeros([9,], dtype = int)

# pg.key.set_repeat()
# kp5_up = 1
# add_KP5 = 1

# dist = .001
# ang = .005
# rot_theta = .005
# inputRate = .05

# Rx = np.array([[1., 0., 0.],[0., np.cos(rot_theta), -np.sin(rot_theta)], [0., np.sin(rot_theta), np.cos(rot_theta)]])
# Ry = np.array([[np.cos(rot_theta), 0., np.sin(rot_theta)], [0., 1., 0.], [-np.sin(rot_theta), 0., np.cos(rot_theta)]])
# Rz = np.array([[np.cos(rot_theta), -np.sin(rot_theta), 0.], [np.sin(rot_theta), np.cos(rot_theta), 0.], [0., 0., 1.]])

# Rxm = np.array([[1., 0., 0.],[0., np.cos(-rot_theta), -np.sin(-rot_theta)], [0., np.sin(-rot_theta), np.cos(-rot_theta)]])
# Rym = np.array([[np.cos(-rot_theta), 0., np.sin(-rot_theta)], [0., 1., 0.], [-np.sin(-rot_theta), 0., np.cos(-rot_theta)]])
# Rzm = np.array([[np.cos(-rot_theta), -np.sin(-rot_theta), 0.], [np.sin(-rot_theta), np.cos(-rot_theta), 0.], [0., 0., 1.]])

# updateT = time.time()

# if logData:
#   dataDir = time.strftime("%Y%m%d")
#   if not os.path.exists(dataDir):
#     os.makedirs(dataDir)
#   trialInd = len(os.listdir(dataDir))
#   fn = dataDir + "/simdata00" + str(trialInd) + ".csv"
      
#   logFile = open(fn, 'w', newline = '')
#   fileObj = csv.writer(logFile)

# while 1:

#   i+=1
#   if (useRealTimeSimulation):
#     dt = datetime.now()
#     t = (dt.second / 60.) * 2. * math.pi
#   else:
#     t = t + 0.01

#   if (useSimulation and useRealTimeSimulation == 0):
#     p.stepSimulation()

#   delta = time.time() - updateT 

#   if delta > inputRate:
#     updateT= time.time()

#     if logData:
#       lsr = p.getLinkState(jacoId, jacoEndEffectorIndex)
#       lsc = p.getBasePositionAndOrientation(cube1Id)

#       ln = [lsr[4][0],lsr[4][1],lsr[4][2], lsr[5][0],lsr[5][1], lsr[5][2], lsr[5][3], fing, lsc[0][0],lsc[0][1],lsc[0][2], lsc[1][0],lsc[1][1], lsc[1][2], lsc[1][3]] 
#       ln_rnd = [round(num, 4) for num in ln]
#       fileObj.writerow(ln_rnd)

#     eulOrn = p.getEulerFromQuaternion(orn)
#     Rrm = R.from_quat(orn)

#     rx = eulOrn[0]
#     ry = eulOrn[1]
#     rz = eulOrn[2]

#     runUI.update()
#     inputMode = runUI.mode
#     inputKey  = runUI.state

#     baseTheta = JP[0]
#     s = math.sin(baseTheta)
#     c = math.cos(baseTheta)

#     c1 = math.cos(ang)
#     s1 = math.sin(ang)

#     c2 = math.cos(-ang)
#     s2 = math.sin(-ang)

#     n = np.sqrt(pos[0]*pos[0] + pos[1]*pos[1])
#     dx = -pos[1]/n
#     dy = pos[0]/n

#     Rnew =  Rrm.as_matrix() 

    
#   if (newPosInput == 1):
#     if (useNullSpace == 1):
#       if (useOrientation == 1):
#         jointPoses = p.calculateInverseKinematics(jacoId, jacoEndEffectorIndex, pos, orn, ll, ul,
#                                                   jr, rp)
#       else:
#         jointPoses = p.calculateInverseKinematics(jacoId,
#                                                   jacoEndEffectorIndex,
#                                                   pos,
#                                                   lowerLimits=ll,
#                                                   upperLimits=ul,
#                                                   jointRanges=jr,
#                                                   restPoses=rp)
#     else:
#       if (useOrientation == 1):
#         jointPoses = p.calculateInverseKinematics(jacoId,
#                                                   jacoEndEffectorIndex,
#                                                   pos,
#                                                   orn,
#                                                   jointDamping=jd,
#                                                   solver=ikSolver,
#                                                   maxNumIterations=100,
#                                                   residualThreshold=.01)
#         JP = list(jointPoses)
#       else:
#         jointPoses = p.calculateInverseKinematics(jacoId,
#                                                   jacoEndEffectorIndex,
#                                                   pos,
#                                                   solver=ikSolver)
#         JP = list(jointPoses)

#   if (useSimulation):
#     JS = p.getJointStates(jacoId, [1, 2, 3, 4, 5, 6, 7, 9, 11, 13])
#     j = 0
#     for i in [2,3,4,5,6,7]:
#       p.setJointMotorControl2(jacoId, i, p.POSITION_CONTROL, JP[j])
#       j = j+1
    
#     for i in  [9, 11, 13]:
#       p.setJointMotorControl2(jacoId, i, p.POSITION_CONTROL, fing)

#   else:
#     j = 0
#     for i in jacoArmJoints:
#       p.resetJointState(jacoId, i, jointPoses[j])
#       j = j+1

#   ls = p.getLinkState(jacoId, jacoEndEffectorIndex)

#   prevPose = tuple(pos)
#   prevPose1 = ls[4]
#   hasPrevPose = 1
#   newPosInput = 0


# file.close()
# p.disconnect()