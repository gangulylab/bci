import pybullet as p
import time
import math
from datetime import datetime
import pybullet_data

import random
import numpy as np
import os
import csv
from scipy.spatial.transform import Rotation as R

class JacoEnv(object):
  def __init__(self):
    p.connect(p.GUI)
    p.setAdditionalSearchPath(pybullet_data.getDataPath())
    # p.setAdditionalSearchPath('../URDFs')
    p.resetSimulation()
    p.configureDebugVisualizer(p.COV_ENABLE_GUI,0)
    p.setGravity(0,0,-10)

    rot_theta   = .005
    self.Rx = np.array([[1., 0., 0.],[0., np.cos(rot_theta), -np.sin(rot_theta)], [0., np.sin(rot_theta), np.cos(rot_theta)]])
    self.Ry = np.array([[np.cos(rot_theta), 0., np.sin(rot_theta)], [0., 1., 0.], [-np.sin(rot_theta), 0., np.cos(rot_theta)]])
    self.Rz = np.array([[np.cos(rot_theta), -np.sin(rot_theta), 0.], [np.sin(rot_theta), np.cos(rot_theta), 0.], [0., 0., 1.]])

    self.Rxm = np.array([[1., 0., 0.],[0., np.cos(-rot_theta), -np.sin(-rot_theta)], [0., np.sin(-rot_theta), np.cos(-rot_theta)]])
    self.Rym = np.array([[np.cos(-rot_theta), 0., np.sin(-rot_theta)], [0., 1., 0.], [-np.sin(-rot_theta), 0., np.cos(-rot_theta)]])
    self.Rzm = np.array([[np.cos(-rot_theta), -np.sin(-rot_theta), 0.], [np.sin(-rot_theta), np.cos(-rot_theta), 0.], [0., 0., 1.]])

    self.useOrientation = 1
    self.useSimulation = 0
    self.useRealTimeSimulation = 1

    self.newPosInput = 0
    self.ikSolver = 0

    self.inputRate = .05
    self.reset()

    self.cubeStep = 3
    self.hasCube = False
    self.dist = 0.01
    self.debuglen = 0.3

    self.wu = [0, 0.2]
    self.wl = [-.7, -.8]
    self.center = np.array([-0.35, -0.3])
    self.bciRate = 0.125

  def set_block_pos(self, pos):
    pos = self.center + pos
    self.cube1Id = p.loadURDF("cube_small.urdf",[pos[0], pos[1], 0], [0,0,0,1])
    d = .08
    c1 = [pos[0] - d, pos[1] - d, 0]
    c2 = [pos[0] - d, pos[1] + d, 0]
    c3 = [pos[0] + d, pos[1] + d, 0]
    c4 = [pos[0] + d, pos[1] - d, 0]
  

    p.addUserDebugLine(c1, c2, [0,0,1], 3, 0)
    p.addUserDebugLine(c2, c3, [0,0,1], 3, 0)
    p.addUserDebugLine(c3, c4, [0,0,1], 3, 0)
    p.addUserDebugLine(c4, c1, [0,0,1], 3, 0)
    # p.resetBasePositionAndOrientation(self.cube1Id, [pos[0], pos[1], 0], [0,0,0,1])

  def updateCommand(self, key):
  
    Rrm = R.from_quat(self.orn)
    Rnew =  Rrm.as_dcm() 

    baseTheta = self.JP[0]
    s = math.sin(baseTheta)
    c = math.cos(baseTheta)

    n = np.sqrt(self.pos[0]*self.pos[0] + self.pos[1]*self.pos[1])
    dx = -self.pos[1]/n
    dy = self.pos[0]/n


    if 0 < key < 13:
      self.newPosInput = 1

    # Position
    if key == 6:
      self.pos[0] = self.pos[0] + self.dist
      self.pos2 = [self.pos[0] + self.debuglen, self.pos[1], self.pos[2]]
    elif key == 4: 
      self.pos[0] = self.pos[0] - self.dist
      self.pos2 = [self.pos[0] - self.debuglen, self.pos[1], self.pos[2]]
    elif key == 8:
      self.pos[1] = self.pos[1] + self.dist
      self.pos2 = [self.pos[0], self.pos[1] + self.debuglen, self.pos[2]]
    elif key == 2:
      self.pos[1] = self.pos[1] - self.dist
      self.pos2 = [self.pos[0], self.pos[1] - self.debuglen, self.pos[2]]

    if self.pos[0] > self.wu[0]:
      self.pos[0] =  self.wu[0]
    if self.pos[0] < self.wl[0]:
      self.pos[0] =  self.wl[0]
    if self.pos[1] > self.wu[1]:
      self.pos[1] =  self.wu[1]
    if self.pos[1] < self.wl[1]:
      self.pos[1] =  self.wl[1]

    if key == 6:
      self.pos2 = [self.pos[0] + self.debuglen, self.pos[1], self.pos[2]]
    elif key == 4: 
      self.pos2 = [self.pos[0] - self.debuglen, self.pos[1], self.pos[2]]
    elif key == 8:
      self.pos2 = [self.pos[0], self.pos[1] + self.debuglen, self.pos[2]]
    elif key == 2:
      self.pos2 = [self.pos[0], self.pos[1] - self.debuglen, self.pos[2]]


  def inverseKin(self):
    if (self.newPosInput == 1):
      self.jointPoses = p.calculateInverseKinematics(self.jacoId,
                                                self.jacoEndEffectorIndex,
                                                self.pos,
                                                self.orn,
                                                jointDamping=self.jd,
                                                solver=self.ikSolver,
                                                maxNumIterations=100,
                                                residualThreshold=.01)

      self.JP = list(self.jointPoses)

  def reset(self):
    p.resetSimulation()

    p.loadURDF("plane.urdf",[0,0,-.65])
    p.loadURDF("table/table.urdf", basePosition=[-0.6,-0.3,-0.65])
    self.jacoId = p.loadURDF("URDFs/jaco/j2n6s300.urdf", [0,0,0],  useFixedBase=True)

    p.resetBasePositionAndOrientation(self.jacoId,[0,0,0],[0,0,0,1])
    # p.resetDebugVisualizerCamera(cameraDistance=0.8, cameraYaw=0, cameraPitch=-89.99, cameraTargetPosition=[-0.35,-0.3,0.0])
    p.resetDebugVisualizerCamera(cameraDistance=0.8, cameraYaw=0, cameraPitch=-30, cameraTargetPosition=[-0.35,-0.3,0.0])
    
    # p.resetDebugVisualizerCamera(cameraDistance=0.8, cameraYaw=30, cameraPitch=-30, cameraTargetPosition=[-0.35,-0.3,0.0])
    
    self.jacoEndEffectorIndex = 8
    self.jacoArmJoints = [2, 3, 4, 5, 6, 7]
    jacofingerJoints = [9, 11, 13]

    self.jd = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    rp = [0,math.pi/4,math.pi,1.0*math.pi, 1.8*math.pi, 0*math.pi, 1.75*math.pi, 0.5*math.pi]

    self.pos =list([-0.35, -0.3, 0.1])
    self.orn = p.getQuaternionFromEuler([0,math.pi,0])

    for i in range(8):
      p.resetJointState(self.jacoId,i, rp[i])
    
    self.newPosInput = 1
    self.inverseKin()
    for i in self.jacoArmJoints:
      p.resetJointState(self.jacoId,i, self.JP[i-2])

    ls = p.getLinkState(self.jacoId, self.jacoEndEffectorIndex)
    p.setRealTimeSimulation(self.useRealTimeSimulation)
    p.configureDebugVisualizer(p.COV_ENABLE_RENDERING,1) 

    self.fing = 0
    self.updateT = time.time()
    self.pos2 = self.pos

  def step(self):

    if self.newPosInput: 
      self.inverseKin()

    if (self.useSimulation):     
      self.JS = p.getJointStates(self.jacoId, [1, 2, 3, 4, 5, 6, 7, 9, 11, 13])
      j = 0
      for i in self.jacoArmJoints:
        p.setJointMotorControl2(self.jacoId, i, p.POSITION_CONTROL, self.JP[j], maxVelocity = 0.5)
        j = j+1
    
      for i in  [9, 11, 13]:
        p.setJointMotorControl2(self.jacoId, i, p.POSITION_CONTROL, self.fing)

    else:
      j = 0
      for i in self.jacoArmJoints:
        p.resetJointState(self.jacoId, i, self.jointPoses[j])
        j = j+1 

    if (self.useRealTimeSimulation):
      dt = datetime.now()
      t = (dt.second / 60.) * 2. * math.pi
    else:
      t = t + 0.01

    if (self.useSimulation and self.useRealTimeSimulation == 0):
      p.stepSimulation()

    self.newPosInput = 0
    
    # p.addUserDebugLine(self.pos, self.pos2, [1,0,0,], 3, 0.3)

    p1 = [self.pos[0] - .02, self.pos[1], 0.002]
    p2 = [self.pos[0] + .02, self.pos[1], 0.002]
    p3 = [self.pos[0], self.pos[1] + .02, 0.002]
    p4 = [self.pos[0], self.pos[1] - .02, 0.002]

    p.addUserDebugLine([self.pos[0], self.pos[1], 0.001], [self.pos2[0], self.pos2[1], 0], [1,0,0,], 8, self.bciRate)

    p.addUserDebugLine(p1, p2, [0,1,0], 6, self.bciRate)
    p.addUserDebugLine(p3, p4, [0,1,0], 6, self.bciRate)
