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
  def __init__(self, mode,angle, dl):
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
    
    self.cubeStep = 3
    self.hasCube = False
    self.dist = 0.0125
    self.dist = 0.0075
    self.distf = .0338
    self.debuglen = 0.3

    self.wu = [0, 0.57, 0.45]
    self.wl = [-0.6, 0, 0.03]
    self.fu = 1.35
    self.fl = 0.
    self.center = np.array([-0.35, 0.3, 0.25])
    self.bciRate = 0.125
    self.mode = mode
    self.angle  =angle
    self.dl  = dl
    # self.distPos = [0.25, 0.05, -0.15, 0.35]
    self.distPos = [-0.10, -.3, -0.5, -0.7]

    p.loadURDF("plane.urdf",[0,0,-.65])
    p.loadURDF("table/table.urdf", basePosition=[-0.6,0.45,-0.65])
    self.jacoId = p.loadURDF("URDFs/jaco/j2n6s300.urdf", [0,0,0],  useFixedBase=True)

    p.resetBasePositionAndOrientation(self.jacoId,[0,0,0],[0,0,0,1])
    # p.resetDebugVisualizerCamera(cameraDistance=0.8, cameraYaw=0, cameraPitch=-89.99, cameraTargetPosition=[-0.35,-0.3,0.0])

    if self.mode == 0 or self.mode == 2: 
      if self.angle == 0:
        p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw=0, cameraPitch=-30, cameraTargetPosition=[-0.35,0.3,0.1])
      else:
        p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw=-45, cameraPitch=-10, cameraTargetPosition=[-0.35,0.3,0.1])
    elif self.mode == 3  or self.mode == 4:
      p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw= 30, cameraPitch=-45, cameraTargetPosition=[-0.35,0.3,0.1])
      # p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw= 0, cameraPitch=-45, cameraTargetPosition=[-0.35,0.3,0.1])
    else: 
      p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw=0, cameraPitch=-30, cameraTargetPosition=[-0.35,0.3,0.1])
    # p.resetDebugVisualizerCamera(cameraDistance=0.8, cameraYaw=30, cameraPitch=-30, cameraTargetPosition=[-0.35,-0.3,0.0])
    
    self.jacoEndEffectorIndex = 8
    self.jacoArmJoints = [2, 3, 4, 5, 6, 7]
    self.jacofingerJoints = [9, 11, 13]
    self.jacoJoints = [2, 3, 4, 5, 6, 7, 9, 11, 13]

    self.jd = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

    self.cube1Id = p.loadURDF("cube_small.urdf",[-1., -1., -1.], [0,0,0,1])
    self.reset()

  def set_block_pos(self, pos):
    pos[0] = self.center[0] + pos[0]
    pos[1] = self.center[1] + pos[1]
    p.resetBasePositionAndOrientation(self.cube1Id, [pos[0], pos[1], 0], [0,0,0,1])


  def drawDistLines(self):
    y1 = -0.05
    y2 = 0.95
    lw = 3
    c = [0,0,1]
    z = -.02

    dist = 0.00
    delta = 0.20

    while dist > -1.0:
      c1 = [dist, y1, z] 
      c2 = [dist, y2, z]

      p.addUserDebugLine(c1, c2, c, lw, 0)
      dist = dist - delta

  def updateDist(self, ds, c):
    x = self.distPos[ds-1];
    z = -.02
    lw = 10

    d = .06;

    pos = [x, 0.2, z]

    c1 = [pos[0] - d, pos[1] - d, z]
    c2 = [pos[0] - d, pos[1] + d, z]
    c3 = [pos[0] + d, pos[1] + d, z]
    c4 = [pos[0] + d, pos[1] - d, z]

    p.addUserDebugLine(c1, c2, c, lw, self.bciRate)
    p.addUserDebugLine(c2, c3, c, lw, self.bciRate)
    p.addUserDebugLine(c3, c4, c, lw, self.bciRate)
    p.addUserDebugLine(c4, c1, c, lw, self.bciRate)

  def drawAxes(self):
    c1 = [0, 0, -0.2] + self.center
    c2 = [0, 0, 0.2] + self.center
    c3 = [0, 0.2, 0.0] + self.center
    c4 = [0, -0.2, 0.0] + self.center
    c5 = [0.2, 0.0, 0.0] + self.center
    c6 = [-0.2, 0.0, 0.0] + self.center

    lw = 2

    p.addUserDebugLine(c1, c2, [0,0,0], lw, 0)
    p.addUserDebugLine(c3, c4, [0,0,0], lw, 0)
    p.addUserDebugLine(c5, c6, [0,0,0], lw, 0)


  def set_cubeTarget(self, pos, c):

      pos[0] = self.center[0] + pos[0]
      pos[1] = self.center[1] + pos[1]
      pos[2] = self.center[2] + pos[2]
      lw = 6
      d = .05
      # z = .1
      self.c1 = [pos[0] - d, pos[1]-d, pos[2]-d]
      self.c2 = [pos[0] + d, pos[1]-d, pos[2]-d]
      self.c3 = [pos[0] + d, pos[1]-d, pos[2]+d]
      self.c4 = [pos[0] - d, pos[1]-d, pos[2]+d]
      self.c5 = [pos[0] - d, pos[1]+d, pos[2]-d]
      self.c6 = [pos[0] + d, pos[1]+d, pos[2]-d]
      self.c7 = [pos[0] + d, pos[1]+d, pos[2]+d]
      self.c8 = [pos[0] - d, pos[1]+d, pos[2]+d]

      self.l1 = p.addUserDebugLine(self.c1, self.c2, c, 6, 0)
      self.l2 = p.addUserDebugLine(self.c2, self.c3, c, 6, 0)
      self.l3 = p.addUserDebugLine(self.c3, self.c4, c, 6, 0)
      self.l4 = p.addUserDebugLine(self.c4, self.c1, c, 6, 0)

      self.l5 = p.addUserDebugLine(self.c5, self.c6, c, 6, 0)
      self.l6 = p.addUserDebugLine(self.c6, self.c7, c, 6, 0)
      self.l7 = p.addUserDebugLine(self.c7, self.c8, c, 6, 0)
      self.l8 = p.addUserDebugLine(self.c8, self.c5, c, 6, 0)

      self.l9 = p.addUserDebugLine(self.c1, self.c5, c, 6, 0)
      self.l10 = p.addUserDebugLine(self.c2, self.c6, c, 6, 0)
      self.l11 = p.addUserDebugLine(self.c3, self.c7, c, 6, 0)
      self.l12 = p.addUserDebugLine(self.c4, self.c8, c, 6, 0)

      self.newPosInput = 1
      self.inverseKin()
      for i in self.jacoArmJoints:
        p.resetJointState(self.jacoId,i, self.JP[i-2])

      for i in  [9, 11, 13]:
        p.setJointMotorControl2(self.jacoId, i, p.POSITION_CONTROL, self.fing)


  def set_cubeColor(self, pos, c):
    lw = 16

    self.l1 = p.addUserDebugLine(self.c1, self.c2, c, lw, 0)
    self.l2 = p.addUserDebugLine(self.c2, self.c3, c, lw, 0)
    self.l3 = p.addUserDebugLine(self.c3, self.c4, c, lw, 0)
    self.l4 = p.addUserDebugLine(self.c4, self.c1, c, lw, 0)

    self.l5 = p.addUserDebugLine(self.c5, self.c6, c, lw, 0)
    self.l6 = p.addUserDebugLine(self.c6, self.c7, c, lw, 0)
    self.l7 = p.addUserDebugLine(self.c7, self.c8, c, lw, 0)
    self.l8 = p.addUserDebugLine(self.c8, self.c5, c, lw, 0)

    self.l9 = p.addUserDebugLine(self.c1, self.c5, c, lw, 0)
    self.l10 = p.addUserDebugLine(self.c2, self.c6, c, lw, 0)
    self.l11 = p.addUserDebugLine(self.c3, self.c7, c, lw, 0)
    self.l12 = p.addUserDebugLine(self.c4, self.c8, c, lw, 0)


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
    # p.resetSimulation()
    p.resetBasePositionAndOrientation(self.cube1Id, [-1., -1., -1.], [0,0,0,1])
    p.removeAllUserDebugItems()
    # self.drawAxes()
    self.drawDistLines()
    rp = [0,math.pi/4,math.pi,1.0*math.pi, 1.8*math.pi, 0*math.pi, 1.75*math.pi, 0.5*math.pi]


    self.pos =list([-0.0, 0.2, 0.2])

    self.orn = p.getQuaternionFromEuler([0,math.pi,math.pi/2])
    self.fing = 0.675

    for i in range(8):
      p.resetJointState(self.jacoId,i, rp[i])
    
    self.newPosInput = 1
    self.inverseKin()
    for i in self.jacoArmJoints:
      p.resetJointState(self.jacoId,i, self.JP[i-2])

    for i in  [9, 11, 13]:
      p.setJointMotorControl2(self.jacoId, i, p.POSITION_CONTROL, self.fing)

    ls = p.getLinkState(self.jacoId, self.jacoEndEffectorIndex)
    p.setRealTimeSimulation(self.useRealTimeSimulation)
    p.configureDebugVisualizer(p.COV_ENABLE_RENDERING,1) 

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

      for i in  [9, 11, 13]:
        p.setJointMotorControl2(self.jacoId, i, p.POSITION_CONTROL, self.fing)


    if (self.useRealTimeSimulation):
      dt = datetime.now()
      t = (dt.second / 60.) * 2. * math.pi
    else:
      t = t + 0.01

    if (self.useSimulation and self.useRealTimeSimulation == 0):
      p.stepSimulation()

    self.newPosInput = 0

    p1 = [self.pos[0] - .02, self.pos[1], 0.002]
    p2 = [self.pos[0] + .02, self.pos[1], 0.002]
    p3 = [self.pos[0], self.pos[1] + .02, 0.002]
    p4 = [self.pos[0], self.pos[1] - .02, 0.002]

    if self.dl:
      if self.mode == 0:
        p.addUserDebugLine([self.pos[0], self.pos[1], 0.001], [self.pos2[0], self.pos2[1], 0], [1,0,0,], 8, self.bciRate)
      elif self.mode == 2:
        p.addUserDebugLine([self.pos[0], self.pos[1], 0.001], [self.pos2[0], self.pos2[1], 0], [1,0,0,], 8, self.bciRate)
      elif self.mode ==1:
        p.addUserDebugLine([self.pos[0], self.pos[1], self.pos[2] + 0.05], [self.pos2[0], self.pos2[1], self.pos2[2] + .05], [1,0,0,], 8, self.bciRate)
      elif self.mode == 3 or self.mode == 4:
        # if self.key == 1 or self.key == 7:
        p.addUserDebugLine([self.pos[0], self.pos[1], self.pos[2] + 0.05], [self.pos2[0], self.pos2[1], self.pos2[2] + .05], [1,0,0,], 8, self.bciRate)
        # else:
        #   p.addUserDebugLine([self.pos[0], self.pos[1], 0.001], [self.pos2[0], self.pos2[1], 0], [1,0,0,], 8, self.bciRate)