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

from scipy.interpolate import interp1d

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

    # self.wu = [0, 0.57, 0.45]
    # self.wl = [-0.6, 0, 0.03]
    self.fu = 1.35
    self.fl = 0.
    self.center = np.array([-0.35, 0.3, 0.25])
    self.bciRate = 0.125
    self.mode = mode
    self.angle  =angle
    self.dl  = dl
    self.tp = [0,0,0]

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
      p.resetDebugVisualizerCamera(cameraDistance=0.25, cameraYaw= 30, cameraPitch=-45, cameraTargetPosition=[-0.25,0.3,0.1])
    else: 
      p.resetDebugVisualizerCamera(cameraDistance=0.3, cameraYaw=0, cameraPitch=-30, cameraTargetPosition=[-0.25,0.3,0.1])
    # p.resetDebugVisualizerCamera(cameraDistance=0.8, cameraYaw=30, cameraPitch=-30, cameraTargetPosition=[-0.35,-0.3,0.0])
    
    self.jacoEndEffectorIndex = 8
    self.jacoArmJoints = [2, 3, 4, 5, 6, 7]
    self.jacofingerJoints = [9, 11, 13]
    self.jacoJoints = [2, 3, 4, 5, 6, 7, 9, 11, 13]

    self.jd = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

    self.cube1Id = p.loadURDF("cube_small.urdf",[-1., -1., -1.], [0,0,0,1])
    self.cube2Id = p.loadURDF("cube_small.urdf",[-1., -1., -1.], [0,0,0,1])

    # self.cube2Id = p.loadURDF("cube_small.urdf",[-0.3, 0.25, 0.0], [0,0,0,1])
    # self.cube3Id = p.loadURDF("cube_small.urdf",[-0.1, 0.2, 0.0], [0,0,0,1])
    self.reset()


  def set_block_pos(self, pos, target):
    pos[0] = self.center[0] + pos[0]
    pos[1] = self.center[1] + pos[1]

    d = .08

    p.resetBasePositionAndOrientation(self.cube1Id, [pos[0], pos[1], 0], [0,0,0,1])
    p.resetBasePositionAndOrientation(self.cube2Id, [pos[0], pos[1], 0.05], [0,0,0,1])
    c1 = [pos[0] - d, pos[1] - d, 0.0]
    c2 = [pos[0] - d, pos[1] + d, 0.0]
    c3 = [pos[0] + d, pos[1] + d, 0.0]
    c4 = [pos[0] + d, pos[1] - d, 0.0]

    
    # if self.dl:
    # Green x
    self.l1 = p.addUserDebugLine(c1, c2, [0,0,1], 3, 0)
    self.l2 = p.addUserDebugLine(c2, c3, [0,0,1], 3, 0)
    self.l3 = p.addUserDebugLine(c3, c4, [0,0,1], 3, 0)
    self.l4 = p.addUserDebugLine(c4, c1, [0,0,1], 3, 0)
    # p.resetBasePositionAndOrientation(self.cube1Id, [pos[0], pos[1], 0], [0,0,0,1])

    n = 20
    a = np.array([-0.15, 0.1])
    b = np.array([pos[0], pos[1]])
    mid  = (b - a)
    c = np.array([mid[0] - mid[1]/2, mid[1], + mid[0]/2])
    
    l = np.linalg.norm(mid)

    x = np.array([0, l/2,l])
    y = np.array([0 , -l/3, 0])

    f = interp1d(x, y, kind='quadratic')


    uv = mid/l

    X = np.linspace(0, l, num = n, endpoint=True)

    unitVec = uv

    c1 = [-0.15, 0.1,0]
    j = 0

    self.tp = pos


  def draw_bound(self, pos):
    col = [0,0,0]
    pos[0] = self.center[0] + pos[0]
    pos[1] = self.center[1] + pos[1]

    d = .08

    c1 = [pos[0] - d, pos[1] - d, 0.0]
    c2 = [pos[0] - d, pos[1] + d, 0.0]
    c3 = [pos[0] + d, pos[1] + d, 0.0]
    c4 = [pos[0] + d, pos[1] - d, 0.0]

    
    # if self.dl:
    # Green x
    self.h1 = p.addUserDebugLine(c1, c2, [0,0,1], 3, 0)
    self.h2 = p.addUserDebugLine(c2, c3, [0,0,1], 3, 0)
    self.h3 = p.addUserDebugLine(c3, c4, [0,0,1], 3, 0)
    self.h4 = p.addUserDebugLine(c4, c1, [0,0,1], 3, 0)

  def update_bound(self, pos, c):
    if c == 1:
      col = [1,0,0]
    else:
      col = [0,1,0] 
    pos[0] = self.center[0] + pos[0]
    pos[1] = self.center[1] + pos[1]

    d = .08
    
    c1 = [pos[0] - d, pos[1] - d, 0.0]
    c2 = [pos[0] - d, pos[1] + d, 0.0]
    c3 = [pos[0] + d, pos[1] + d, 0.0]
    c4 = [pos[0] + d, pos[1] - d, 0.0]
    
    # if self.dl:
    # Green x
    self.h1 = p.addUserDebugLine(c1, c2, col, 3, 0,replaceItemUniqueId=self.h1)
    self.h2 = p.addUserDebugLine(c2, c3, col, 3, 0, replaceItemUniqueId=self.h2)
    self.h3 = p.addUserDebugLine(c3, c4, col, 3, 0, replaceItemUniqueId=self.h3)
    self.h4 = p.addUserDebugLine(c4, c1, col, 3, 0,replaceItemUniqueId=self.h4)

  def set_bound_color(self, pos, c):

    if c == 1:

      col = [1,0,0]
    else:
      col = [0,1,0] 
    pos[0] = self.center[0] + pos[0]
    pos[1] = self.center[1] + pos[1]

    d = .08

    p.resetBasePositionAndOrientation(self.cube1Id, [pos[0], pos[1], 0], [0,0,0,1])
    c1 = [pos[0] - d, pos[1] - d, 0.0]
    c2 = [pos[0] - d, pos[1] + d, 0.0]
    c3 = [pos[0] + d, pos[1] + d, 0.0]
    c4 = [pos[0] + d, pos[1] - d, 0.0]

    
    # if self.dl:
    # Green x
    self.l1 = p.addUserDebugLine(c1, c2, col, 3, 0,replaceItemUniqueId=self.l1)
    self.l2 = p.addUserDebugLine(c2, c3, col, 3, 0, replaceItemUniqueId=self.l2)
    self.l3 = p.addUserDebugLine(c3, c4, col, 3, 0, replaceItemUniqueId=self.l3)
    self.l4 = p.addUserDebugLine(c4, c1, col, 3, 0,replaceItemUniqueId=self.l4)
    # p.resetBasePositionAndOrientation(self.cube1Id, [pos[0], pos[1], 0], [0,0,0,1])


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

      d1 = [pos[0], pos[1], pos[2]]
      d2 = [self.center[0],self.center[1] ,self.center[2] ]
      self.l13 = p.addUserDebugLine(d1, d2, [0,0,0], 4, 0)

      self.tp = pos

  def set_cubeColor(self, pos, c, lw):
    # lw = 16
    # print(c)

    self.l1 = p.addUserDebugLine(self.c1, self.c2, c, lw, 0, replaceItemUniqueId=self.l1)
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



  def updateCommand(self, key):
    self.key = key
    Rrm = R.from_quat(self.orn)
    Rnew =  Rrm.as_dcm() 

    baseTheta = self.JP[0]
    s = math.sin(baseTheta)
    c = math.cos(baseTheta)

    n = np.sqrt(self.pos[0]*self.pos[0] + self.pos[1]*self.pos[1])
    dx = -self.pos[1]/n
    dy = self.pos[0]/n

    if 0 < key < 30:
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
    elif key == 7:
      self.pos[2] = self.pos[2] + self.dist
      self.pos2 = [self.pos[0], self.pos[1], self.pos[2] + self.debuglen]
    elif key == 1:
      self.pos[2] = self.pos[2] - self.dist
      self.pos2 = [self.pos[0], self.pos[1], self.pos[2] - self.debuglen]
    elif key == 16:
      self.fing = self.fing - self.distf
      self.pos2 = [self.pos[0] + self.debuglen, self.pos[1], self.pos[2]]
      if self.dl: 
        p.addUserDebugText('OPEN', (self.pos2[0] + 0.05, self.pos2[1], self.pos2[2] + .01),  [0,0,0], 4, self.bciRate)
    elif key == 14: 
      self.fing = self.fing + self.distf
      self.pos2 = [self.pos[0] - self.debuglen, self.pos[1], self.pos[2]]
      if self.dl:  
        p.addUserDebugText('CLOSE', (self.pos2[0] - 0.3, self.pos2[1], self.pos2[2] + .01),  [0,0,0], 4, self.bciRate)
    elif key == 18:
      self.pos[2] = self.pos[2] + self.dist
      self.pos2 = [self.pos[0], self.pos[1], self.pos[2] + self.debuglen]
    elif key == 12:
      self.pos[2] = self.pos[2] - self.dist
      self.pos2 = [self.pos[0], self.pos[1], self.pos[2] - self.debuglen]
    elif key == 0:
      self.pos2 = self.pos
    elif key == 26:
      self.pos[0] = self.pos[0] - self.dist*dx
      self.pos[1] = self.pos[1] - self.dist*dy
      self.pos2 = [self.pos[0] - self.debuglen*dx, self.pos[1] - self.debuglen*dy]

    elif key == 24: 
      self.pos[0] = self.pos[0] + self.dist*dx
      self.pos[1] = self.pos[1] + self.dist*dy
      self.pos2 = [self.pos[0] + self.debuglen*dx, self.pos[1] + self.debuglen*dy]

    elif key == 28:
      self.pos[0] = self.pos[0] + self.dist*c
      self.pos[1] = self.pos[1] - self.dist*s
      self.pos2 = [self.pos[0] + self.debuglen*c, self.pos[1] - self.debuglen*s]

    elif key == 22:
      self.pos[0] = self.pos[0] - self.dist*c
      self.pos[1] = self.pos[1] + self.dist*s
      self.pos2 = [self.pos[0] - self.debuglen*c, self.pos[1] + self.debuglen*s]

    if key == 46:
      # self.pos[0] = self.pos[0] + self.dist
      self.pos2 = [self.pos[0] + self.debuglen, self.pos[1], self.pos[2]]
    elif key == 44: 
      # self.pos[0] = self.pos[0] - self.dist
      self.pos2 = [self.pos[0] - self.debuglen, self.pos[1], self.pos[2]]
    elif key == 48:
      # self.pos[1] = self.pos[1] + self.dist
      self.pos2 = [self.pos[0], self.pos[1] + self.debuglen, self.pos[2]]
    elif key == 42:
      # self.pos[1] = self.pos[1] - self.dist
      self.pos2 = [self.pos[0], self.pos[1] - self.debuglen, self.pos[2]]
    elif key == 47:
      # self.pos[2] = self.pos[2] + self.dist
      self.pos2 = [self.pos[0], self.pos[1], self.pos[2] + self.debuglen]
    elif key == 41:
      # self.pos[2] = self.pos[2] - self.dist
      self.pos2 = [self.pos[0], self.pos[1], self.pos[2] - self.debuglen]


    if self.fing > self.fu:
      self.fing =  self.fu
    if self.fing < self.fl:
      self.fing =  self.fl

  def set_robotPos(self, rp, key):
    self.pos[0] = self.center[0] + rp[0]
    self.pos[1] = self.center[1] + rp[1]
    self.pos[2] = self.center[2]  + rp[2]
    # print(key)
    # Position
    if key == 6:
      self.pos2 = [self.pos[0] + self.debuglen, self.pos[1], self.pos[2]]
    elif key == 4: 
      self.pos2 = [self.pos[0] - self.debuglen, self.pos[1], self.pos[2]]
    elif key == 8:
      self.pos2 = [self.pos[0], self.pos[1] + self.debuglen, self.pos[2]]
    elif key == 2:
      self.pos2 = [self.pos[0], self.pos[1] - self.debuglen, self.pos[2]]
    elif key == 7:
      self.pos2 = [self.pos[0], self.pos[1], self.pos[2] + self.debuglen]
    elif key == 1:
      self.pos2 = [self.pos[0], self.pos[1], self.pos[2] - self.debuglen]
    else:
      self.pos2 = [self.pos[0], self.pos[1], self.pos[2]] 
    
    if key == 100:

      c1 = [self.pos[0] + .02, self.pos[1], self.pos[2] - .05]
      c2 = [self.pos[0] - .02, self.pos[1], self.pos[2] - .05]
      c3 = [self.pos[0], self.pos[1] - .02, self.pos[2] - .05]
      c4 = [self.pos[0], self.pos[1] + .02, self.pos[2] - .05]

      p.addUserDebugLine(c1,c2, [0,1,1], 8, self.bciRate)
      p.addUserDebugLine(c3,c4, [0,1,1], 8, self.bciRate)



    Ez = np.arctan2(self.tp[1] -self.pos[1], self.tp[0] - self.pos[0])
    dist = np.sqrt((self.tp[1] -self.pos[1])*(self.tp[1] -self.pos[1]) + (self.tp[0] -self.pos[0])*(self.tp[0] -self.pos[0]))
    # print((self.tp[1] -self.pos[1]))
    # print((self.tp[0] -self.pos[0]))
    print(dist)
    # self.orn = p.getQuaternionFromEuler([math.pi/2,3*math.pi/2,3*math.pi/2 + 0.2])
    if dist > 0.05:
      self.orn = p.getQuaternionFromEuler([math.pi/2,3*math.pi/2,Ez + math.pi/2])

    self.newPosInput = 1

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
    # p.resetBasePositionAndOrientation(self.cube2Id, [-0.3, 0.25, 0.0], [0,0,0,1])
    # p.resetBasePositionAndOrientation(self.cube3Id, [-0.1, 0.2, 0.0], [0,0,0,1])
    p.removeAllUserDebugItems()
    # self.drawAxes()
    rp = [0,math.pi/4,math.pi,1.0*math.pi, 1.8*math.pi, 0*math.pi, 1.75*math.pi, 0.5*math.pi]

   


    self.pos =list([-0.15, 0.15, 0.2])
    self.pos =list([-0.4, 0.4, 0.05])

    # self.orn = p.getQuaternionFromEuler([0,math.pi,math.pi/2])
    self.orn = p.getQuaternionFromEuler([math.pi/2,3*math.pi/2,3*math.pi/2 + 0.2])
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
    
    # p.addUserDebugLine(self.pos, self.pos2, [1,0,0,], 3, 0.3)

    p1 = [self.pos[0] - .02, self.pos[1], 0.000]
    p2 = [self.pos[0] + .02, self.pos[1], 0.000]
    p3 = [self.pos[0], self.pos[1] + .02, 0.000]
    p4 = [self.pos[0], self.pos[1] - .02, 0.000]

    p.addUserDebugLine(p1, p2, [0,1,0], 6, self.bciRate)
    p.addUserDebugLine(p3, p4, [0,1,0], 6, self.bciRate)

    if self.dl:
      if self.mode == 0:
        p.addUserDebugLine([self.pos[0], self.pos[1], 0.001], [self.pos2[0], self.pos2[1], 0], [1,0,0,], 8, self.bciRate)
      elif self.mode == 2:
        p.addUserDebugLine([self.pos[0], self.pos[1], 0.001], [self.pos2[0], self.pos2[1], 0], [1,0,0,], 8, self.bciRate)
      elif self.mode ==1:
        p.addUserDebugLine([self.pos[0], self.pos[1], self.pos[2] + 0.05], [self.pos2[0], self.pos2[1], self.pos2[2] + .05], [1,0,0,], 8, self.bciRate)
      elif self.mode == 3 or self.mode == 4:
        # if self.key == 1 or self.key == 7:
        p.addUserDebugLine([self.pos[0], self.pos[1], self.pos[2] ], [self.pos2[0], self.pos2[1], self.pos2[2] ], [1,0,0,], 8, self.bciRate)
        # else:
        #   p.addUserDebugLine([self.pos[0], self.pos[1], 0.001], [self.pos2[0], self.pos2[1], 0], [1,0,0,], 8, self.bciRate)