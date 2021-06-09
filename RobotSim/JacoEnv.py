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

    # self.wu = [0, 0.57, 0.45]
    # self.wl = [-0.6, 0, 0.03]
    self.fu = 1.35
    self.fl = 0.
    self.center = np.array([-0.35, 0.3, 0.25])
    self.bciRate = 0.125
    self.mode = mode
    self.angle  =angle
    self.dl  = dl
    self.key = 0

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
    elif self.mode == 3  or self.mode == 4 or self.mode == 6 or self.mode == 8:
      p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw= 30, cameraPitch=-45, cameraTargetPosition=[-0.35,0.3,0.1])
    elif self.mode == 5:
      p.resetDebugVisualizerCamera(cameraDistance=0.4, cameraYaw= 30, cameraPitch=-45, cameraTargetPosition=[-0.35,0.3,0.1])
    else: 
      p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw=0, cameraPitch=-30, cameraTargetPosition=[-0.35,0.3,0.1])

    self.jacoEndEffectorIndex = 8
    self.jacoArmJoints = [2, 3, 4, 5, 6, 7]
    self.jacofingerJoints = [9, 11, 13]
    self.jacoJoints = [2, 3, 4, 5, 6, 7, 9, 11, 13]

    self.jd = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    # self.jd = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
        # self.jd = [01.01, 01.01, 01.01, 01.01, 01.01, 01.01, 01.01, 01.01, 01.01, 0.01]

    # self.cube1Id = p.loadURDF("cube_small.urdf",[0.2, 0, -0.2] + self.center, [0,0,0,1])
    # self.cube2Id = p.loadURDF("cube_small.urdf",[-0.2, 0., -0.2] + self.center, [0,0,0,1])

    # self.cube3Id = p.loadURDF("cube_small.urdf",[0, 0.2, -0.2] + self.center, [0,0,0,1])
    # self.cube4Id = p.loadURDF("cube_small.urdf",[-0, -0.2, -0.2] + self.center, [0,0,0,1])


    c  = [0, 1, 0]
    pos = np.array([0,0, -2])
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

    self.l13 = p.addUserDebugLine([0,0,0], [0,0,0], [0,0,0], 4, 0)

    self.LetterMode  = 0
    self.reset()
    self.newPosInput = 1
    self.inverseKin()
    for i in self.jacoArmJoints:
      p.resetJointState(self.jacoId,i, self.JP[i-2])

    for i in  [9, 11, 13]:
      p.setJointMotorControl2(self.jacoId, i, p.POSITION_CONTROL, self.fing)

    ls = p.getLinkState(self.jacoId, self.jacoEndEffectorIndex)
    p.setRealTimeSimulation(self.useRealTimeSimulation)
    p.configureDebugVisualizer(p.COV_ENABLE_RENDERING,1) 
    

  def write_letter(self, dpos, c, targetID):

    pos = dpos.copy()
    pos[0] = self.center[0] + pos[0]
    pos[1] = self.center[1] + pos[1]
    pos[2] = self.center[2] + pos[2]

    # print(letter)

    letters = ["A", "F", "B", "E", "C", "D"]
    offset = [[0,0,0], [0, -0.05, -0.05],  [0, -0.05, -0.05], [-.1, -0.05, -0.1], [0,0,-0.1], [0.1,-0.2,0.02], [0.05,-0.15,-0.1]]

    letter = letters[targetID-1]
    pos = pos + offset[targetID]
    print(pos)
    # p.addUserDebugText(letter, ([-0.15, 0.3, 0.15]),  [0,0,0], 10)
    # p.addUserDebugText(letter, (pos),  c, 12,replaceItemUniqueId=self.letter)
    self.letter = p.addUserDebugText(letter, (pos),  c, 12)
    self.LetterMode = 1

  def set_letterColor(self, dpos, c, targetID):

    pos = dpos.copy()
    pos[0] = self.center[0] + pos[0]
    pos[1] = self.center[1] + pos[1]
    pos[2] = self.center[2] + pos[2]

    # print(letter)

    letters = ["A", "F", "B", "E", "C", "D"]
    offset = [[0,0,0], [0, -0.05, -0.0],  [0, -0.05, -0.05], [-.1, -0.05, -0.1], [0,0,-0.1], [0.1,-0.2,0.02], [0.05,-0.15,-0.1]]

    letter = letters[targetID-1]
    pos = pos + offset[targetID]
    print(pos)
    # p.addUserDebugText(letter, ([-0.15, 0.3, 0.15]),  [0,0,0], 10)
    p.addUserDebugText(letter, (pos),  c, 12,replaceItemUniqueId=self.letter)

  def set_block_pos(self, pos, target):
    pos[0] = self.center[0] + pos[0]
    pos[1] = self.center[1] + pos[1]

    d = .05

    # p.resetBasePositionAndOrientation(self.cube1Id, [pos[0], pos[1], 0], [0,0,0,1])
    c1 = [pos[0] - d, pos[1] - d, 0.0]
    c2 = [pos[0] - d, pos[1] + d, 0.0]
    c3 = [pos[0] + d, pos[1] + d, 0.0]
    c4 = [pos[0] + d, pos[1] - d, 0.0]

    # if self.dl:
    # Green x
    self.l1 = p.addUserDebugLine(c1, c2, [0,1,0], 6, 0)
    self.l2 = p.addUserDebugLine(c2, c3, [0,1,0], 6, 0)
    self.l3 = p.addUserDebugLine(c3, c4, [0,1,0], 6, 0)
    self.l4 = p.addUserDebugLine(c4, c1, [0,1,0], 6, 0)

  def set_bound_color(self, pos, c):

    if c == 1:
      col = [1,0,0]
    elif c == 2:
      col = [0,0,1]
    else:
      col = [0,1,0] 
    pos[0] = self.center[0] + pos[0]
    pos[1] = self.center[1] + pos[1]

    d = .05

    # p.resetBasePositionAndOrientation(self.cube1Id, [pos[0], pos[1], 0], [0,0,0,1])
    c1 = [pos[0] - d, pos[1] - d, 0.0]
    c2 = [pos[0] - d, pos[1] + d, 0.0]
    c3 = [pos[0] + d, pos[1] + d, 0.0]
    c4 = [pos[0] + d, pos[1] - d, 0.0]

    
    # if self.dl:
    # Green x
    self.l1 = p.addUserDebugLine(c1, c2, col, 6, 0,replaceItemUniqueId=self.l1)
    self.l2 = p.addUserDebugLine(c2, c3, col, 6, 0, replaceItemUniqueId=self.l2)
    self.l3 = p.addUserDebugLine(c3, c4, col, 6, 0, replaceItemUniqueId=self.l3)
    self.l4 = p.addUserDebugLine(c4, c1, col, 6, 0,replaceItemUniqueId=self.l4)
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

      p.addUserDebugLine(self.c1, self.c2, c, lw, 0, replaceItemUniqueId=self.l1)
      p.addUserDebugLine(self.c2, self.c3, c, lw, 0, replaceItemUniqueId=self.l2)
      p.addUserDebugLine(self.c3, self.c4, c, lw, 0, replaceItemUniqueId=self.l3)
      p.addUserDebugLine(self.c4, self.c1, c, lw, 0, replaceItemUniqueId=self.l4)

      p.addUserDebugLine(self.c5, self.c6, c, lw, 0, replaceItemUniqueId=self.l5)
      p.addUserDebugLine(self.c6, self.c7, c, lw, 0, replaceItemUniqueId=self.l6)
      p.addUserDebugLine(self.c7, self.c8, c, lw, 0, replaceItemUniqueId=self.l7)
      p.addUserDebugLine(self.c8, self.c5, c, lw, 0, replaceItemUniqueId=self.l8)

      p.addUserDebugLine(self.c1, self.c5, c, lw, 0, replaceItemUniqueId=self.l9)
      p.addUserDebugLine(self.c2, self.c6, c, lw, 0, replaceItemUniqueId=self.l10)
      p.addUserDebugLine(self.c3, self.c7, c, lw, 0, replaceItemUniqueId=self.l11)
      p.addUserDebugLine(self.c4, self.c8, c, lw, 0, replaceItemUniqueId=self.l12)

      d1 = [pos[0], pos[1], pos[2]]
      d2 = [self.center[0],self.center[1] ,self.center[2] ]

  def set_cubeColor(self, pos, c, lw):
    p.addUserDebugLine(self.c1, self.c2, c, lw, 0, replaceItemUniqueId=self.l1)
    p.addUserDebugLine(self.c2, self.c3, c, lw, 0, replaceItemUniqueId=self.l2)
    p.addUserDebugLine(self.c3, self.c4, c, lw, 0, replaceItemUniqueId=self.l3)
    p.addUserDebugLine(self.c4, self.c1, c, lw, 0, replaceItemUniqueId=self.l4)

    p.addUserDebugLine(self.c5, self.c6, c, lw, 0, replaceItemUniqueId=self.l5)
    p.addUserDebugLine(self.c6, self.c7, c, lw, 0, replaceItemUniqueId=self.l6)
    p.addUserDebugLine(self.c7, self.c8, c, lw, 0, replaceItemUniqueId=self.l7)
    p.addUserDebugLine(self.c8, self.c5, c, lw, 0, replaceItemUniqueId=self.l8)

    p.addUserDebugLine(self.c1, self.c5, c, lw, 0, replaceItemUniqueId=self.l9)
    p.addUserDebugLine(self.c2, self.c6, c, lw, 0, replaceItemUniqueId=self.l10)
    p.addUserDebugLine(self.c3, self.c7, c, lw, 0, replaceItemUniqueId=self.l11)
    p.addUserDebugLine(self.c4, self.c8, c, lw, 0, replaceItemUniqueId=self.l12)

  def grabCube(self):
    cp = p.getBasePositionAndOrientation(self.cube1Id)
    cpos = cp[0]
    self.set_robotPos([cpos[0] - self.center[0], cpos[1] - self.center[1], cpos[2] - self.center[2] + .02] , 0)

    self.setFing(1.05)

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

  def setFing(self, fp):
    if self.mode == 5:
      if self.TargetID >6:
        self.fing = fp
      else:

        orn = p.getEulerFromQuaternion(self.orn)
        ornNew = [math.pi, fp, 0]
        self.orn = p.getQuaternionFromEuler(ornNew)
        self.newPosInput = 1
    else:
      self.fing = fp


  def displayCue(self, cue, c):
    p.removeAllUserDebugItems()

    if c == 0:
      color = [0,0,0]
    elif c == 1:
      color = [0,0,1]
    elif c == 2:
      color = [0, 1, 0]
    elif c == 3:
      color = [1, 0, 0]

    if cue == 7:
        p.addUserDebugText('OPEN', [-0.45,0.3,0.01],  color, 12)
    elif cue == 8:
        p.addUserDebugText('CLOSE', [-0.47,0.3,0.01],  color, 12)
    elif cue == 3:
        p.addUserDebugText('UNITED STATES', [-0.65,0.3,0.01],  color, 12)
        self.orn = p.getQuaternionFromEuler([math.pi,0,0])
    elif cue == 6:
        p.addUserDebugText('SOUTH AMERICA', [-0.65,0.3,0.01], color, 12)
        self.orn = p.getQuaternionFromEuler([math.pi,math.pi/2,0])

    self.TargetID = cue
    self.newPosInput = 1
    self.step()
    self.newPosInput = 1
    self.step()


  def set_robotPos(self, rp, key):
    self.key = key;
    self.pos[0] = self.center[0] + rp[0]
    self.pos[1] = self.center[1] + rp[1]
    self.pos[2] = self.center[2]  + rp[2]

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
    
    # if key == 100:

    #   c1 = [self.pos[0] + .02, self.pos[1], self.pos[2] - .05]
    #   c2 = [self.pos[0] - .02, self.pos[1], self.pos[2] - .05]
    #   c3 = [self.pos[0], self.pos[1] - .02, self.pos[2] - .05]
    #   c4 = [self.pos[0], self.pos[1] + .02, self.pos[2] - .05]

    #   p.addUserDebugLine(c1,c2, [0,1,1], 8, self.bciRate)
    #   p.addUserDebugLine(c3,c4, [0,1,1], 8, self.bciRate)

    self.newPosInput = 1

  def inverseKin(self):
    if (self.newPosInput == 1):
      self.jointPoses = p.calculateInverseKinematics(self.jacoId,
                                                self.jacoEndEffectorIndex,
                                                self.pos,
                                                self.orn,
                                                jointDamping=self.jd,
                                                solver=self.ikSolver,
                                                maxNumIterations=1000,
                                                residualThreshold=.001)

      self.JP = list(self.jointPoses)

  def reset(self):
    # p.resetSimulation()
    # p.resetBasePositionAndOrientation(self.cube1Id, [-1., -1., -1.], [0,0,0,1])
    
    if self.LetterMode:
      p.removeAllUserDebugItems()

    if self.mode == 0 or self.mode == 2: 
      if self.angle == 0:
        p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw=0, cameraPitch=-30, cameraTargetPosition=[-0.35,0.3,0.1])
      else:
        p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw=-45, cameraPitch=-10, cameraTargetPosition=[-0.35,0.3,0.1])
    elif self.mode == 3  or self.mode == 4 or self.mode == 8:
      p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw= 30, cameraPitch=-45, cameraTargetPosition=[-0.35,0.3,0.1])
    elif self.mode == 10:
      p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw= 30, cameraPitch=-45, cameraTargetPosition=[-0.35,0.3,0.1])
    elif self.mode == 5:
      p.resetDebugVisualizerCamera(cameraDistance=0.25, cameraYaw= 0, cameraPitch=0, cameraTargetPosition=[-0.35,0.3,0.2])  
    elif self.mode == 7:
      p.resetDebugVisualizerCamera(cameraDistance=0.5, cameraYaw= 25, cameraPitch=-20, cameraTargetPosition=[-0.35,0.3,0.2])    
    else: 
      p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw=0, cameraPitch=-0, cameraTargetPosition=[-0.35,0.3,0.1])

    pos = np.array([0,0, -2])
    lw = 6
    d = .05
    c = [0,0,1]
    # z = .1
    self.c1 = [pos[0] - d, pos[1]-d, pos[2]-d]
    self.c2 = [pos[0] + d, pos[1]-d, pos[2]-d]
    self.c3 = [pos[0] + d, pos[1]-d, pos[2]+d]
    self.c4 = [pos[0] - d, pos[1]-d, pos[2]+d]
    self.c5 = [pos[0] - d, pos[1]+d, pos[2]-d]
    self.c6 = [pos[0] + d, pos[1]+d, pos[2]-d]
    self.c7 = [pos[0] + d, pos[1]+d, pos[2]+d]
    self.c8 = [pos[0] - d, pos[1]+d, pos[2]+d]

    p.addUserDebugLine(self.c1, self.c2, c, lw, 0, replaceItemUniqueId=self.l1)
    p.addUserDebugLine(self.c2, self.c3, c, lw, 0, replaceItemUniqueId=self.l2)
    p.addUserDebugLine(self.c3, self.c4, c, lw, 0, replaceItemUniqueId=self.l3)
    p.addUserDebugLine(self.c4, self.c1, c, lw, 0, replaceItemUniqueId=self.l4)

    p.addUserDebugLine(self.c5, self.c6, c, lw, 0, replaceItemUniqueId=self.l5)
    p.addUserDebugLine(self.c6, self.c7, c, lw, 0, replaceItemUniqueId=self.l6)
    p.addUserDebugLine(self.c7, self.c8, c, lw, 0, replaceItemUniqueId=self.l7)
    p.addUserDebugLine(self.c8, self.c5, c, lw, 0, replaceItemUniqueId=self.l8)

    p.addUserDebugLine(self.c1, self.c5, c, lw, 0, replaceItemUniqueId=self.l9)
    p.addUserDebugLine(self.c2, self.c6, c, lw, 0, replaceItemUniqueId=self.l10)
    p.addUserDebugLine(self.c3, self.c7, c, lw, 0, replaceItemUniqueId=self.l11)
    p.addUserDebugLine(self.c4, self.c8, c, lw, 0, replaceItemUniqueId=self.l12)


    d1 = [0,0,0]
    d2 = [0,0,0]
    p.addUserDebugLine(d1, d2, [0,0,0], 4, 0,replaceItemUniqueId=self.l13)

    if self.mode < 5:
        self.drawAxes()
    rp = [0,math.pi/4,math.pi,1.0*math.pi, 1.8*math.pi, 0*math.pi, 1.75*math.pi, 0.5*math.pi]

    if self.mode == 0:
      self.pos =list([-0.35, 0.3, 0.2])
    elif self.mode ==3:
      self.pos =list([-0.35, 0.3, 0.25])
    else:
      self.pos =list([-0.35, 0.3, 0.2])

    self.orn = p.getQuaternionFromEuler([math.pi,0,math.pi/2])
    if self.mode == 7:
      self.fing = 0.0
    else:
      self.fing = 0.675;

    # for i in range(8):
    #   p.resetJointState(self.jacoId,i, rp[i])
    
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

    p1 = [self.pos[0] - .02, self.pos[1], 0.002]
    p2 = [self.pos[0] + .02, self.pos[1], 0.002]
    p3 = [self.pos[0], self.pos[1] + .02, 0.002]
    p4 = [self.pos[0], self.pos[1] - .02, 0.002]

    if self.mode == 7:
      p.addUserDebugLine(p1, p2, [0,1,0], 6, self.bciRate)
      p.addUserDebugLine(p3, p4, [0,1,0], 6, self.bciRate) 

    if self.dl:
      if self.mode == 0:
        p.addUserDebugLine([self.pos[0], self.pos[1], 0.001], [self.pos2[0], self.pos2[1], 0], [1,0,0,], 8, self.bciRate)
      elif self.mode == 2:
        p.addUserDebugLine([self.pos[0], self.pos[1], 0.001], [self.pos2[0], self.pos2[1], 0], [1,0,0,], 8, self.bciRate)
      elif self.mode ==1:
        p.addUserDebugLine([self.pos[0], self.pos[1], self.pos[2] + 0.05], [self.pos2[0], self.pos2[1], self.pos2[2] + .05], [1,0,0,], 8, self.bciRate)
      elif self.mode == 3 or self.mode == 4 or self.mode == 6 or self.mode == 7 or self.mode == 8:
        if self.key == 100:
          p1 = [self.pos[0] - .15, self.pos[1], self.pos2[2]]
          p2 = [self.pos[0] + .15, self.pos[1], self.pos2[2]]
          p3 = [self.pos[0], self.pos[1] + .15, self.pos2[2]]
          p4 = [self.pos[0], self.pos[1] - .15, self.pos2[2]]

          p.addUserDebugLine(p1, p2, [0,1,1], 8, self.bciRate)
          p.addUserDebugLine(p3, p4, [0,1,1], 8, self.bciRate)
        else:
          p.addUserDebugLine([self.pos[0], self.pos[1], self.pos[2] + 0.05], [self.pos2[0], self.pos2[1], self.pos2[2] + .05], [1,0,0,], 8, self.bciRate)
        # else:
        #   p.addUserDebugLine([self.pos[0], self.pos[1], 0.001], [self.pos2[0], self.pos2[1], 0], [1,0,0,], 8, self.bciRate)