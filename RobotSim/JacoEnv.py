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
    p.resetSimulation()
    p.configureDebugVisualizer(p.COV_ENABLE_GUI,0)
    p.setGravity(0,0,-10)

    self.useOrientation         = 1
    self.useSimulation          = 1
    self.useRealTimeSimulation  = 1

    self.mode           = mode
    self.angle          = angle
    self.dl             = dl
    self.newPosInput    = 0
    self.ikSolver       = 0
    self.dist           = 0.0075
    self.distf          = .0338
    self.debuglen       = 0.3
    self.fu             = 1.35
    self.fl             = 0.
    self.center         = np.array([-0.35, 0.3, 0.25])
    self.bciRate        = 0.125
    self.key            = 0
    self.robotTargetRad = .05
    self.opMode         = 0

    p.loadURDF("plane.urdf",[0,0,-.65])
    p.loadURDF("table/table.urdf", basePosition=[-0.6,0.45,-0.65])
    self.jacoId = p.loadURDF("URDFs/jaco/j2n6s300.urdf", [0,0,0],  useFixedBase=True)

    p.resetBasePositionAndOrientation(self.jacoId,[0,0,0],[0,0,0,1])

    # Camera settings
    if self.mode == 0 or self.mode == 2: 
      if self.angle == 0:
        p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw=0, cameraPitch=-30, cameraTargetPosition=[-0.35,0.3,0.1])
      else:
        p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw=-45, cameraPitch=-10, cameraTargetPosition=[-0.35,0.3,0.1])
    elif self.mode == 3  or self.mode == 4 or self.mode == 6 or self.mode == 8:
      p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw= 30, cameraPitch=-45, cameraTargetPosition=[-0.35,0.3,0.1])
    elif self.mode == 5:
      p.resetDebugVisualizerCamera(cameraDistance=0.4, cameraYaw= 30, cameraPitch=-45, cameraTargetPosition=[-0.35,0.3,0.1])
    elif self.mode == 9:
      p.resetDebugVisualizerCamera(cameraDistance=5., cameraYaw= 0, cameraPitch= 80, cameraTargetPosition=[-0.35,0.3,0.3])
    elif self.mode == 10:
      p.resetDebugVisualizerCamera(cameraDistance=0.7, cameraYaw= 0, cameraPitch=-10, cameraTargetPosition=[-0.35,0.3,0.1])
    elif self.mode == 12:
      p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw= 180, cameraPitch=-30, cameraTargetPosition=[-0.35,0.3,0.1])
    elif self.mode == 14:
      p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw= 180 + 30, cameraPitch=-45, cameraTargetPosition=[-0.35,0.3,0.1])  
    else: 
      p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw=0, cameraPitch=-30, cameraTargetPosition=[-0.35,0.3,0.1])

    self.jacoEndEffectorIndex = 8
    self.jacoArmJoints        = [2, 3, 4, 5, 6, 7]
    self.jacofingerJoints     = [9, 11, 13]
    self.jacoJoints           = [2, 3, 4, 5, 6, 7, 9, 11, 13]
    self.jd                   = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    self.JP                   = [0,0,0,0,0,0]

    self.RP = [0,2.4,-1.4,1.0*math.pi, 1.8*math.pi, 0*math.pi, 1.75*math.pi, 0.5*math.pi,0]

    self.cube1Id = p.loadURDF("cube_small.urdf",[-0.2, 0, -0.2] + self.center, [0,0,0,1],  globalScaling=3)
    # self.cube1Id = p.loadURDF("domino/domino.urdf",[0.2, 1.0, -2] + self.center, [ 0.4996018, 0.4999998, 0.4999998, 0.5003982 ],   globalScaling=3)

    c       = [0, 1, 0]
    pos     = np.array([0,0, -2])
    pos[0]  = self.center[0] + pos[0]
    pos[1]  = self.center[1] + pos[1]
    pos[2]  = self.center[2] + pos[2]
    lw      = 6
    d       = .05

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
    self.m1 = p.addUserDebugLine(self.c1, self.c2, c, lw, 0)

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

    letters = ["A", "F", "B", "E", "C", "D"]
    offset  = [[0,0,0], [0, -0.05, -0.05],  [0, -0.05, -0.05], [-.1, -0.05, -0.1], [0,0,-0.1], [0.1,-0.2,0.02], [0.05,-0.15,-0.1]]

    letter  = letters[targetID-1]
    pos     = pos + offset[targetID]
    self.letter = p.addUserDebugText(letter, (pos),  c, 12)
    self.LetterMode = 1

  def set_letterColor(self, dpos, c, targetID):

    pos = dpos.copy()
    pos[0] = self.center[0] + pos[0]
    pos[1] = self.center[1] + pos[1]
    pos[2] = self.center[2] + pos[2]

    letters = ["A", "F", "B", "E", "C", "D"]
    offset  = [[0,0,0], [0, -0.05, -0.0],  [0, -0.05, -0.05], [-.1, -0.05, -0.1], [0,0,-0.1], [0.1,-0.2,0.02], [0.05,-0.15,-0.1]]
    letter  = letters[targetID-1]
    pos     = pos + offset[targetID]
    p.addUserDebugText(letter, (pos),  c, 12,replaceItemUniqueId=self.letter)

    self.JP = [0,0,0,0,0,0]

  def set_block_pos(self, pos, target):
    pos[0] = self.center[0] + pos[0]
    pos[1] = self.center[1] + pos[1]

    d = self.robotTargetRad

    p.resetBasePositionAndOrientation(self.cube1Id, [pos[0], pos[1], 0], [ 0.4996018, 0.4999998, 0.4999998, 0.5003982 ])
    c1 = [pos[0] - d, pos[1] - d, -0.02]
    c2 = [pos[0] - d, pos[1] + d, -0.02]
    c3 = [pos[0] + d, pos[1] + d, -0.02]
    c4 = [pos[0] + d, pos[1] - d, -0.02]

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

    c1 = [pos[0] - d, pos[1] - d, 0.0]
    c2 = [pos[0] - d, pos[1] + d, 0.0]
    c3 = [pos[0] + d, pos[1] + d, 0.0]
    c4 = [pos[0] + d, pos[1] - d, 0.0]

    self.l1 = p.addUserDebugLine(c1, c2, col, 6, 0,replaceItemUniqueId=self.l1)
    self.l2 = p.addUserDebugLine(c2, c3, col, 6, 0, replaceItemUniqueId=self.l2)
    self.l3 = p.addUserDebugLine(c3, c4, col, 6, 0, replaceItemUniqueId=self.l3)
    self.l4 = p.addUserDebugLine(c4, c1, col, 6, 0,replaceItemUniqueId=self.l4)

  def drawLine(self, p1, p2):
    p.addUserDebugLine(p1, p2, [0,1,0],200,0)

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



  def draw2DAxes(self):

    c3 = [0, 0.2, -0.26] + self.center
    c4 = [0, -0.2, -0.26] + self.center
    c1 = [0.2, 0.0, -0.26] + self.center
    c2 = [-0.2, 0.0, -0.26] + self.center

    lw = 2

    p.addUserDebugLine(c1, c2, [0,0,0], lw, 0)
    p.addUserDebugLine(c3, c4, [0,0,0], lw, 0)

  def set_ringTarget(self, pos, c):

      pos[0] = self.center[0] + pos[0]
      pos[1] = self.center[1] + pos[1]
      pos[2] = self.center[2] + pos[2]
      lw = 6

      d = self.robotTargetRad
      dz = .03

      self.c1 = [pos[0] - d, pos[1]-d, pos[2]-d - dz]
      self.c2 = [pos[0] + d, pos[1]-d, pos[2]-d - dz]
      self.c3 = [pos[0] + d, pos[1]-d, pos[2]+d - dz]
      self.c4 = [pos[0] - d, pos[1]-d, pos[2]+d - dz]
      self.c5 = [pos[0] - d, pos[1]+d, pos[2]-d - dz]
      self.c6 = [pos[0] + d, pos[1]+d, pos[2]-d - dz]
      self.c7 = [pos[0] + d, pos[1]+d, pos[2]+d - dz]
      self.c8 = [pos[0] - d, pos[1]+d, pos[2]+d - dz]

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

  def set_ringColor(self, pos, c, lw):
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
  
  def set_cubeTarget(self, pos, c):

    pos[0] = self.center[0] + pos[0]
    pos[1] = self.center[1] + pos[1]
    pos[2] = self.center[2] + pos[2]
    lw = 6
    d = self.robotTargetRad
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


  def drawBetaLine(self, betaScalar):
    print(betaScalar)
    p1 = [-0.8, 0.1, 0]
    p2 = [-0.8, 0.1, betaScalar*0.5+ 0.1]

    c = [0, 0, 1]
    lw = 10

    self.l13 = p.addUserDebugLine(p1, p2, c, lw, 0, replaceItemUniqueId=self.l13)

  def set_triangleTarget(self, pos, c, target):

    if target == 9:
      xd = 1
    else:
      xd = -1

    pos[0] = self.center[0] + pos[0]
    pos[1] = self.center[1] + pos[1]
    pos[2] = self.center[2] + pos[2]
    lw = 6
    d = self.robotTargetRad

    self.c1 = [pos[0], pos[1], pos[2]-d]
    self.c2 = [pos[0] - 2*xd*d, pos[1] - d, pos[2]-d]
    self.c3 = [pos[0] - 2*xd*d, pos[1] + d, pos[2]-d]

    p.addUserDebugLine(self.c1, self.c2, c, lw, 0, replaceItemUniqueId=self.l1)
    p.addUserDebugLine(self.c2, self.c3, c, lw, 0, replaceItemUniqueId=self.l2)
    p.addUserDebugLine(self.c3, self.c1, c, lw, 0, replaceItemUniqueId=self.l3)

  def set_triangleColor(self, pos, c, lw):

    p.addUserDebugLine(self.c1, self.c2, c, lw, 0, replaceItemUniqueId=self.l1)
    p.addUserDebugLine(self.c2, self.c3, c, lw, 0, replaceItemUniqueId=self.l2)
    p.addUserDebugLine(self.c3, self.c1, c, lw, 0, replaceItemUniqueId=self.l3)

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


  def set_robotRotation(self, dim, rot_theta):
    Rrm = R.from_quat(self.orn)     # convert current orn to rot matrix
    Rnew =  Rrm.as_matrix() 

    if dim == 0:    #x
      Rx = np.array([[1., 0., 0.],[0., np.cos(rot_theta), -np.sin(rot_theta)], [0., np.sin(rot_theta), np.cos(rot_theta)]])
      Rnew = Rrm.as_matrix() @ Rx
      print("Update X")
    if dim == 1:  #y
      Ry = np.array([[np.cos(rot_theta), 0., np.sin(rot_theta)], [0., 1., 0.], [-np.sin(rot_theta), 0., np.cos(rot_theta)]])
      Rnew = Rrm.as_matrix() @ Ry
    if dim == 2:  #z
      Rz = np.array([[np.cos(rot_theta), -np.sin(rot_theta), 0.], [np.sin(rot_theta), np.cos(rot_theta), 0.], [0., 0., 1.]])
      Rnew = Rrm.as_matrix() @ Rz

    Rn = R.from_matrix(Rnew)
    self.orn = Rn.as_quat()
    self.newPosInput = 1

  def set_robotOrn(self, orn):
    # ornNew    = [math.pi, 0, orn]
    self.orn  = p.getQuaternionFromEuler(orn)
    self.newPosInput = 1
    

  def set_robotPos(self, rp, key):
    self.key    = key;
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

    self.newPosInput = 1

  def drawMode(self, mode):

    if self.mode == 12:
      c1 = [-0.45, 0.6, 0.];
      c2 = [-0.15, 0.6, 0.];
    else:
      c1 = [-0.45, 0, 0.];
      c2 = [-0.15, 0, 0.];
    if mode == 1:
      self.m1 = p.addUserDebugLine(c1,c2, [1,0,0], 20, replaceItemUniqueId=self.m1)
    else:
      self.m1 = p.addUserDebugLine(c1,c2, [0,0,1], 20, replaceItemUniqueId=self.m1)

  def setGoCue(self, val):
    p.removeUserDebugItem(self.t1)
    if val == 1:
      self.t1 = p.addUserDebugText('o o o', [-1.2,0.0,0.3],  [0,1,0], 16)
    elif val == 0:
      self.t1 = p.addUserDebugText('o o o', [-1.2,0.0,0.3],  [1,0,0], 16)

  def inverseKin(self):
    if (self.newPosInput == 1):
      self.jointPoses = p.calculateInverseKinematics(self.jacoId,
                                                self.jacoEndEffectorIndex,
                                                self.pos,
                                                self.orn,
                                                jointDamping=self.jd,
                                                solver=self.ikSolver,
                                                maxNumIterations=1000,
                                                residualThreshold=.0001,
                                                lowerLimits = [-3.14,-6.28,-6.28,-3.14,-3.14,-3.14,-3.14, -3.14,-3.14],
                                                upperLimits = [3.14,-1.6,6.28,3.14,3.14,3.14,3.14,3.14,3.14],
                                                jointRanges = [6,6,6,6,6,6,6,6,6],
                                                restPoses = self.RP)

      self.JP = list(self.jointPoses)
      # print(self.JP)

  def removeDebug(self):
    p.removeAllUserDebugItems()
    pos = np.array([0,0, -2])
    lw = 6
    d = .05
    c = [0,0,1]

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

    self.m1 = p.addUserDebugLine(self.c1, self.c2, c, lw, 0)

    d1 = [0,0,0]
    d2 = [0,0,0]
    p.addUserDebugLine(d1, d2, [0,0,0], 4, 0,replaceItemUniqueId=self.l13)

  def reset(self):
    # p.resetSimulation()
    # p.resetBasePositionAndOrientation(self.cube1Id, [-1., -1., -1.], [0,0,0,1])
    
    p.removeAllUserDebugItems()
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
      p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw= 15, cameraPitch=-30, cameraTargetPosition=[-0.35,0.3,0.1])
    elif self.mode == 5:
      p.resetDebugVisualizerCamera(cameraDistance=0.25, cameraYaw= 0, cameraPitch=0, cameraTargetPosition=[-0.35,0.3,0.2])  
    elif self.mode == 7:
      p.resetDebugVisualizerCamera(cameraDistance=0.5, cameraYaw= 25, cameraPitch=-20, cameraTargetPosition=[-0.35,0.3,0.2])    
    elif self.mode == 9:
      p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw= 30, cameraPitch=-45, cameraTargetPosition=[-0.35,0.3,0.1])
    elif self.mode == 11:
      p.resetDebugVisualizerCamera(cameraDistance=0.7, cameraYaw= 15, cameraPitch=-45, cameraTargetPosition=[-0.35,0.3,0.1])
    elif self.mode == 12:
      p.resetDebugVisualizerCamera(cameraDistance=0.55, cameraYaw= 180, cameraPitch=-30, cameraTargetPosition=[-0.35,0.3,0.1])
    elif self.mode == 14:
      p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw= 180 - 20, cameraPitch=-30, cameraTargetPosition=[-0.3,0.25,0.1])  
    elif self.mode == 15:
      p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw= 15, cameraPitch=-30, cameraTargetPosition=[-0.35,0.3,0.1])
    elif self.mode == 16: # RobotBeta
      p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw= 0, cameraPitch=0, cameraTargetPosition=[-0.35,0.3,0.3])
    else: 
      p.resetDebugVisualizerCamera(cameraDistance=0.6, cameraYaw=0, cameraPitch=-0, cameraTargetPosition=[-0.35,0.3,0.1])


    self.t1 = p.addUserDebugText('o o o', [-1.2,0.2,-1.4],  [0,1,0], 1)
    pos = np.array([0,0, -2])
    lw = 6
    d = .05
    c = [0,0,1]

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

    self.m1 = p.addUserDebugLine(self.c1, self.c2, c, lw, 0)

    d1 = [0,0,0]
    d2 = [0,0,0]
    p.addUserDebugLine(d1, d2, [0,0,0], 4, 0,replaceItemUniqueId=self.l13)

    if self.mode == 10 or self.mode == 12:
      self.draw2DAxes()
    elif self.mode < 5 or self.mode == 14:
        self.drawAxes()
    if self.mode == 9:
      self.drawAxes()
    rp = [0,math.pi/4,math.pi,1.0*math.pi, 1.8*math.pi, 0*math.pi, 1.75*math.pi, 0.5*math.pi]

    if self.mode == 0:
      self.pos =list([-0.35, 0.3, 0.2])
    elif self.mode ==3 or self.mode == 9:
      self.pos =list([-0.35, 0.3, 0.25])
    else:
      self.pos =list([-0.35, 0.3, 0.2])

    self.orn = p.getQuaternionFromEuler([math.pi,0,math.pi])
    # self.pos   = list([-0.5, 0.0, 0.25])
    # self.orn  = p.getQuaternionFromEuler([0, -math.pi/2, 0])

    if self.mode == 7:
      self.fing = 0.0
    elif self.mode == 9:
      self.fing = 1.35
    else:
      self.fing = 0.675
      
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

    if self.mode == 10 or self.mode == 12:
      p.resetBasePositionAndOrientation(self.cube1Id, [-0.25, 0, -0.2] + self.center, [0,0,0,1])
    else:
      p.resetBasePositionAndOrientation(self.cube1Id, [-0.25, 0, -2] + self.center, [0,0,0,1])

  def step(self):

    if self.newPosInput: 
      self.inverseKin()

    if self.opMode == 0:
      j = 0
      for i in self.jacoArmJoints:
        p.resetJointState(self.jacoId, i, self.jointPoses[j])
        j = j+1 
    elif self.opMode == 1:
      j = 1
      for i in [3,4,5,6,7]:
        p.resetJointState(self.jacoId, i, self.jointPoses[j])
        j = j+1 

    for i in  [9, 11, 13]:
      p.setJointMotorControl2(self.jacoId, i, p.POSITION_CONTROL, self.fing)

    self.newPosInput = 0
    
    # Debug lines for visualization

    p1 = [self.pos[0] - .02, self.pos[1], 0.002]
    p2 = [self.pos[0] + .02, self.pos[1], 0.002]
    p3 = [self.pos[0], self.pos[1] + .02, 0.002]
    p4 = [self.pos[0], self.pos[1] - .02, 0.002]

    if self.mode == 7 or self.mode == 10 or self.mode == 12:
      p.addUserDebugLine(p1, p2, [0,1,0], 6, self.bciRate)
      p.addUserDebugLine(p3, p4, [0,1,0], 6, self.bciRate) 

    p1 = [self.pos[0] - .02, self.pos[1], self.pos[2]]
    p2 = [self.pos[0] + .02, self.pos[1], self.pos[2]]
    p3 = [self.pos[0], self.pos[1] + .02, self.pos[2]]
    p4 = [self.pos[0], self.pos[1] - .02, self.pos[2]]

    if self.mode == 11 or self.mode == 15:
      p.addUserDebugLine(p1, p2, [0,0,1], 6, self.bciRate)
      p.addUserDebugLine(p3, p4, [0,0,1], 6, self.bciRate) 

    if self.dl:
      if self.mode == 0:
        p.addUserDebugLine([self.pos[0], self.pos[1], 0.001], [self.pos2[0], self.pos2[1], 0], [1,0,0,], 8, self.bciRate)
      elif self.mode == 2:
        p.addUserDebugLine([self.pos[0], self.pos[1], 0.001], [self.pos2[0], self.pos2[1], 0], [1,0,0,], 8, self.bciRate)
      elif self.mode ==1:
        p.addUserDebugLine([self.pos[0], self.pos[1], self.pos[2] + 0.05], [self.pos2[0], self.pos2[1], self.pos2[2] + .05], [1,0,0,], 8, self.bciRate)
      elif self.mode == 15:
        if self.key == 100:
          p1 = [self.pos[0] - .15, self.pos[1], self.pos2[2]]
          p2 = [self.pos[0] + .15, self.pos[1], self.pos2[2]]
          p3 = [self.pos[0], self.pos[1] + .15, self.pos2[2]]
          p4 = [self.pos[0], self.pos[1] - .15, self.pos2[2]]

          p.addUserDebugLine(p1, p2, [0,1,1], 8, self.bciRate)
          p.addUserDebugLine(p3, p4, [0,1,1], 8, self.bciRate)
        elif self.key == 101:
          p1 = [self.pos[0], self.pos[1], self.pos2[2]]
          p2 = [self.pos[0], self.pos[1], self.pos2[2] - 0.15]
          p3 = [self.pos[0], self.pos[1] + .15, self.pos2[2] - 0.15]
          p.addUserDebugLine(p1, p2, [0,1,1], 8, self.bciRate)
          p.addUserDebugLine(p2, p3, [0,1,1], 8, self.bciRate)
          print("Here")
        else:
          p.addUserDebugLine([self.pos[0], self.pos[1], self.pos[2]], [self.pos2[0], self.pos2[1], self.pos2[2]], [1,0,0,], 8, self.bciRate)  
      elif self.mode == 3 or self.mode == 4 or self.mode > 5:
        if self.key == 100:
          p1 = [self.pos[0] - .15, self.pos[1], self.pos2[2]]
          p2 = [self.pos[0] + .15, self.pos[1], self.pos2[2]]
          p3 = [self.pos[0], self.pos[1] + .15, self.pos2[2]]
          p4 = [self.pos[0], self.pos[1] - .15, self.pos2[2]]

          p.addUserDebugLine(p1, p2, [0,1,1], 8, self.bciRate)
          p.addUserDebugLine(p3, p4, [0,1,1], 8, self.bciRate)

        elif self.key == 101:
          p1 = [self.pos[0], self.pos[1], self.pos[2] + 0.1]
          p2 = [self.pos[0], self.pos[1], self.pos[2]-0.05]
          p3 = [self.pos[0] + 0.075, self.pos[1], self.pos[2] - 0.1]
          p4 = [self.pos[0] + 0.15, self.pos[1], self.pos[2] - 0.05]
          p5 = [self.pos[0] + 0.15, self.pos[1], self.pos[2] + 0.05]
          p.addUserDebugLine(p1, p2, [1,1,0], 8, self.bciRate)
          p.addUserDebugLine(p2, p3, [1,1,0], 8, self.bciRate)
          p.addUserDebugLine(p3, p4, [1,1,0], 8, self.bciRate)
          p.addUserDebugLine(p4, p5, [1,1,0], 8 , self.bciRate)

        elif self.key == 102:
          p1 = [self.pos[0], self.pos[1], self.pos[2] + 0.1]
          p2 = [self.pos[0], self.pos[1], self.pos[2]-0.05]
          p3 = [self.pos[0] - 0.075, self.pos[1], self.pos[2] - 0.1]
          p4 = [self.pos[0] - 0.15, self.pos[1], self.pos[2] - 0.05]
          p5 = [self.pos[0] - 0.15, self.pos[1], self.pos[2] + 0.05]
          p.addUserDebugLine(p1, p2, [1,1,0], 8, self.bciRate)
          p.addUserDebugLine(p2, p3, [1,1,0], 8, self.bciRate)
          p.addUserDebugLine(p3, p4, [1,1,0], 8, self.bciRate)
          p.addUserDebugLine(p4, p5, [1,1,0], 8 , self.bciRate)
    
        else:
          p.addUserDebugLine([self.pos[0], self.pos[1], self.pos[2] + 0.05], [self.pos2[0], self.pos2[1], self.pos2[2] + .05], [1,0,0,], 8, self.bciRate)