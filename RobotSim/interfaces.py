# from jacoEnvDebug import JacoEnv

from JacoEnv import JacoEnv

import numpy as np

class DiscreteActionsRobot():
    def __init__(self, *args, **kwargs):  
        
        for _k, _v in kwargs.items():
            if hasattr(self, _k):
                setattr(self, _k, _v)

    def open(self):
        self.robotenv = JacoEnv(self.mode, self.angle, self.debugLines)
        self.robotenv.center_pos = np.array([-0.35, -0.3])
        self.pos = self.robot_to_bci_transform([self.robotenv.pos[0], self.robotenv.pos[1]])
        self.path =  [[None]*3 for i in range(4)]
    def create_target(self, cpos):
        # self.targetPos = self.bci_to_robot_transform(pos)
        pos = cpos.copy()
        if (pos[0] > 0) & (pos[1] == 0):
            target = 1
        elif (pos[0]  == 0) & (pos[1] > 0):
            target = 2
        elif (pos[0] < 0) & (pos[1] == 0):
            target = 3
        else:
            target = 4
        print("TARGET:", target)
        self.target = target
        self.robotenv.set_block_pos(pos, target)

    def create_target3D(self, cpos, st):
        # self.targetPos = self.bci_to_robot_transform(pos)
        pos = cpos.copy()
        if st == 1:
            color = [0,1,0];
            self.robotenv.set_cubeTarget(pos, color)
            self.newTarget = 1
        elif st == 0:
            color = [1,0,0];
            # if self.newTarget ==1:
            self.robotenv.set_cubeColor(pos, color, 16)
                # self.newTarget = 0
        elif st == 2:
            color = [0,0,1];
            self.robotenv.set_cubeColor(pos, color, 22)
        else:
            color = [0,1,0];
            self.robotenv.set_cubeColor(pos, color, 22)

    def create_targetRing(self, cpos, st):
        # self.targetPos = self.bci_to_robot_transform(pos)
        pos = cpos.copy()
        if st == 1:
            color = [0,1,0];
            self.robotenv.set_ringTarget(pos, color)
            self.newTarget = 1
        elif st == 0:
            color = [1,0,0];
            # if self.newTarget ==1:
            self.robotenv.set_ringColor(pos, color, 16)
                # self.newTarget = 0
        elif st == 2:
            color = [0,0,1];
            self.robotenv.set_ringColor(pos, color, 22)
        else:
            color = [0,1,0];
            self.robotenv.set_ringColor(pos, color, 22)

    def create_targetLetter(self, cpos, st):
        # self.targetPos = self.bci_to_robot_transform(pos)
        pos = cpos.copy()
        if st == 1:
            color = [0,0,1];
            self.robotenv.write_letter(pos, color, self.targetID)
            self.newTarget = 1
        elif st == 0:
            color = [1,0,0];
            # if self.newTarget ==1:
            self.robotenv.set_letterColor(pos, color, self.targetID)
                # self.newTarget = 0
        elif st == 2:
            color = [0,0,1];
            self.robotenv.set_letterColor(pos, color, self.targetID)
        else:
            color = [0,1,0];
            self.robotenv.set_letterColor(pos, color, self.targetID)
    
    def render(self):
        self.robotenv.step()

    def updateRefresh(self, rate):
        self.robotenv.bciRate = rate

    def updateMode(self, mode):
        self.mode = mode
        self.robotenv.mode = self.mode
        self.robotenv.reset()

    def setTargetRad(self, rad):
        print('setRad')
        self.robotenv.robotTargetRad = rad

    def updateDebugLines(self, dl):
        self.debugLines = dl
        self.robotenv.dl = self.debugLines

    def updateRobotPos(self, rp, key):
        if key == 1:
            self.key = 6
        elif key == 4:
            self.key = 8
        elif key == 3:
            self.key = 4
        elif key == 2:
            self.key = 2
        elif key == 5:
            self.key = 7
        elif key == 6:
            self.key = 1
        elif key == 7:
            self.key = 100
        else:
            self.key = 0
        self.robotenv.set_robotPos(rp, self.key)

    def bci_to_robot_transform(self, pos):
        robot_pos = pos / 1200.0 + self.robotenv.center_pos
        return robot_pos

    def robot_to_bci_transform(self, pos):
        bci_pos = (pos - self.robotenv.center_pos)*1200.0
        return bci_pos
        
    def reset(self, arrow=True, target=False, all_objects=False):
        self.robotenv.reset()

    def displayCue(self, cue, c):
        self.robotenv.displayCue(cue, c)

    def setFing(self, fing):
        self.robotenv.setFing(fing)

    def setMode(self, mode):
        self.robotenv.mode = mode

    def grabCube(self):
        self.robotenv.grabCube()
        self.render()

    def update_color(self, cpos, c):
        # self.targetPos = self.bci_to_robot_transform(pos)
        pos = cpos.copy()
        if (pos[0] > 0) & (pos[1] == 0):
            target = 1
        elif (pos[0]  == 0) & (pos[1] > 0):
            target = 2
        elif (pos[0] < 0) & (pos[1] == 0):
            target = 3
        else:
            target = 4
        # print("TARGET:", target)
        self.target = target
        self.robotenv.set_bound_color(pos, c)

    def setPath(self, ind):

        # horizontal 1 (left, forward, left forward)
        p1 = [[0.2, 0.2, 0], [0., 0.2, 0],[0,0,0], [-0.2, 0,0], [-0.2, -0.2, 0]]

        # vertical 1 (up, left, down, left)
        p2 = [[0.2,0,-0.1], [0.2, 0.0, 0.1], [0.0, 0.0, 0.1],[0., 0., -0.1], [-0.2,0,-0.1]]

        # diagonal 1
        p3 = [[0.2,0.1, 0],[0, 0.1, 0],[-0.1414, -.0414,0],[0, -.18,0],[0.2, -.18,0]]

        p = [p1,p1,p2,p2,p3,p3]
        self.path = p[ind]


    def setPath_ES(self, p, ind):
        print(self.path[ind])
        self.path[ind][0] = p[0]
        self.path[ind][1] = p[1]
        self.path[ind][2] = p[2]
        print(self.path)



    def drawPath(self):

        # p = [[0,0,0], [0.3, 0,0], [0.3, 0.2, 0],[0.3, 0.2, -0.2]] 
        p = self.path
        print(self.path)
        for i in range(len(p) - 1):
            self.robotenv.drawLine(p[i] + self.robotenv.center, p[i+1]+ self.robotenv.center)


    def setGoCue(self, c):
        self.robotenv.setGoCue(c)