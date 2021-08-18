from openJacoEnv import JacoEnv
import numpy as np

class DiscreteActionsRobot():
    """Interface for discrete action classification."""
    # TODO: Optimize for srpites. Currently drawing a new arrow when changing size/rotation.
    def __init__(self, *args, **kwargs):  
        
        for _k, _v in kwargs.items():
            if hasattr(self, _k):
                setattr(self, _k, _v)

    def open(self):
        self.robotenv = JacoEnv(self.mode, self.angle, self.debugLines)
        self.robotenv.center_pos = np.array([-0.35, -0.3])
        self.pos = self.robot_to_bci_transform([self.robotenv.pos[0], self.robotenv.pos[1]])
        # self.robotenv.drawAxes()

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
        # print("TARGET:", target)
        self.target = target
        self.robotenv.set_block_pos(pos, target)

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

    def draw_bound(self, cpos):
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
        self.robotenv.draw_bound(pos)

    def update_bound(self, cpos, col):
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
        self.robotenv.update_bound(pos, col)



    def create_target3D(self, cpos, st):

        # self.targetPos = self.bci_to_robot_transform(pos)
        pos = cpos.copy()
        if st == 1:
            color = [0,1,0];
            self.robotenv.set_cubeTarget(pos, color)
            self.newTarget = 1
        elif st == 0:
            color = [1,0,0];
            if self.newTarget ==1:
                self.robotenv.set_cubeColor(pos, color, 16)
                self.newTarget = 0
        elif st == 2:
            color = [0,0,1];
            self.robotenv.set_cubeColor(pos, color, 22)
        else:
            color = [0,1,0];
            self.robotenv.set_cubeColor(pos, color, 22)


        

    def render(self):
        self.robotenv.step()

    def updateRefresh(self, rate):
        self.robotenv.bciRate = rate

    def updateMode(self, mode):
        self.mode = mode
        self.robotenv.mode = self.mode

    def updateDebugLines(self, dl):
        self.debugLines = dl
        self.robotenv.dl = self.debugLines

    def update_joystick(self, key):
        if self.mode == 0:
            if key == 1:
                self.key = 6
            elif key == 4:
                self.key = 8
            elif key == 3:
                self.key = 4
            elif key == 2:
                self.key = 2
            else:
                self.key = 0
        elif self.mode == 1:
            if key == 1:
                self.key = 16
            elif key == 4:
                self.key = 18
            elif key == 3:
                self.key = 14
            elif key == 2:
                self.key = 12
            else:
                self.key = 0
        elif self.mode == 2:
            if key == 1:
                self.key = 26
            elif key == 2:
                self.key = 28
            elif key == 3:
                self.key = 24
            elif key == 4:
                self.key = 22
            else:
                self.key = 0
        if self.mode == 3:
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
        if self.mode == 4:
            if key == 1:
                self.key = 46
            elif key == 4:
                self.key = 48
            elif key == 3:
                self.key = 44
            elif key == 2:
                self.key = 42
            elif key == 5:
                self.key = 47
            elif key == 6:
                self.key = 41
            elif key == 7:
                self.key = 100
            else:
                self.key = 0
                
        # print(self.key)
        self.robotenv.updateCommand(self.key)
        # self.pos =self.robot_to_bci_transform([self.robotenv.pos[0], self.robotenv.pos[1]])
        # print(self.key)

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



