from JacoEnv import JacoEnv
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

    def create_target(self, pos):
        # self.targetPos = self.bci_to_robot_transform(pos)
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

    def render(self):
        self.robotenv.step()

    def updateRefresh(self, rate):
        self.robotenv.bciRate = rate

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
                
        # print(self.key)
        self.robotenv.updateCommand(self.key)
        # self.pos =self.robot_to_bci_transform([self.robotenv.pos[0], self.robotenv.pos[1]])
        # print(self.key)

    def bci_to_robot_transform(self, pos):
        robot_pos = pos / 1200.0 + self.robotenv.center_pos
        return robot_pos

    def robot_to_bci_transform(self, pos):
        bci_pos = (pos - self.robotenv.center_pos)*1200.0
        return bci_pos
        
    def reset(self, arrow=True, target=False, all_objects=False):
        self.robotenv.reset()



