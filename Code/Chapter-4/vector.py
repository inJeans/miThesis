import numpy as np


class Vector(object):
    def __init__(self, args):
        """ Create a vector, example: v = Vector([1,2]) """
        if len(args) == 0:
            self.x = 0.
            self.y = 0.
            self.values = np.array([0., 0.])
        elif len(args) == 1:
            self.x = args[0]
            self.values = np.array([self.x])
        elif len(args) == 2:
            self.x = args[0]
            self.y = args[1]
            self.values = np.array([self.x, self.y])
        elif len(args) == 3:
            self.x = args[0]
            self.y = args[1]
            self.z = args[2]
            self.values = np.array([self.x, self.y, self.z])
        elif len(args) == 4:
            self.x = args[0]
            self.y = args[1]
            self.z = args[2]
            self.w = args[3]
            self.values = np.array([self.x, self.y, self.z, self.w])
        else:
            self.values = np.array(args).astype(np.float64)

    def __floordiv__(self, scalar):
        return Vector(self.values // scalar)

    def __truediv__(self, scalar):
        return Vector(self.values / scalar)

    def __repr__(self):
        return "{0}".format(self.values)

    def __str__(self):
        return "{0}".format(self.values)
