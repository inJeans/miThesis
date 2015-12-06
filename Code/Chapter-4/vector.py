import numpy as np


class Vector(object):
    def __init__(self, args):
        """ Create a vector, example: v = Vector([1,2]) """
        self._values = np.zeros(len(args))
        for i, element in enumerate(args):
            self._values[i] = element

    def norm(self):
        return np.linalg.norm(self._values)

    def unit(self):
        return Vector(self._values / self.norm())

    @property
    def x(self):
        return self._values[0]

    @x.setter
    def x(self, x):
        self._values[0] = x

    @property
    def y(self):
        return self._values[1]

    @y.setter
    def y(self, y):
        self._values[1] = y

    @property
    def z(self):
        return self._values[2]

    @z.setter
    def z(self, z):
        self._values[2] = z

    @property
    def w(self):
        return self._values[3]

    @w.setter
    def w(self, w):
        self._values[3] = w

    def __floordiv__(self, scalar):
        return Vector(self._values // scalar)

    def __truediv__(self, scalar):
        return Vector(self._values / scalar)

    def __repr__(self):
        return "{0}".format(self._values)

    def __str__(self):
        return "{0}".format(self._values)
