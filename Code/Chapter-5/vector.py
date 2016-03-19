import numpy as np


class Vector(np.ndarray):
    def __new__(cls, input_array):
        """ Create a vector, example: v = Vector([1,2]) """
        obj = np.asarray(input_array).view(cls)
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        # initialise class variables

    # def __init__(self, vec):
        
    #     self._values = self.__new__(Vector, vec)

    def __array_wrap__(self, out_arr, context=None):
        return np.ndarray.__array_wrap__(self, out_arr, context)

    def norm(self):
        return np.linalg.norm(self)

    def unit(self):
        return Vector(self / self.norm())

    @property
    def x(self):
        return self[0]

    @x.setter
    def x(self, x):
        self[0] = x

    @property
    def y(self):
        return self[1]

    @y.setter
    def y(self, y):
        self[1] = y

    @property
    def z(self):
        return self[2]

    @z.setter
    def z(self, z):
        self[2] = z

    @property
    def w(self):
        return self[3]

    @w.setter
    def w(self, w):
        self[3] = w

    # def __floordiv__(self, scalar):
    #     return Vector(self.data // scalar)

    # def __truediv__(self, scalar):
    #     return Vector(self.data / scalar)

    # def __repr__(self):
    #     return "{0}".format(self.data)

    # def __str__(self):
    #     return "{0}".format(self.data)
