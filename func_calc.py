from math import sin, cos, tan, log, exp, pi, sqrt

def func(fun):
    def und_func(x):
        return eval(fun)

    return und_func

