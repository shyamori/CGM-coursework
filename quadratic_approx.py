import copy
from func_calc import func
import numpy as np


def method_powell(func_str, extr, dot1, delta, eps):
    func_str_ = func_str
    if extr == 'max':
        func_str_ = f'-({func_str_})'
    function = func(func_str_)
    dot = [dot1, 0, 0]

    while True:
        dot[1] = dot[0] + delta

        if function(dot[0]) > function(dot[1]):
            dot[2] = dot[0] + 2 * delta
        else:
            dot[2] = dot[0] - delta

        func_dot = [function(dot[0]), function(dot[1]), function(dot[2])]
        index_min = np.argmin(func_dot)
        dot_min = dot[index_min]

        # a1 = (func_dot[1] - func_dot[0]) / (dot[1] - dot[0])
        # a2 = 1 / (dot[2] - dot[1]) * ((func_dot[2] - func_dot[0]) / (dot[2] - dot[0]) - (func_dot[1] - func_dot[0])
        #                               / (dot[1] - dot[0]))
#
        # dot_line = ((dot[1] + dot[1]) / 2) - (a1 / (2 * a2))

        dot_line = 0.5 * (((dot[1]**2 - dot[2]**2) * func_dot[0] + (dot[2]**2 - dot[0]**2) * func_dot[1] +
                           (dot[0]**2 - dot[1]**2) * func_dot[2]) / ((dot[1] - dot[2]) * func_dot[0] +
                                                                     (dot[2] - dot[0]) * func_dot[1] +
                                                                     (dot[0] - dot[1]) * func_dot[2]))

        if abs((function(dot_line) - function(dot_min)) / function(dot_line)) < eps and \
                abs((dot_line - dot_min) / dot_line) < eps:
            return dot_line

        if dot[0] < dot_line < dot[2]:
            if function(dot_line) < function(dot_min):
                if dot[1] < dot_line:
                    dot = [dot[1], dot_line, dot[2]]
                else:
                    dot = [dot[0], dot_line, dot[1]]
            else:
                if dot[1] < dot_min:
                    dot = [dot[1], dot_min, dot_min + delta]
                else:
                    dot = [dot_min - delta, dot_min, dot[1]]
        else:
            dot[0] = dot_line


def method_powell_old(func, interval, extr, x1, dx, eps1, eps2):
    int = copy.deepcopy(interval)
    x = [x1, 0, 0]
    y = [func(x1), 0, 0]
    x[1] = x[0] + dx
    y[1] = func(x[1])

    if extr == 'min':
        if y[0] > y[1]:
            x[2] = x[0] + 2 * dx
        else:
            x[2] = x[0] - dx
        y[2] = func(x[2])

        while True:
            y_extr = y[0]
            i_min = 0
            for i in range(3):
                if y[i] < y_extr:
                    y_extr = y[i]
                    i_min = i
            x_extr = x[i_min]

            a1 = (y[1] - y[0]) / (x[1] - x[0])
            a2 = 1 / (x[2] - x[1]) * ((y[2] - y[0]) / (x[2] - x[0]) - (y[1] - y[0]) / (x[1] - x[0]))
            x_min_ = ((x[0] + x[1]) / 2) - (a1 / (2 * a2))
            y_min_ = func(x_min_)

            if abs((y_extr - y_min_) / y_min_) <= eps1 and abs((x_extr - x_min_) / x_min_) <= eps2:
                break

            if y_extr < y_min_:
                x1 = x_extr - dx
                x3 = x_extr + dx
                for i in range(3):
                    if x[i] < x_extr:
                        if x[i] > x1:
                            x1 = x[i]
                    if x[i] > x_extr:
                        if x[i] < x3:
                            x3 = x[i]
                x = [x1, x_extr, x3]
                y = [func(x1), func(x_extr), func(x3)]

            else:
                x1 = x_min_ - dx
                x3 = x_min_ + dx
                for i in range(3):
                    if x[i] < x_min_:
                        if x[i] > x1:
                            x1 = x[i]
                    if x[i] > x_min_:
                        if x[i] < x3:
                            x3 = x[i]
                x = [x1, x_min_, x3]
                y = [func(x1), func(x_min_), func(x3)]

    else:
        if y[0] < y[1]:
            x[2] = x[0] + 2 * dx
        else:
            x[2] = x[0] - dx
        y[2] = func(x[2])

        while True:
            y_extr = y[0]
            i_max = 0
            for i in range(3):
                if y[i] > y_extr:
                    y_extr = y[i]
                    i_max = i
            x_extr = x[i_max]

            a1 = (y[1] - y[0]) / (x[1] - x[0])
            a2 = 1 / (x[2] - x[1]) * ((y[2] - y[0]) / (x[2] - x[0]) - (y[1] - y[0]) / (x[1] - x[0]))
            x_max_ = (x[0] + x[1]) / 2 - a1 / 2 * a2
            y_max_ = func(x_max_)

            if abs((y_extr - y_max_) / y_max_) <= eps1 and abs((x_extr - x_max_) / x_max_) <= eps2:
                break

            if y_extr > y_max_:
                x1 = x_extr - dx
                x3 = x_extr + dx
                for i in range(3):
                    if x[i] < x_extr:
                        if x[i] > x1:
                            x1 = x[i]
                    if x[i] > x_extr:
                        if x[i] < x3:
                            x3 = x[i]
                x = [x1, x_extr, x3]
                y = [func(x1), func(x_extr), func(x3)]

            else:
                x1 = x_max_ - dx
                x3 = x_max_ + dx
                for i in range(3):
                    if x[i] < x_max_:
                        if x[i] > x1:
                            x1 = x[i]
                    if x[i] > x_max_:
                        if x[i] < x1:
                            x3 = x[i]
                x = [x1, x_max_, x3]
                y = [func(x1), func(x_max_), func(x3)]

    if x_extr > int[1]:
        x_extr = int[1]
    elif x_extr < int[0]:
        x_extr = int[0]

    return x_extr


if __name__ == '__main__':
    # func = func('5+3*(-2+x*18.999995997859287)-4*(-4+x*7.999998999252966)+(-2+x*18.999995997859287)**2*(-4+x*7.999998999252966)-(-4+x*7.999998999252966)**2')
    # interval = [0, 10]
    # extr = 'max'
    # x1 = 1
    # dx = 1
    # eps1 = 0.1
    # eps2 = 0.1
    # X, f_X = method_powell(func, interval, extr, x1, dx, eps1, eps2)
#
    # print('Метод Пауэлла: \n x = {}, f(x) = {}'.format(X, f_X))

    func_str = '2*x**2-12*x'
    extr = 'min'
    dot1 = 1
    dx = 0.001
    eps1 = 0.001
    print(method_powell(func_str, extr, dot1, dx, eps1))
