import func_calc as f
from copy import deepcopy


def get_func_point(func, point):
    value = func[0]
    for i in range(len(func)-1, 0, -1):
        value += func[i] * point**i
    return value


def method_div_segments(func, interval, epsilon, extr):
    int_ = deepcopy(interval)
    x_m = (int_[0] + int_[1]) / 2
    L = int_[1] - int_[0]

    while abs(L) >= epsilon:
        x_1 = int_[0] + L / 4
        x_2 = int_[1] - L / 4
        y_1 = func(x_1)
        y_2 = func(x_2)
        y_m = func(x_m)

        if extr == 'min':
            if y_1 < y_m:
                int_[1] = x_m
                x_m = x_1
            else:
                if y_m > y_2:
                    int_[0] = x_m
                    x_m = x_2
                else:
                    int_[1] = x_2
                    int_[0] = x_1

        else:
            if y_1 < y_m:
                int_[0] = x_m
                x_m = x_2
            else:
                if y_m > y_2:
                    int_[1] = x_m
                    x_m = x_1
                else:
                    int_[1] = x_1
                    int_[0] = x_2

        L = int_[1] - int_[0]

    return x_m, func(x_m)


def method_fibonacci(func, interval, n_op, extr):
    int_ = deepcopy(interval)
    L = int_[1] - int_[0]
    epsilon = round(L / get_fibonacci(n_op+1), 3) - 0.001
    L_2 = get_fibonacci(n_op-1) / get_fibonacci(n_op) * L + (-1)**n_op * epsilon / get_fibonacci(n_op)
    x_1 = int_[0] + L_2  # x_2
    x_2 = int_[0] + int_[1] - x_1  # x_4

    for i in range(1, n_op):
        y_1 = func(x_1)
        y_2 = func(x_2)

        if extr == 'min':
            if x_2 < x_1 and y_2 <= y_1:
                int_[1] = x_1
                x_1 = int_[0] + int_[1] - x_2
                x_m = x_2
            elif x_2 < x_1 and y_2 > y_1:
                int_[0] = x_2
                x_2 = int_[1] + int_[0] - x_1
                x_m = x_1
            elif x_2 > x_1 and y_2 <= y_1:
                int_[0] = x_1
                x_1 = int_[1] + int_[0] - x_2
                x_m = x_2
            elif x_2 > x_1 and y_2 > y_1:
                int_[1] = x_2
                x_2 = int_[0] + int_[1] - x_1
                x_m = x_1
        else:
            if x_2 < x_1 and y_2 <= y_1:
                int_[0] = x_2
                x_2 = int_[0] + int_[1] - x_1
                x_m = x_1
            elif x_2 < x_1 and y_2 > y_1:
                int_[1] = x_1
                x_1 = int_[1] + int_[0] - x_2
                x_m = x_2
            elif x_2 > x_1 and y_2 <= y_1:
                int_[1] = x_2
                x_2 = int_[1] + int_[0] - x_1
                x_m = x_1
            elif x_2 > x_1 and y_2 > y_1:
                int_[0] = x_1
                x_1 = int_[0] + int_[1] - x_2
                x_m = x_2

        L = int_[1] - int_[0]
        if abs(L) < epsilon:
            break

    return x_m, func(x_m)


def method_golden_section(func_str, interval, extr, epsilon):
    func_str_ = func_str
    if extr == 'max':
        func_str_ = f'-({func_str_})'
    function = f.func(func_str_)

    a, b = interval
    x1 = a + 0.382 * (b - a)
    x2 = b - 0.382 * (b - a)

    while True:
        y1 = function(x1)
        y2 = function(x2)

        if y1 < y2:
            b = x2
            if b - a > epsilon:
                x2 = deepcopy(x1)
                x1 = a + 0.382 * (b - a)
            else:
                return (a + b) / 2

        elif y1 > y2:
            a = x1
            if b - a > epsilon:
                x1 = deepcopy(x2)
                x2 = b - 0.382 * (b - a)
            else:
                return (a + b) / 2


def get_fibonacci(n):
    F = [1, 1]
    for i in range(1, n):
        F.append(F[i] + F[i-1])
    return F[n]


if __name__ == '__main__':
    func_str = '(4+x*(-39.0016907658719))**3+(3+x*(-14.999557151895715))**3-3*(4+x*(-39.0016907658719))*(3+x*(-14.999557151895715))'
    func = f.func(func_str)
    interval = [0, 2]
    extr = 'max'

    epsilon = 0.01
    # X, f_X = method_div_segments(func, interval, epsilon, extr)
    # print('Метод деления отрезков пополам: \n x = {}, f(x) = {}'.format(X, f_X))
#
    # n_op = 20
    # X, f_X = method_fibonacci(func, interval, n_op, extr)
    # print('Метод Фибоначчи: \n x = {}, f(x) = {}'.format(X, f_X))

    X = method_golden_section(func_str, interval, extr, epsilon)
    print('Метод золотого сечения: \n x = {}'.format(X))
