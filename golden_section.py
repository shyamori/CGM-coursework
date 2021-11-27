from func_calc import func
from copy import deepcopy


def method_golden_section(func_str, interval, extr, epsilon):
    func_str_ = func_str
    if extr == 'max':
        func_str_ = f'-({func_str_})'
    function = func(func_str_)

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
