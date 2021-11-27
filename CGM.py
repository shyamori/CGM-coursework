from func_calc import func
import numpy as np
from math import sqrt
from copy import deepcopy
from sympy import *
from golden_section import method_golden_section


def get_gradient(function, x_0):
    """
    Вычисление приближенного значения градиента функции function в точке x_0 с шагом dx

    :param function: функция, для которой нужно вычислить градиент
    :param x_0: точка, в которой нужно вычислить градиент
    :return: значение градиента функции function в точке x_0
    """

    dx = 1e-12
    grad = []
    for i in range(len(x_0)):
        x_1 = deepcopy(x_0)
        x_1[i] += dx
        grad.append((function(x_1) - function(x_0)) / dx)

    return grad


def get_norm(x_1, x_2=None):
    """
    Вычисление евклидовой нормы разности векторов x_1 и x_2, либо евклидовой нормы вектора x_1

    :param x_1: координаты первого вектора
    :param x_2: координаты второго вектора (если имеется)
    :return: значение нормы
    """

    norm = 0
    if x_2 is not None:
        for i in range(len(x_1)):
            norm += (x_1[i] - x_2[i]) ** 2
    else:
        for x_1 in x_1:
            norm += x_1 ** 2
    norm = sqrt(norm)
    return norm


def get_extremum(function_str, extremum):
    """
    Нахождение точки экстремума функции f_str из необходимого и достаточного условия экстремума функции

    :param function_str: функция в виде строки
    :param extremum: критерий поиска ('min' / 'max')
    :return: координаты найденной точки экстремума
    """

    x = Symbol('x')
    # первая производная приравнивается к 0, решается уравнение
    dots_eq = solve(diff(function_str, x), x)
    # нахождение второй производной
    f_second_derivative = func(str(diff(diff(function_str, x), x)))
    for dot_eq in dots_eq:
        if type(dot_eq) != Add:
            # если в найденной точке dot_eq вторая производная > 0, то dot_eq является точкой минимума
            if extremum == 'min' and f_second_derivative(dot_eq) > 0:
                return dot_eq
            # если в найденной точке dot_eq вторая производная < 0, то dot_eq является точкой максимума
            if extremum == 'max' and f_second_derivative(dot_eq) < 0:
                return dot_eq


def conjugate_gradient_method(function_str, extremum, dot0, max_iter, quadratic, eps1=1e-4, eps2=1e-4,
                              golden_section=False,
                              gs_interval=None):
    """
    Реализация метода сопряженных градиентов

    :param function_str: целевая функция в виде строки
    :param extremum: критерий поиска ('min' / 'max')
    :param dot0: начальная точка
    :param eps1, eps2: значения погрешностей
    :param max_iter: максимальное число итераций
    :param quadratic: True - целевая функция квадратичная,
                      False - целевая функция неквадратичная
    :param golden_section: True - использование метода золотого сечения при поиске t,
                           False - поиск t с использованием необходимого и достаточного условия экстремума функции
    :param gs_interval: интервал поиска для метода золотого сечения
    :return: найденная точка экстремума, значение функции в этой точке,
                        массив mas_results со значениями приближений точки и функции на каждой итерации;
             либо сообщение в виде строки
    """

    if gs_interval is None:
        gs_interval = (0, 1)
    count = 0
    iter_num = 0
    mem_er = 0
    direct0 = []
    grad0 = []
    mas_results = []
    function = func(function_str)

    while True:
        mas_results.append([dot0, function(dot0)])

        # нахождение градиента grad1 в точке dot0
        grad1 = get_gradient(function, dot0)

        # проверка условий остановки
        if get_norm(grad1) < eps1:
            return dot0, function(dot0), mas_results
        if iter_num >= max_iter:
            return dot0, function(dot0), mas_results

        # вычисление направления direct1 на первом шаге
        if iter_num == 0:
            if extremum == 'min':
                direct1 = [-grad1[i] for i in range(len(dot0))]
            else:
                direct1 = grad1
        # вычисление направления direct1 на последующих шагах
        else:
            # если функция квадратичная, то бета ищется по формуле Флетчера-Ривса
            if quadratic is True:
                beta = get_norm(grad1) ** 2 / get_norm(grad0) ** 2
            # если функция неквадратичная, то бета ищется по формуле Полака-Рибьера
            else:
                if iter_num % len(dot0) != 0:
                    temp = [grad1[i] - grad0[i] for i in range(len(grad0))]
                    beta = np.dot(grad1, temp) / get_norm(grad0) ** 2
                else:
                    beta = 0
            # вычисление направления с найденным значением бета
            if extremum == 'min':
                direct1 = [-grad1[i] + beta * direct0[i] for i in range(len(dot0))]
            else:
                direct1 = [grad1[i] + beta * direct0[i] for i in range(len(dot0))]

        # нахождение величины шага t из условия f(dot0 + t * direct1) -> extr
        # запись целевой функции f(x) в виде f(dot0 + t * direct1)
        fi_str = deepcopy(function_str)
        for k in range(len(dot0)):
            fi_str = fi_str.replace(f'x[{k}]', f'({dot0[k]} + x*({direct1[k]}))')

        # решение задачи одномерной оптимизации
        if golden_section is False and mem_er == 0:
            try:
                t = get_extremum(fi_str, extremum)
            except MemoryError:
                mem_er = 1
                print('Использован метод золотого сечения')
        if golden_section is True or mem_er == 1:
            t = method_golden_section(fi_str, gs_interval, extremum, 0.001)

        if t is not None:
            if t < 0:
                return f'Не было найдено значение t{iter_num}, удовлетворяющее условиям.' \
                       f'\nBозможно, неверно выбрана начальная точка.'
        else:
            return f'Не было найдено значение t{iter_num}, удовлетворяющее условиям.' \
                   f'\nBозможно, неверно выбрана начальная точка.'

        # вычисление новой точки
        dot1 = [dot0[i] + t * direct1[i] for i in range(len(dot0))]
        print(f't{iter_num} = {round(t, 4)}\t\tx{iter_num + 1} = {dot1}\tf(x{iter_num + 1}) = {function(dot1)}')

        # проверка условия остановки
        if get_norm(dot0, dot1) < eps2 and abs(function(dot1) - function(dot0)) < eps2:
            count += 1
            if count == 2:
                return dot1, function(dot1), mas_results

        iter_num += 1
        dot0 = dot1
        direct0 = direct1
        grad0 = grad1


