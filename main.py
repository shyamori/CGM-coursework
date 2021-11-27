from CGM import conjugate_gradient_method, get_norm
import math

f_str = 'sin(x0+x1)+(x0-x1)**2-1.5*x0+2.5*x1+1'
extremum = 'min'
dot0 = [-9, 8]
epsilon1 = 1e-4
epsilon2 = 1e-4
max_iter = 30
quadratic = True
golden_section = True
gs_interval = (0, 1)


i = 0
while True:
    if f'x{i}' in f_str:
        f_str = f_str.replace(f'x{i}', f'x[{i}]')
        i += 1
    else:
        break

output = conjugate_gradient_method(function_str=f_str, extremum=extremum, dot0=dot0, eps1=epsilon1, eps2=epsilon2,
                                   max_iter=max_iter, quadratic=quadratic, golden_section=golden_section,
                                   gs_interval=gs_interval)

dot_orig = [-0.54719, -1.54719]
func_orig = -1.91322

if type(output) is str:
    print(output)
else:
    X, fX, approx = output
    print(f'\nx* = {X}\nf(x*) = {fX}')
    error_dot = [get_norm(approx[i][0], dot_orig) for i in range(len(approx))]
    error_func = [abs(approx[i][1] - func_orig) for i in range(len(output[2]))]
    print('Погрешность точки', error_dot, '\nПогрешность функции', error_func)


