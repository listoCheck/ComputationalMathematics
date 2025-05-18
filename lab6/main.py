from math import sin, exp, cos, inf
import numpy as np
from matplotlib import pyplot as plt
from tabulate import tabulate

MAX_ITERS = 20

def improved_euler(f, xs_orig, y0, eps, exact_y, n):
    name = "Усовершенствованный метод Эйлера"
    print('\n', name)

    inaccuracy = inf
    iters = 1
    ni = n
    x_start, x_end = xs_orig[0], xs_orig[-1]

    while inaccuracy > eps:
        h = (x_end - x_start) / ni
        xs = [x_start + i * h for i in range(ni + 1)]
        ys = [y0]

        for i in range(ni):
            f_pred = f(xs[i], ys[i])
            euler = ys[i] + h * f_pred
            y_next = ys[i] + h / 2 * (f_pred + f(xs[i + 1], euler))
            ys.append(y_next)

        h_half = h / 2
        xs_half = [x_start + i * h_half for i in range(2 * ni + 1)]
        ys_half = [y0]
        for i in range(2 * ni):
            f_pred = f(xs_half[i], ys_half[i])
            euler = ys_half[i] + h_half * f_pred
            y_next = ys_half[i] + h_half / 2 * (f_pred + f(xs_half[i + 1], euler))
            ys_half.append(y_next)

        runge_errors = []
        for i in range(1, ni + 1):
            y_h = ys[i]
            y_h2 = ys_half[i * 2]
            R = abs(y_h2 - y_h) / (2 ** 2 - 1)
            runge_errors.append(R)

        inaccuracy = max(runge_errors)
        if inaccuracy <= eps or iters >= MAX_ITERS:
            break
        iters += 1
        ni *= 2

    table = []
    for i in range(ni + 1):
        xi = xs[i]
        yi = ys[i]
        f_val = f(xi, yi)
        exact = exact_y(xi, x_start, y0)
        row = [
            i,
            f"{xi:.5f}",
            f"{yi:.5f}",
            f"{f_val:.5f}",
            f"{exact:.5f}"
        ]
        if i > 0:
            row.append(f"{runge_errors[i - 1]:.5e}")
        else:
            row.append("-")
        table.append(row)

    headers = ["i", "xi", "yi", "f(xi, yi)", "Точное решение", "Оценка погрешности (Рунге)"]
    print(tabulate(table, headers=headers, tablefmt="grid"))

    draw_approx_vs_exact(xs, ys, exact_y, name, x_start, y0)
    print(
        f"\nДля точности eps={eps} интервал был разбит на n={ni} частей с шагом h={round((xs[-1] - xs[0]) / ni, 6)} за {iters} итераций.\n")
    print("=" * 100, '\n')


def fourth_order_runge_kutt(f, xs_orig, y0, eps, exact_y, n):
    name = "Метод Рунге-Кутты 4-го порядка"
    print('\n', name)

    inaccuracy = inf
    iters = 1
    ni = n
    x_start, x_end = xs_orig[0], xs_orig[-1]

    while inaccuracy > eps:
        h = (x_end - x_start) / ni
        xs = [x_start + i * h for i in range(ni + 1)]
        ys = [y0]

        for i in range(ni):
            k1 = h * f(xs[i], ys[i])
            k2 = h * f(xs[i] + h / 2, ys[i] + k1 / 2)
            k3 = h * f(xs[i] + h / 2, ys[i] + k2 / 2)
            k4 = h * f(xs[i] + h, ys[i] + k3)
            y_next = ys[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
            ys.append(y_next)

        h_half = h / 2
        xs_half = [x_start + i * h_half for i in range(2 * ni + 1)]
        ys_half = [y0]
        for i in range(2 * ni):
            k1 = h_half * f(xs_half[i], ys_half[i])
            k2 = h_half * f(xs_half[i] + h_half / 2, ys_half[i] + k1 / 2)
            k3 = h_half * f(xs_half[i] + h_half / 2, ys_half[i] + k2 / 2)
            k4 = h_half * f(xs_half[i] + h_half, ys_half[i] + k3)
            y_next = ys_half[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
            ys_half.append(y_next)

        runge_errors = []
        for i in range(1, ni + 1):
            y_h = ys[i]
            y_h2 = ys_half[i * 2]
            R = abs(y_h2 - y_h) / (2 ** 4 - 1)
            runge_errors.append(R)

        inaccuracy = max(runge_errors)
        if inaccuracy <= eps or iters >= MAX_ITERS:
            break
        iters += 1
        ni *= 2

    table = []
    for i in range(ni + 1):
        xi = xs[i]
        yi = ys[i]
        f_val = f(xi, yi)
        exact = exact_y(xi, x_start, y0)
        row = [
            i,
            f"{xi:.5f}",
            f"{yi:.5f}",
            f"{f_val:.5f}",
            f"{exact:.5f}"
        ]
        if i > 0:
            row.append(f"{runge_errors[i - 1]:.5e}")
        else:
            row.append("-")
        table.append(row)

    headers = ["i", "xi", "yi", "f(xi, yi)", "Точное решение", "Оценка погрешности (Рунге)"]
    print(tabulate(table, headers=headers, tablefmt="grid"))

    draw_approx_vs_exact(xs, ys, exact_y, name, x_start, y0)
    print("=" * 100, '\n')


def milne(f, xs_orig, y0, eps, exact_y, n):
    global max_abs_error
    name = "Метод Милна"
    print("\n", name)

    inaccuracy = inf
    iters = 0
    ni = n
    x_start, x_end = xs_orig[0], xs_orig[-1]

    while inaccuracy > eps and iters < MAX_ITERS:
        h = (x_end - x_start) / ni
        xs = [x_start + i * h for i in range(ni + 1)]

        y = [y0]
        for i in range(1, 4):
            k1 = h * f(xs[i - 1], y[i - 1])
            k2 = h * f(xs[i - 1] + h / 2, y[i - 1] + k1 / 2)
            k3 = h * f(xs[i - 1] + h / 2, y[i - 1] + k2 / 2)
            k4 = h * f(xs[i - 1] + h, y[i - 1] + k3)
            y.append(y[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6)

        for i in range(4, len(xs)):
            yp = y[i - 4] + 4 * h * (
                2 * f(xs[i - 3], y[i - 3])
                - f(xs[i - 2], y[i - 2])
                + 2 * f(xs[i - 1], y[i - 1])
            ) / 3

            y_next = yp
            while True:
                yc = y[i - 2] + h * (
                    f(xs[i - 2], y[i - 2])
                    + 4 * f(xs[i - 1], y[i - 1])
                    + f(xs[i], y_next)
                ) / 3
                if abs(yc - y_next) < eps:
                    y_next = yc
                    break
                y_next = yc
            y.append(y_next)

        max_abs_error = 0
        for i in range(len(xs)):
            y_exact = exact_y(xs[i], x_start, y0)
            abs_error = abs(y_exact - y[i])
            max_abs_error = max(max_abs_error, abs_error)

        inaccuracy = max_abs_error
        if inaccuracy <= eps or iters >= MAX_ITERS:
            break
        iters += 1
        ni *= 2

    table = []
    for i in range(len(xs)):
        y_exact = exact_y(xs[i], x_start, y0)
        abs_error = abs(y_exact - y[i])
        row = [
            i,
            f"{xs[i]:.5f}",
            f"{y[i]:.5f}",
            f"{f(xs[i], y[i]):.5f}",
            f"{y_exact:.5f}",
            f"{abs_error:.5e}"
        ]
        table.append(row)

    headers = ["i", "xi", "yi", "f(xi, yi)", "Точное решение", "Абсолютная ошибка"]
    print(tabulate(table, headers=headers, tablefmt="grid"))
    if max_abs_error <= eps:
        print(f"Точность достигнута: погрешность {max_abs_error:.5e} <= {eps}")
    else:
        print(f"Точность НЕ достигнута: погрешность {max_abs_error:.5e} > {eps}")

    draw_approx_vs_exact(xs, y, exact_y, name, xs[0], y0)
    print("=" * 100, '\n')




"""def draw_plot(a, b, func, name, dx=0.001):
    xs, ys = [], []
    a -= dx
    b += dx
    x = a
    while x <= b:
        xs.append(x)
        ys.append(func(x))
        x += dx
    plt.plot(xs, ys, 'g', label=name)"""


def draw_approx_vs_exact(xs, ys, exact_y_func, name, x0, y0):
    plt.figure(figsize=(10, 6))
    plt.title(name)

    smooth_xs = np.linspace(xs[0], xs[-1], 500)
    smooth_ys = [exact_y_func(x, x0, y0) for x in smooth_xs]
    plt.plot(smooth_xs, smooth_ys, label='Точное решение', color='blue')

    plt.plot(xs, ys, 'o-', label='Приближенное решение', color='red')

    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.grid(True)
    plt.show()


def solve(f, x0, xn, y0, n, eps, exact_y):
    xs = [x0 + i * ((xn - x0) / n) for i in range(n + 1)]
    methods = [("Усовершенствованный метод Эйлера", improved_euler),
               ("Метод Рунге-Кутта 4-го порядка", fourth_order_runge_kutt),
               ("Метод Милна", milne)]
    for method in methods:
        ys = method[1](f, xs, y0, eps, exact_y, n)


def main():
    print("ОДУ:")
    print('1. y + (1 + x)*y^2')
    print('2. sin(x) - y')
    print('3. e^x\n')
    while True:
        try:
            input_func = int(input('> Выберите ОДУ [1/2/3]: '))
            if input_func == 1:
                f = lambda x, y: y + (1 + x) * y ** 2
                exact_y = lambda x, x0, y0: -exp(x) / (x * exp(x) - (x0 * exp(x0) * y0 + exp(x0)) / y0)
                break
            elif input_func == 2:
                f = lambda x, y: sin(x) - y
                exact_y = lambda x, x0, y0: (2 * exp(x0) * y0 - exp(x0) * sin(x0) + exp(x0) * cos(x0)) / (
                        2 * exp(x)) + (sin(x)) / 2 - (cos(x)) / 2
                break
            elif input_func == 3:
                f = lambda x, y: exp(x)
                exact_y = lambda x, x0, y0: y0 - exp(x0) + exp(x)
                break
            else:
                print("! Некорректный ввод. Попробуйте еще раз\n")
        except:
            print("! Некорректный ввод. Попробуйте еще раз\n")
    while True:
        try:
            x0 = float(input('> Введите первый элемент интервала x0: '))
            xn = float(input('> Введите последний элемент интервала xn: '))
            n = int(input('> Введите количество элементов в интервале n: '))

            y0 = float(input('> Введите y0: '))
            eps = float(input('> Введите точность eps: '))

            if xn <= x0:
                print('! xn должен быть больше x0. Введите еще раз.')
            elif n <= 1:
                print('! Количество элементов n должно быть > 1. Введите еще раз.')
            else:
                break
        except:
            print("! Некорректный ввод. Попробуйте еще раз\n")
    solve(f, x0, xn, y0, n, eps, exact_y)


if __name__ == "__main__":
    main()
