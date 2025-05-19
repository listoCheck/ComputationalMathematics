from math import sin, sqrt
from functools import reduce
from math import factorial
import numpy as np
from matplotlib import pyplot as plt


def lagrange(xs, ys, n):
    return lambda x: sum([
        ys[i] * reduce(
            lambda a, b: a * b,
            [(x - xs[j]) / (xs[i] - xs[j])
             for j in range(n) if i != j])
        for i in range(n)])


def divided_differences(x, y):
    n = len(y)
    coef = np.copy(y).astype(float)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) / (x[i] - x[i - j])
    return coef


def newton_divided_difference_polynomial(xs, ys, n):
    coef = divided_differences(xs, ys)
    return lambda x: ys[0] + sum([
        coef[k] * reduce(lambda a, b: a * b, [x - xs[j] for j in range(k)]) for k in range(1, n)
    ])


def finite_differences(y):
    n = len(y)
    delta_y = np.zeros((n, n))
    delta_y[:, 0] = y
    for j in range(1, n):
        for i in range(n - j):
            delta_y[i, j] = delta_y[i + 1, j - 1] - delta_y[i, j - 1]
    return delta_y


def print_finite_differences_table(delta_y):
    n = delta_y.shape[0]
    print("Таблица конечных разностей:")
    for i in range(n):
        row = [f"{delta_y[i, j]:.4f}" if i + j < n else "" for j in range(n)]
        print("\t".join(row))


def find_fin_difs(ys, n):
    fin_difs = [ys[:]]
    for i in range(1, n):
        prev = fin_difs[-1]
        fin_difs.append([prev[j + 1] - prev[j] for j in range(len(prev) - 1)])
    return fin_difs


def gauss(xs, ys, n):
    n = len(xs) - 1
    alpha_ind = n // 2
    fin_difs = [ys[:]]

    for k in range(1, n + 1):
        last = fin_difs[-1][:]
        fin_difs.append(
            [last[i + 1] - last[i] for i in range(n - k + 1)])

    h = xs[1] - xs[0]
    dts1 = [0, -1, 1, -2, 2, -3, 3, -4, 4]
    f1 = lambda x: ys[alpha_ind] + sum([
        reduce(lambda a, b: a * b,
               [(x - xs[alpha_ind]) / h + dts1[j] for j in range(k)])
        * fin_difs[k][len(fin_difs[k]) // 2] / factorial(k)
        for k in range(1, n + 1)])

    f2 = lambda x: ys[alpha_ind] + sum([
        reduce(lambda a, b: a * b,
               [(x - xs[alpha_ind]) / h - dts1[j] for j in range(k)])
        * fin_difs[k][len(fin_difs[k]) // 2 - (1 - len(fin_difs[k]) % 2)] / factorial(k)
        for k in range(1, n + 1)])

    return lambda x: f1(x) if x > xs[alpha_ind] else f2(x)


def stirling_polynomial(xs, ys, n):
    print("Делаем интерполяцию Стирлинга, количество точек нечетное")
    assert n % 2 == 1, "Для интерполяции Стирлинга нужно нечетное число точек"
    h = xs[1] - xs[0]
    m = n // 2
    a = xs[m]

    fin_difs = find_fin_difs(ys, n)

    def stirling(x):
        t = (x - a) / h
        result = ys[m]
        fact = 1
        t2 = 1

        for k in range(1, n):
            fact *= k

            if k % 2 == 1:
                idx = m - k // 2
                term = (t * t2 * (fin_difs[k][idx] + fin_difs[k][idx - 1]) / 2) / fact
                t2 *= (t ** 2 - ((k // 2) ** 2))
            else:
                idx = m - k // 2 - 1
                term = (t2 * fin_difs[k][idx]) / fact
                t2 *= (t ** 2 - ((k // 2 - 1) ** 2))
            result += term

        return result

    return stirling


def bessel_polynomial(xs, ys, n):
    print("Делаем интерполяцию Бесселя, колиество точек четное")
    assert n % 2 == 0, "Для интерполяции Бесселя нужно чётное число точек"
    h = xs[1] - xs[0]
    m = n // 2
    a = (xs[m - 1] + xs[m]) / 2

    fin_difs = find_fin_difs(ys, n)

    def bessel(x):
        t = (x - a) / h
        result = (ys[m] + ys[m - 1]) / 2
        fact = 1
        t2 = 1
        for k in range(1, n):
            fact *= k
            if k % 2 == 1:
                idx = m - (k // 2) - 2
                delta1 = fin_difs[k][idx]
                delta2 = fin_difs[k][idx + 1]
                term = (t * t2 * (delta1 + delta2) / 2) / fact
                t2 *= (t ** 2 - (k // 2) ** 2)
            else:
                idx = m - (k // 2)
                term = (t2 * fin_difs[k][idx]) / fact
                t2 *= (t ** 2 - (k // 2) ** 2)
            result += term
        return result

    return bessel


def draw_plot(a, b, func, name, dx=0.001):
    xs, ys = [], []
    a -= dx
    b += dx
    x = a
    while x <= b:
        xs.append(x)
        ys.append(func(x))
        x += dx
    plt.plot(xs, ys, 'g', label=name)


def solve(xs, ys, x, n):
    delta_y = finite_differences(ys)
    print_finite_differences_table(delta_y)

    print('\n' + '-' * 60)

    methods = [("Многочлен Лагранжа", lagrange),
               ("Многочлен Ньютона с разделенными разностями", newton_divided_difference_polynomial),
               ("Многочлен Гаусса", gauss),
               ("Многочлен Стирлинга", stirling_polynomial),
               ("Многочлен Бесселя", bessel_polynomial)]

    for name, method in methods:
        if (method is stirling_polynomial) and len(xs) % 2 == 0:
            print("убрать последнюю точку, чтобы посчитать метод Стирлинга? (y/n)")
            if 'y' == input():
                newx = xs
                newy = ys
                p2 = stirling_polynomial(newx, newx, n - 1)
                h = newx[1] - newx[0]
                alpha_ind = (n + 1) // 2
                newt = (x - newx[alpha_ind]) / h
                print("t: ", f'{newt:.5f}')
                print(name)
                print(f'P({x}) = {p2(x):.6f}')
                print('-' * 60)

                plt.title(name)
                draw_plot(newx[0], newx[-1], p2, name)
                plt.xlabel("X")
                plt.ylabel("Y")
                plt.scatter(x, p2(x), c='r')
                for i in range(len(newx)):
                    plt.scatter(newx[i], newy[i], c='b')

                plt.show()
                continue
            else:
                continue
        if method is bessel_polynomial and len(xs) % 2 == 1:
            print("убрать последнюю точку, чтобы посчитать метод Бесселя? (y/n)")
            if 'y' == input():
                newx = xs
                newy = ys
                p2 = bessel_polynomial(newx, newy, n - 1)
                h = newx[1] - newx[0]
                alpha_ind = (n + 1) // 2
                newt = (x - newx[alpha_ind]) / h
                print("t: ", f'{newt:.5f}')
                print(name)
                print(f'P({x}) = {p2(x):.6f}')
                print('-' * 60)
                plt.title(name)
                draw_plot(newx[0], newx[-1], p2, name)
                plt.xlabel("X")
                plt.ylabel("Y")
                plt.scatter(x, p2(x), c='r')
                for i in range(len(newx)):
                    plt.scatter(newx[i], newy[i], c='b')
                plt.show()
                continue
            else:
                continue

        h = xs[1] - xs[0]
        alpha_ind = n // 2
        t = (x - xs[alpha_ind]) / h
        print("t: ", f'{t:.5f}')

        print(name)
        P = method(xs, ys, n)
        print(f'P({x}) = {P(x):.6f}')
        print('-' * 60)

        plt.title(name)
        draw_plot(xs[0], xs[-1], P, name)
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.scatter(x, P(x), c='r')
        for i in range(len(xs)):
            plt.scatter(xs[i], ys[i], c='b')

        plt.show()


def read_data_from_file(filename):
    try:
        with open(filename, 'r') as file:
            xs = []
            ys = []
            x_read = False
            for line in file:
                if not x_read:
                    x_read = True
                    x = float(line.strip())
                else:
                    point = line.strip().split()
                    if len(point) == 2:
                        xs.append(float(point[0]))
                        ys.append(float(point[1]))

        return x, xs, ys, None
    except IOError as err:
        return None, None, None, "! Невозможно прочитать файл {0}: {1}".format(filename, err)


def read_data_from_input():
    x = float(input("Введите точку интерполяции: "))
    str = ''
    xs = []
    ys = []
    print("Введите 'quit', чтобы закончить ввод.")
    print("Введите узлы интерполяции:")
    while str != 'quit':
        str = input()
        point = str.strip().split()
        if len(point) == 2:
            xs.append(float(point[0]))
            ys.append(float(point[1]))
        else:
            if str != 'quit':
                print("! Неправильный ввод. Введенная точка не будет использована.")
    return x, xs, ys


def read_data_from_example():
    x = 0.502
    # x = 0.645
    xs = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8]
    ys = [1.5320, 2.5356, 3.5406, 4.5462, 5.5504, 6.5559, 7.5594]
    return x, xs, ys


def read_data_from_function():
    print('Функции: ')
    print('1. 2*x^2 - 5*x')
    print('2. x^5')
    print('3. sin(x)')
    print('4. sqrt(x)')

    while True:
        input_func = int(input('Выберите функцию [1/2/3]: '))
        if input_func == 1:
            f = lambda x: 2 * x ** 2 - 5 * x
            break
        elif input_func == 2:
            f = lambda x: x ** 5
            break
        elif input_func == 3:
            f = lambda x: sin(x)
            break
        elif input_func == 4:
            f = lambda x: sqrt(x)
            break
        else:
            print("! Некорректный ввод. Попробуйте еще раз\n")

    n = int(input("Введите число узлов: "))
    x0 = float(input('Введите x0: '))
    xn = float(input('Введите xn: '))

    h = (xn - x0) / (n - 1)
    xs = [x0 + h * i for i in range(n)]
    ys = list(map(f, xs))

    x = float(input('Введите точку интерполяции: '))
    return x, xs, ys


def main():
    while True:
        while True:
            print(
                "Введите: 'fi' для ввода из файла; 'e' для задания вычислительной части; 't' для ввода с терминала; 'fu' для задания функции.")
            option = input("Ваш ввод: ")
            if option == 'fi':
                while True:
                    filename = input("Введите имя файла: ")
                    x, xs, ys, error = read_data_from_file(filename)
                    if error != None:
                        print(error)
                        one_more_time = input("Вы хотите попробовать другое имя файла? [y/n]: ")
                        if one_more_time == 'y':
                            continue
                        else:
                            print('Ввод с клавиатуры:')
                            x, xs, ys = read_data_from_input()
                            break
                    else:
                        break
                n = len(xs)
                break
            elif option == 't':
                x, xs, ys = read_data_from_input()
                n = len(xs)
                break
            elif option == 'fu':
                x, xs, ys = read_data_from_function()
                n = len(xs)
                break
            elif option == 'e':
                x, xs, ys = read_data_from_example()
                n = len(xs)
                break
            else:
                print("! Некорректный ввод. Попробуйте еще раз\n")

        if len(set(xs)) != len(xs):
            print('! Узлы интерполции не должны совпадать. Введите еще раз.')
        elif xs != sorted(xs):
            print('! X интерполяции должны быть отсортированы. Введите еще раз.')
        else:
            break
    print("иксы:   ", *xs, sep="  ")
    print("угреки: ", *ys, sep="  ")
    solve(xs, ys, x, n)


if __name__ == "__main__":
    main()
