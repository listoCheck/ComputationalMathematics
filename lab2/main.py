import math
import numpy as np
import matplotlib.pyplot as plt


def equation_1(x):
    return -1.38 * x ** 3 - 5.42 * x ** 2 + 2.57 * x + 10.95


def equation_2(x):
    return x ** 3 - 1.89 * x ** 2 - 2 * x + 1.76


def equation_3(x):
    return 2 * x ** 3 + 3.41 * x ** 2 - 23.74 * x + 2.95


def derivative_1(x):
    return -4.14 * x ** 2 - 10.84 * x + 2.57


def derivative_2(x):
    return 3 * x ** 2 - 3.78 * x - 2


def derivative_3(x):
    return 6 * x ** 2 + 6.82 * x - 23.74


def second_derivative_1(x):
    return -8.28 * x - 10.84


def second_derivative_2(x):
    return 6 * x - 3.78


def second_derivative_3(x):
    return 12 * x + 6.82


def system_1(xy):
    x, y = xy
    return [x ** 2 + y ** 2 - 1, x ** 2 - y - 0.5]


def phi1(x, y):
    return math.sqrt(1 - y ** 2)


def phi2(x, y):
    return x ** 2 - 0.5


def solve_system(system, phi1, phi2, x0, epsilon, max_iterations=1000):
    x = np.array(x0, dtype=float)

    for iterations in range(max_iterations):
        x1 = phi1(x[0], x[1])
        x2 = phi2(x[0], x[1])
        x_next = np.array([x1, x2])

        print(f'{iterations}. x1={x1}, x2={x2}, |xk+1 - xk|={np.linalg.norm(x_next - x)}')

        if abs(system(x_next)[0]) < epsilon and abs(system(x_next)[1]) < epsilon:
            return x_next, iterations

        x = x_next

    print("Метод не сошелся за заданное число итераций.")
    return None, None


def get_user_input_system():
    x0, y0 = map(float, input("Введите начальные приближения x0, y0: ").split())
    epsilon = float(input("Введите погрешность вычисления: "))
    return (x0, y0), epsilon


def plot_function(f, a, b):
    x = np.linspace(a, b, 400)
    y = [f(xi) for xi in x]
    plt.plot(x, y, label='Функция')
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(True, linestyle='--', linewidth=0.5)
    plt.legend()
    plt.show()


def bisection_method(f, a, b, eps):
    if f(a) * f(b) > 0:
        print("Ошибка: на данном интервале нет корня или их несколько.")
        return
    plot_function(f, a, b)
    n = 0
    while abs(b - a) > eps:
        x = (a + b) / 2
        if f(x) == 0:
            break
        elif f(a) * f(x) < 0:
            b = x
        else:
            a = x
        n += 1
        print(f"Найденный корень: {x}, Значение функции в корне: {f(x)}, Число итераций: {n}")

    return


def newton_method(f, df, a, b, eps, ddf):
    if f(a) * ddf(a) > 0:
        x0 = a
    elif f(b) * ddf(b) > 0:
        x0 = b
    else:
        print("Ошибка: на данном интервале нет корня или их несколько.")
        return
    plot_function(f, a, b)
    n = 0
    while True:
        x = x0 - f(x0) / df(x0)
        print(
            f"Найденный корень: {x}, Значение функции в x: {f(x)}, Значение производной функции в x: {df(x)}, Число итераций: {n}")
        if (x - x0) <= eps or abs(f(x) / df(x)) <= eps or abs(f(x)) <= eps:
            break
        n += 1
        x0 = x
    return


def simple_iteration_method(f, df, a, b, eps, ddf):
    def phi(x):
        return x - f(x) / max(df(a), df(b))

    if f(a) * ddf(a) > 0:
        x0 = a
    elif f(b) * ddf(b) > 0:
        x0 = b
    else:
        print("Ошибка: не удается выбрать начальное приближение.")
        return

    print(f"Метод простой итерации. Начальное приближение: x0 = {x0}")

    n = 0
    # for i in range(5):
    while True:
        x1 = phi(x0)
        print(f"Итерация {n}: x0 = {x0}, x1 = {x1}, |x1 - x0| = {abs(x1 - x0)}")
        if abs(x1 - x0) < eps:
            break
        x0 = x1
        n += 1

    print(f"Найденный корень: {x1}, Число итераций: {n}")
    return


def get_input_from_file(filename):
    try:
        with open(filename, 'r') as file:
            data = file.readlines()
            return [float(value.strip()) for value in data]
    except Exception as e:
        print(f"Ошибка чтения файла: {e}")
        return None


def get_user_input():
    choice = input("Введите имя файла для загрузки данных или оставьте пустым для ручного ввода: ").strip()
    if choice:
        return get_input_from_file(choice)
    else:
        return [
            float(input("Введите левую границу интервала: ")),
            float(input("Введите правую границу интервала: ")),
            float(input("Введите погрешность вычисления: "))
        ]


def main():
    while True:
        print("Выберите тип программы:")
        print("1: Нелинейное уравнение")
        print("2: Система нелинейных уравнений")
        print("3: Выход")
        prog_type = input("Введите номер типа: ")
        if prog_type == "3":
            break
        if prog_type == "2":
            (x0, y0), epsilon = get_user_input_system()
            solution, iterations = solve_system(system_1, phi1, phi2, (x0, y0), epsilon)
            if solution is not None:
                print(f"Решение: x = {solution[0]:.5f}, y = {solution[1]:.5f}, итерации: {iterations}")
            continue
        print("Выберите уравнение:")
        print("1: -1.38*x^3 - 5.42*x^2 + 2.57*x + 10.95")
        print("2: x^3 - 1.89*x^2 - 2*x + 1.76")
        print("3: 2*x^3 + 3.41*x^2 - 23.74*x + 2.95")
        eq_choice = input("Введите номер уравнения: ")

        equations = {"1": equation_1, "2": equation_2, "3": equation_3}
        derivatives = {"1": derivative_1, "2": derivative_2, "3": derivative_3}
        second_derivatives = {"1": second_derivative_1, "2": second_derivative_2, "3": second_derivative_3}
        f = equations.get(eq_choice)
        df = derivatives.get(eq_choice)
        ddf = second_derivatives.get(eq_choice)

        if not f:
            print("Некорректный выбор уравнения.")
            continue

        print("Выберите метод:")
        print("1: Метод половинного деления")
        print("2: Метод Ньютона")
        print("3: Метод простой итерации")
        method_choice = input("Введите номер метода: ")

        input_data = get_user_input()
        if not input_data:
            continue

        if method_choice == "1":
            bisection_method(f, input_data[0], input_data[1], input_data[2])
        elif method_choice == "2":
            newton_method(f, df, input_data[0], input_data[1], input_data[2], ddf)
        elif method_choice == "3":
            simple_iteration_method(f, df, input_data[0], input_data[1], input_data[2], ddf)
        else:
            print("Некорректный выбор метода.")
            continue

        again = input("Еще раз? [y/n]: ")
        if again.lower() != "y":
            break


main()
