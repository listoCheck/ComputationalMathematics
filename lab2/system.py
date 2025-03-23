import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def system1(xy):
    x, y = xy
    return [x ** 2 + y ** 2 - 1, x ** 2 - y - 0.5]

def system2(xy):
    x, y = xy
    return [math.tan(x * y + 0.1) - x ** 2, x ** 2 + 2 * y ** 2 - 1]

def f1(x, y):
    return math.tan(x * y + 0.1) - x ** 2

def f2(x, y):
    return x ** 2 + 2 * y ** 2 - 1

def plot_system(system):
    x = np.linspace(-2, 2, 400)
    y = np.linspace(-2, 2, 400)
    X, Y = np.meshgrid(x, y)

    Z1 = np.array([system([x_, y_])[0] for x_, y_ in zip(np.ravel(X), np.ravel(Y))]).reshape(X.shape)
    Z2 = np.array([system([x_, y_])[1] for x_, y_ in zip(np.ravel(X), np.ravel(Y))]).reshape(X.shape)

    plt.contour(X, Y, Z1, levels=[0], colors='r')
    plt.contour(X, Y, Z2, levels=[0], colors='b')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

def solve(system, phi1, phi2, x0, epsilon, max_iterations=1000):
    x = np.array(x0, dtype=float)

    try:
        for iteration in range(max_iterations):
            x_next = np.array([phi1(x[0], x[1]), phi2(x[0], x[1])])

            print(f'{iteration}. x1={x_next[0]}, x2={x_next[1]}')

            if abs(system(x_next)[0]) < epsilon and abs(system(x_next)[1]) < epsilon:
                return x_next, iteration

            x = x_next

        print(f"Метод простой итерации не сошелся за {max_iterations} итераций!")
        return x, max_iterations
    except ValueError as e:
        print(f'Ошибка: {e}')
        return None, None

def choose_system_of_equations(systems):
    print("Выберите систему уравнений:")
    for key, (system, description, _, _) in systems.items():
        print(f"{key}: {description}")

    while True:
        try:
            equations_number = int(input("Введите номер системы: "))
            if equations_number in systems:
                return equations_number
            else:
                print("(!) Такого номера нет.")
        except ValueError:
            print('(!) Вы ввели не число')

def phi1_system1(x, y):
    return math.sqrt(1 - y ** 2)

def phi2_system1(x, y):
    return x ** 2 - 0.5

def phi1_system2(x, y):
    if x == 0:
        return 0  # Избегаем деления на ноль
    return (math.atan(x ** 2) - 0.1) / x

def phi2_system2(x, y):
    value = (1 - 2 * y ** 2)**2
    if value < 0:
        return 0  # Избегаем вычисления квадратного корня из отрицательного числа
    return math.sqrt(value)

def run():
    systems = {
        1: [system1, "x^2 + y^2 = 1, x^2 - y = 0.5", phi1_system1, phi2_system1],
        2: [system2, "tan(xy + 0.1) = x^2, x^2 + 2y^2 = 1", phi1_system2, phi2_system2]
    }

    equations_number = choose_system_of_equations(systems)
    plot_system(systems[equations_number][0])

    x0, y0 = map(float, input("Введите начальные приближения x0, y0: ").split())
    epsilon = float(input('Введите погрешность вычисления: '))

    xy_solution, iterations = solve(systems[equations_number][0], systems[equations_number][2],
                                    systems[equations_number][3], (x0, y0), epsilon)

    if iterations is not None:
        print(f"\nНеизвестные: x = {xy_solution[0]:.5f}, y = {xy_solution[1]:.5f}")
        print(f"Количество итераций: {iterations}")


if __name__ == "__main__":
    run()
