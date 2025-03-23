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


def system3(xy):
    x, y = xy
    return [math.exp(x) + y ** 2 - 3, x ** 2 + y - 1]


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


def phi1_system1(x, y):
    return math.sqrt(1 - y ** 2)


def phi2_system1(x, y):
    return x ** 2 - 0.5


def phi1_system2(x, y):
    return math.sqrt((1 - 2 * y ** 2) ** 2)


def phi2_system2(x, y):
    return (math.atan(x ** 2) - 0.1) / x if x != 0 else 0


def phi1_system3(x, y):
    return math.sqrt(3 - math.pow(math.e, x))


def phi2_system3(x, y):
    return math.sqrt(1 - y)


def check_convergence(phi1, phi2, x0, y0):
    try:
        diff1_x = abs((phi1(x0 + 1e-5, y0) - phi1(x0, y0)) / 1e-5)
        diff1_y = abs((phi1(x0, y0 + 1e-5) - phi1(x0, y0)) / 1e-5)
        diff2_x = abs((phi2(x0 + 1e-5, y0) - phi2(x0, y0)) / 1e-5)
        diff2_y = abs((phi2(x0, y0 + 1e-5) - phi2(x0, y0)) / 1e-5)
        q = max(diff1_x + diff1_y, diff2_x + diff2_y)
        return q < 1
    except:
        return False


def solve(system, phi1, phi2, x0, epsilon, max_iterations=1000):
    x = np.array(x0, dtype=float)
    if not check_convergence(phi1, phi2, x[0], x[1]):
        print("Метод простой итерации может не сойтись с выбранными начальными условиями!")
        print("продолжить? (Y/N)")
        if (input() == "N"):
            return None, None

    for iteration in range(max_iterations):
        x_next = np.array([phi1(x[0], x[1]), phi2(x[0], x[1])])
        print(f'{iteration}. x1={x_next[0]}, x2={x_next[1]}')
        if abs(system(x_next)[0]) < epsilon and abs(system(x_next)[1]) < epsilon:
            return x_next, iteration
        x = x_next
    print(f"Метод не сошелся за {max_iterations} итераций!")
    return x, max_iterations


def run():
    systems = {
        1: [system1, "x^2 + y^2 = 1, x^2 - y = 0.5", phi1_system1, phi2_system1],
        2: [system2, "tan(xy + 0.1) = x^2, x^2 + 2y^2 = 1", phi1_system2, phi2_system2],
        3: [system3, "e^x + y^2 = 3, x^2 + y = 1", phi1_system3, phi2_system3]
    }

    print("Доступные системы уравнений:")
    for key, (_, description, _, _) in systems.items():
        print(f"{key}: {description}")

    equations_number = int(input("Выберите систему уравнений (1, 2, 3): "))
    if equations_number not in systems:
        print("Неверный номер системы.")
        return

    plot_system(systems[equations_number][0])
    x0, y0 = map(float, input("Введите начальные приближения x0, y0: ").split())
    epsilon = float(input("Введите погрешность вычисления: "))
    xy_solution, iterations = solve(systems[equations_number][0], systems[equations_number][2],
                                    systems[equations_number][3], (x0, y0), epsilon)
    if xy_solution is not None:
        print(f"\nРешение: x = {xy_solution[0]:.5f}, y = {xy_solution[1]:.5f}, Итераций: {iterations}")


if __name__ == "__main__":
    run()
