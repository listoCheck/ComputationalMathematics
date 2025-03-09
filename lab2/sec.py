def func(n):
    return -1.38 * (n ** 3) - 5.42 * (n ** 2) + 2.57 * n + 10.95


x_p = -1.5
x_n = -1.0
for i in range(10):
    x_f = x_n - (x_n - x_p)/(func(x_n) - func(x_p))*func(x_n)
    print(
        f"{round(x_p, 10):.7f}",
        f"{round(x_n, 10):.7f}",
        f"{round(x_f, 10):.7f}",
        f"{round(func(x_f), 10):.7f}",
        f"{round(abs(x_f - x_n), 10):.7f}"
    )
    x_p, x_n = x_n, x_f
