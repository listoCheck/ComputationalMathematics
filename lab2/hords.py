def func(n):
    return -1.38 * (n ** 3) - 5.42 * (n ** 2) + 2.57 * n + 10.95


a = -4.0
b = -2.0
for i in range(6):
    f_a = func(a)
    f_b = func(b)
    x = (a * f_b - b * f_a) / (f_b - f_a)
    print(
        f"{round(a, 10):.5f}",
        f"{round(b, 10):.5f}",
        f"{round(x, 10):.5f}",
        f"{round(f_a, 10):.5f}",
        f"{round(f_b, 10):.5f}",
        f"{round(func(x), 10):.5f}",
        f"{round(x - b, 10):.5f}"
    )
    b = x
