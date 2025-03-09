num = float(input())
for i in range(5):
    func = -1.38 * (num ** 3) - 5.42 * (num ** 2) + 2.57 * num + 10.95
    new_num = num + func / 23.005
    func = -1.38 * (new_num ** 3) - 5.42 * (new_num ** 2) + 2.57 * new_num + 10.95
    print(
        f"{round(num, 10):.10f}",
        f"{round(new_num, 10):.10f}",
        f"{round(func, 10):.10f}",
        f"{round(abs(new_num - num), 10):.10f}"
    )
    num = new_num
