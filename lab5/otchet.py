a = [1.5320, 2.5356, 3.5406, 4.5462, 5.5504, 6.5559, 7.5594]
ans = []
for _ in range(6):
    newArray = []
    for i in range(len(a) - 1):
        newArray.append(a[i + 1] - a[i])
    a = newArray
    ans.append(a)
for i in range(6):
    for j in range(6):
        try:
            print(f'{ans[j][i]:.4f}',  end=";")
        except:
            print(";", end="")
    print()
