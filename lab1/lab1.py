def sordDiag():
    for i in range(h):
        elements = []
        for j in range(h):
            if j < i:
                elements.append(0)
            else:
                elements.append(abs(matrix[j][i]))
        matrix[i], matrix[elements.index(max(elements))] = matrix[elements.index(max(elements))], matrix[i]


def firstTheorem():
    for i in range(h):
        if 2 * abs(matrixC[i][i]) < abs(sum(matrixC[i])):
            return False
    return True


def toDiag():
    for i in range(h):
        delitel = matrixC[i][i]
        for j in range(h):
            if i == j:
                matrixC[i][j] = 0
            else:
                matrixC[i][j] /= -delitel
        matrixD[i] /= delitel


def secondTheorem():
    for i in range(h):
        if abs(sum(matrixC[i])) <= eps:
            print("False")
            return False
    return True


def printMatrix():
    for i in range(h):
        print(matrix[i])


def makeMatrixs():
    for i in range(h):
        forC = []
        for j in range(h):
            forC.append(matrix[i][j])
        matrixC.append(forC)
        matrixD.append(matrix[i][-1])


def start():
    for i in range(h):
        matrix.append([float(j) for j in input().split()])

    print("Исходная матрица:")
    printMatrix()

    sordDiag()

    print("Преобразованная матрица:")
    printMatrix()

    makeMatrixs()

    if not firstTheorem():
        print("Диагональное преобладание невозможно!")
        return False

    toDiag()

    if not secondTheorem():
        print("Метод не сходится!")
        return False

    matrixK = []
    prevMatrixK = [i for i in matrixD]
    iterations = 1

    while iterations < max_iter:
        for i in range(h):
            row = 0
            for j in range(h):
                row += matrixC[i][j] * prevMatrixK[j]
            row += matrixD[i]
            matrixK.append(row)

        delt = 0
        for i in range(h):
            delt = max(abs(matrixK[i] - prevMatrixK[i]), delt)

        if delt < eps:
            print(f"Решение найдено за {iterations} итераций")
            print("Вектор неизвестных:", matrixK)
            print("Вектор погрешностей:", [abs(matrixK[i] - prevMatrixK[i]) for i in range(h)])
            return True
        else:
            prevMatrixK = matrixK
            matrixK = []
            iterations += 1

    print("Достигнуто максимальное число итераций. Решение могло не сойтись.")
    return False


eps = float(input("Введите точность: "))
max_iter = int(input("Введите максимальное число итераций: "))
h = int(input("Введите размерность матрицы: "))

matrix = []
matrixC = []
matrixD = []

if start():
    print("Метод успешно завершился.")
else:
    print("Решение не найдено.")
