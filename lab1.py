def sordDiag():
    for i in range(h):
        elements = []
        for j in range(h):
            if j < i:
                elements.append(0)
            else:
                elements.append(abs(matrix[j][i]))
        print("el", elements)
        matrix[i], matrix[elements.index(max(elements))] = matrix[elements.index(max(elements))], matrix[i]
        print(elements.index(max(elements)))

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
    sordDiag()
    printMatrix()
    makeMatrixs()
    if firstTheorem():
        print("ok")
    else:
        print("bad")
        return False
    toDiag()
    if secondTheorem():
        print("ok")
    else:
        print("bad")
        return False
    matrixK = []
    prevMatrixK = [i for i in matrixD]
    while True:
        for i in range(h):
            row = 0
            for j in range(h):
                row += matrixC[i][j] * prevMatrixK[j]
            row += matrixD[i]
            matrixK.append(row)
        delt = 0
        print("K", matrixK)
        print("PK", prevMatrixK)
        for i in range(h):
            delt = max(abs(matrixK[i] - prevMatrixK[i]), delt)
        print("delt", delt)
        if delt < eps:
            return True
        else:
            print("-------------------------------<next>-----------------------------------")
            prevMatrixK = matrixK
            matrixK = []



eps = float(input())
h = int(input())
matrix = []
matrixC = []
matrixD = []
if start():
    print("ok")
else:
    print("нас рано в значениях")
