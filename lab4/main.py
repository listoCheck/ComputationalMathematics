import math
from sys import flags


def func(x):
    return (15 * x) / (x ** 4 + 2)


def func2(answer):
    ansStr = ""
    for i in answer:
        ansStr += str(i) + ";"
    print(ansStr)


def sxay(answer):
    sum = 0
    for i in answer:
        sum += float(i)
    print(f'{sum:.3f}')


def sxx(answer):
    sum = 0
    for i in answer:
        sum += float(i) ** 2
    print(f'{sum:.3f}')

def sxxx(answer):
    sum = 0
    for i in answer:
        sum += float(i) ** 3
    print(f'{sum:.3f}')

def sxxxx(answer):
    sum = 0
    for i in answer:
        sum += float(i) ** 4
    print(f'{sum:.3f}')


def sxy(answer):
    sum = 0
    for i in range(len(answer)):
        sum += float(answers1[i]) * float(answers2[i])
    print(f'{sum:.3f}')

def sxxy(answer):
    sum = 0
    for i in range(len(answer)):
        sum += float(answers1[i])**2 * float(answers2[i])
    print(f'{sum:.3f}')

def fi(x):
    return 0.141 * x + 1.492

def gamma(answer):
    sum = 0
    for i in answer:
        sum += float(i)
    print(f'{math.sqrt(sum / len(answer)):.3f}')

def fi2(x):
    return -3.223 * x ** 2 - 7.868 * x + 2.551

x = 0.0
answers0 = []
answers1 = []
answers2 = []
answersfi0 = []
answersfi1 = []
answersfi2 = []
answersfi3 = []
iter = 1
while x <= 4:
    answers0.append(iter)
    answers1.append(f'{x:.3f}')
    answers2.append(f'{func(x):.3f}')
    answersfi0.append(f'{fi(x):.3f}')
    answersfi1.append(f'{(float(answersfi0[-1]) - float(answers2[-1])) ** 2:.3f}')
    #answersfi2.append(f'{fi2(x):.3f}')
    #answersfi3.append(f'{(float(answersfi2[-1]) - float(answers2[-1])) ** 2:.3f}')
    x += 0.4
    iter += 1
func2(answers0)
func2(answers1)
func2(answers2)
func2(answersfi0)
func2(answersfi1)
#func2(answersfi2)
#func2(answersfi3)

sxay(answers1)
sxx(answers1)
sxxx(answers1)
sxxxx(answers1)

sxay(answers2)
sxy(answers2)
sxxy(answers2)

gamma(answersfi1)
