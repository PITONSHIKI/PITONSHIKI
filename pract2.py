"""
Пакет pract2 решает задачи практической 2.

FUNCTIONS:
    main() - основная работа программы

Latipov, Amelin, Zueva, Babayan
# License: BSD
"""

import numpy as np
import pandas as pd
from datetime import datetime
import random
import csv
import sys

def main():
    """
    Сравнение скорости выполнения операций с матрицами.
    
    Вычисление определителя собственной функцией и транспонирование.
    Вычисление определителя и транспонирование встроенной функцией.
    Сравнение результатов
    """
    #Генерирует матрицу с размерностями axb и элементами в диапазоне [c;d]
    def matrix_generator(a, b, c, d):
        """Генерирует матрицу с размерностями axb и элементами в диапазоне [c;d]"""
        help = []
        for i in range(a):
            st = []
            for j in range(b):
                st.append(random.random()*(d - c) + c)#Пользуемся генератором
            help.append(st)
        return np.array(help)

    #Функция вычисления ранга
    def rank(a):
        """Вычисление ранга"""
        return np.linalg.matrix_rank(a)

    #Транспонирование матрицы
    def transpon(A):
        """Транспонирование матрицы"""
        return np.array([[A[j][i] for j in range(len(A))] for i in range(len(A[0]))])

    #Вычисление определителя
    def det(A):
        """Вычисление определителя"""
        #Частный случай
        if(len(A) != len(A[0])):
            print("В определитель передана неквадратная матрица")
            sys.exit()
            
        #Выход из рекурсии
        if(rank(A) < 2): 
            if(len(A) != 1): 
                return 0
            else:
                return A[0][0]
        
        #Вычисление определителя
        ans = 0
        for k in range(len(A)):
            help = []
            for i in range(1, len(A[0])):
                h = []
                for j in range(len(A)):
                    if(j != k):
                        h.append(A[i][j])
                help.append(h)
            ans += (-1)**k * det(help) * A[0][k]
        return ans

    #Ввод размерностей с проверкой
    print("Введите два числа, размерности матрицы")
    a, b = 0, 0
    while(a < 1 or b < 1):
        a, b = list(map(int, input().split()))
        if(a < 1 or b < 1):
            print("Размерности некорректны, они должны быть\
    целыми положительными числами, повторите попытку")
        
    #Ввод матрицы с проверками
    k = 0
    print("Введите матрицу:")
    while(k == 0):
        k = 1
        A = []
        for i in range(a):
            A.append(list(map(float, input().split())))

        for i in range(a):
            if(len(A[i]) != b):
                k = 0
        if(k == 0):
            print("Ошибка ввода матрицы, возможно кол-во чисел неравное\
    в разных строках, или их длина не соответстует заданной, повторите попытку:")

    A = np.array(A)

    #Ввод матрицы из файла с проверками
    B = []
    print("В введенном файле матрица должна иметь разделители - запятые")
    file = input("Введите название файла с расширением csv: ")
    with open(file) as File:
        reader = csv.reader(File, delimiter=',', quotechar=',',
                            quoting=csv.QUOTE_MINIMAL)
        for row in reader:
            B.append(list(map(float, row)))

    k = 1
    for i in range(len(B)):
        if(len(B[i]) != len(B[0])):
            k = 0
        if(k == 0):
            print("Ошибка ввода матрицы, возможно кол-во\
    чисел неравное в разных строках, некорректный файл:")
            sys.exit()
    B = np.array(B)

    #Генерация матрицы
    print("Введите размерности матрицы и диапазон внутренних\
    значений, чтобы сгенерировать матрицу")
    a, b, c, d = list(map(int, input().split()))
    C = matrix_generator(a, b, c, d)

    #Размер таблицы
    nmax = int(input("Введите кол-во тестов: "))

    print("Сравним скорость работы системного и нашего траспонирования:")
    #Формируем тесты
    TEST = []
    for i in range(nmax):
        TEST.append(
            matrix_generator(random.randint(100, 200),
            random.randint(100, 200), c, d)
            )

    #Генерируем списки с записанным временем работы
    lrun = []
    rrun = []
    drun = []
    for i in range(nmax):
        x = datetime.now()
        transpon(TEST[i])
        lrun.append(datetime.now() - x)
        x = datetime.now()
        TEST[i].transpose()
        rrun.append(datetime.now() - x)
        drun.append(lrun[-1] - rrun[-1])

    #Создаем таблицу
    date = pd.DataFrame({
        "Our tranpose" : lrun,
        "System transpose" : rrun,
        "delta" : drun
        })
    print(date)
    print()

    #Аналогичные действия
    print("Сравним скорость работы системного и нашего вычисления определителя:")
    TEST = []
    for i in range(nmax):
        k = random.randint(1, 5)
        TEST.append(matrix_generator(k, k, c, d))

    lrun = []
    rrun = []
    drun = []
    for i in range(nmax):
        x = datetime.now()
        det(TEST[i])
        lrun.append(datetime.now() - x)
        x = datetime.now()
        np.linalg.det(TEST[i])
        rrun.append(datetime.now() - x)
        drun.append(lrun[-1] - rrun[-1])

    date = pd.DataFrame({
        "Our det" : lrun,
        "System det" : rrun,
        "delta" : drun
        })
    print(date)
    print("A = \n", A)
    print()
    print("B = \n", B)
    print()
    print("C = \n", C)
    print("Введите матричное выражение c А, B и C, используя такие сокращения:")
    print("rank(A) - ранг матрицы")
    print("transpon(A) - транспонирование матрицы")
    print("det(A) - определитель матрицы")
    s = input()
    print(eval(s))