"""
Пакет pract3 решает задачи практической 3.

FUNCTIONS:
    main() - основная работа программы

Latipov, Amelin, Zueva, Babayan
# License: BSD
"""

import copy
import numpy
import fractions
import random
import csv


#Тестовые значения функций
MATRIX = [
    [20, 2, 3, 7],
    [1, 12, -2, -5],
    [5, -3, 13, 0],
    [0, 0, -1, 15]
    ]
FREE = [5, 4, -3, 7]


def print_matrix(matrix):
    """Форматный вывод на экран"""
    for row in matrix:
        print(row)

def get_input_from_keyboard():
    """Получение значений с клавиатуры"""
    while True:
        length = input('Введите размерность матрицы')
        try:
            length = int(length)
            break
        except:
            print("Неверный ввод")
    matrix = [[None for _ in range(length)] for __ in range(length)]
    while True:
        try:
            for i in range(length):
                for j in range(length):
                    element = input(f"Введите a{i+1}{j+1}: ")
                    matrix[i][j] = complex(element)
            break
        except:
            print("Неверный ввод")
    free = [None for _ in range(length)]
    while True:
        try:     
            for i in range(length):
                element = input(f"Введите b{i+1}: ")
                free[i] = complex(element)
            break
        except:
            print("Неверный ввод")
    return matrix, free

def get_random():
    """Сгенерировать случайные значения матрицы"""
    length = random.randint(2,10)
    matrix = [[None for _ in range(length)] for __ in range(length)]
    free = [None for _ in range(length)]
    for i in range(length):
        free[i] = random.random()*10-5
        for j in range(length):
            matrix[i][j] = random.random()*10-5
    return matrix, free

def get_from_csv(filename): 
    """
    Получить значения из csv файла.
    
    Обязательный аргумент - имя файла без .csv"""
    matrix = []
    col = 0
    modify = float
    with open(f'{filename}.csv', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=';')
        for row in spamreader:
            if col:
                if len(row) != col:
                    print("Ошибка считывания файла")
                    return None
            elif len(row)>1:
                col = len(row)
            else:
                print("Ошибка в матрице 1")
                return None
            
            matrix.append(row)
            for item in row:
                if "j" in item:
                    modify = complex
                    
    if len(matrix) != 0 and (len(matrix) + 1 == len(matrix[0])):
        free = []
        try:
            for i in range(len(matrix)):
                free.append(modify(matrix[i].pop()))
                for j in range(len(matrix[i])):
                    matrix[i][j] = modify(matrix[i][j])
        except:
            print("Ошибка в матрице 2")
            return None 
    else:
        print("Ошибка в матрице")
        return None
    return matrix, free

def jacobi(matrix = MATRIX, free = FREE):
    """
    Решение СЛАУ методом простых итераций Якоби.
    
    Передается матрица и массив свобдных членов
    matrix : list of lists
        матрица коэффициентов
    free : list
        массив свободных членов
    """
    matrix = copy.deepcopy(matrix)
    free = copy.deepcopy(free)
    def required_precision(x_old, x_new):
        eps = 0.0001
        sum_up = 0
        sum_low = 0
        for k in range(0, len(x_old)):
            sum_up += ( x_new[k] - x_old[k] ) ** 2
            sum_low += ( x_new[k] ) ** 2
        sum_low = sum_low if sum_low!=0 else 10**(-100)
        return abs(( sum_up / sum_low )**(1/2)) < eps

    def check_null(matrix, free):
        length = len(matrix)
        for row in range(length):
            if( matrix[row][row] == 0 ):
                for row_2 in range(length):
                    if matrix[row_2][row]!=0:
                        for i in range(length):
                            matrix[row][i] += matrix[row_2][i]
                        free[row] += free[row_2]
                        break
                else:
                    print('Неисправимо')
                    return None, None
        return matrix, free
    matrix, free = check_null(matrix, free)
    if not matrix:
        return None
    count = len(free) # количество корней
    x = copy.deepcopy(free) # начальное приближение корней
    iter_numb = 0  # подсчет количества итераций
    max_iter = 100    # максимально допустимое число итераций
    while( iter_numb < max_iter ):
        x_prev = copy.deepcopy(x)
        for k in range(count):
            summ = 0
            for j in range(count):
                if( j != k ):
                    summ = summ + matrix[k][j] * x[j]
            x[k] = free[k]/matrix[k][k] - summ / matrix[k][k]
            
        if required_precision(x_prev, x) : # проверка на выход
            break
        iter_numb += 1
        if iter_numb + 1 == max_iter:
            print(f"Требуемая точность не достигнута за {max_iter} итераций.")
    return matrix, numpy.linalg.inv(matrix), x

def jordan_gauss(matrix = MATRIX, free = FREE, fraction_convert = False):
    """
    Решение СЛАУ методом Гаусса-Жордана
    
    Передается матрица и массив свобдных членов
    matrix : list of lists
        матрица коэффициентов
    free : list
        массив свободных членов
    """
    length = len(matrix)
    free = copy.deepcopy(free)
    matrix = copy.deepcopy(matrix)
    returned = copy.deepcopy(matrix)
    returned2 = (numpy.linalg.inv(returned)).tolist()
    try:
        if fraction_convert:
            for i in range(length):
                free[i] = fractions.Fraction(free[i])
                for j in range(length):
                    matrix[i][j] = fractions.Fraction(matrix[i][j])
                    returned[i][j] = fractions.Fraction(returned[i][j])
                    returned2[i][j] = fractions.Fraction(returned2[i][j])
    except TypeError:
        print("Комплексные нельзя представить в виде дроби")
        return None
                
    for col in range(length):
        for row in range(length):
            if matrix[row][col] != 0 and sum(matrix[row][:col]) == 0:
                ratio = matrix[row][col]
                free[row] /= ratio
                for i in range(length):
                    matrix[row][i] /= ratio
                for i in range(length):
                    if i == row:
                        continue
                    free[i] -= free[row]*matrix[i][col]
                    ratio = matrix[i][col]
                    for j in range(length):
                        matrix[i][j] -= matrix[row][j]*ratio
                break
        else:
            print("Неисправимо")
            return None
        
    x = [None for _ in range(length)]
    for i in range(length):
        for j in range(length):
            if matrix[i][j] != 0:
                x[j] = free[i]
    return returned, returned2, x

def main(a = MATRIX, b = FREE):
    """
    Сравнение работы методов Якоби и Жордана-Гаусса
    
    Передается матрица и массив свобдных членов
    matrix : list of lists
        матрица коэффициентов
    free : list
        массив свободных членов
    """
    
    try:
        results = jacobi(a, b)
        cord = numpy.linalg.cond(a)
        print("Прямая матрица коэффициентов")
        print_matrix(results[0])
        print("Обратная матрица коэффициентов")
        print_matrix(results[1])
        print("Решение СЛАУ")
        print_matrix(results[2])
        print("Cord(A):",cord,"\nАлгоритм простых итераций Якоби\n", "-"*100)
    except OverflowError:
        cord = 101
        print("Не могу посчитать большие числа")
    if cord >=100:
        results = jordan_gauss(a, b)
        cord = numpy.linalg.cond(a)
        print("Прямая матрица коэффициентов")
        print_matrix(results[0])
        print("Обратная матрица коэффициентов")
        print_matrix(results[1])
        print("Решение СЛАУ")
        print_matrix(results[2])
        print("Cord(A):",cord,"\nПрямой алгоритм Гаусса-Жордана\n", "-"*100)
        if cord >=100:
            results = jordan_gauss(a, b, True)
            cord = numpy.linalg.cond(a)
            print("Прямая матрица коэффициентов")
            print_matrix(results[0])
            print("Обратная матрица коэффициентов")
            print_matrix(results[1])
            print("Решение СЛАУ")
            print_matrix(results[2])
            print("Cord(A):",cord,"\nПрямой алгоритм\
Гаусса-Жордана с правильными дробями\n", "-"*100)