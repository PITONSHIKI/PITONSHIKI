"""
Пакет pract4 решает задачи практической 4.

FUNCTIONS:
    lagrange(csv_file) - вычисление методом Лагранжа
    newton(csv_file) - вычисление методом Ньютона
    linear(csv_file) - Апрксимация линейной функцией
    quadr(csv_file) - Апроксимация квадратичной функцией
    norm(scv_file) - Апроксимация функцией нормального распределния
    comparision (csv_file) - сравнение всех методов

Latipov, Amelin, Zueva, Babayan
# License: BSD
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import array, dot, linalg, arange
from numpy import polyfit
from pandas import read_csv
from scipy.linalg import solve
from scipy import interpolate
from math import log, exp


def csv(path, delimiter = ','):
    """
    path = Путь к файлу
    delimiter = разделитель в файлe
    Возвращает два списка x, y 
    """
    points = read_csv(path, sep = delimiter)
    return [i for i in points[list(points)[0]]], [i for i in points[list(points)[1]]]

def factor(n):
    """Вычисление факториала"""
    r = 1
    for i in range(1,n+1):
        r = r*i   
    return r

def check_step(x):
    """Проверка шага на интервале"""
    s = x[1] - x[0]
    return s

def sol(x, y):
    """
    x, y = точки для решения многочлена
    Выдает коэффициенты многочлна Лагранжа 
    """
    
    main_mtrx = []
    for i in x:
        main_mtrx.append([i**(el) for el in range(len(x))])
    main_mtrx = array(main_mtrx)
    y = [[i] for i in y]
    y = array(y)
    answ = solve(main_mtrx, y)
    return answ    

def answer(x, y, answ, flag):
    """
    x, y = точки для апроксимации
    answ = Коэффициенты многочлена
    flag = вывести или выдать интерполяционную функцию 
    """
    h = ''
    answ = [float(el[0]) for el in answ]
    
    for ind, el in enumerate(answ[::-1]):
        h += f'+{el}*(x**{len(answ) - ind-1})'
    x1 = []
    y1 = []
    for i in range(int(min(x))*10, int(max(x))*10 + 1):
        x1.append(i/10)
        y1.append(eval(h.replace('x', str(i/10))))
    if flag:
        print('Многочлен: ' + h)
        print()
        for ind, el in enumerate(x):
            print('x: '+ str(el)+' '+'y: '+ str(y[ind])+ ' ' + 'fi: ' + str(round(eval(h.replace('x', str(el))))))
        
        plt.plot(x1,y1,color = 'black', label = 'interpol')
        plt.scatter(x, y, color = 'orange', label = 'dots')
        plt.legend()
    else:
        return [x1, y1] 
    
def lagrange(path, header = False, deli = ',', flag = True):
    """
    Объединение предыдущих двух методов
    """
    x, y = csv(path = path, delimiter = deli)
    return answer(x, y, sol(x, y), flag = flag)
    
def table(x, y):
    """
    x, y = точки для интерполяции
    Вывод, Таблица, для метода ньютона
    """
    table1 = [[0 for j in range(len(x))] for i in range(len(x)+1)]
    table1[0], table1[1] = x, y 
    for ind, el in enumerate(table1[2:]):
        for i in range(len(el)-ind-1):
            table1[ind+2][i] = -(table1[ind+1][i] - table1[ind+1][i+1])
    table1 = np.array(table1).T
    return table1

def newton_F(x, y, table):
    """
    x, y = Точки
    table = Таблица Ньютона
    Выход - Интерполирующая функция Ньютоном вперед"""
    
    P = ''
    h = check_step(x)
    
    for i in range(1, len(y) + 1):
        skobochki = ''
        for j in range(1, i):
            skobochki += f'*(x - {table[j-1][0]})'
            
        helpi = f'{table[0][i]}/({factor(i-1)}*{h**(i-1)})'
        if i == 1:
            P += f'({helpi}){skobochki}'
        else:
            P += f'+({helpi}){skobochki}'
    return P

def newton_S(x, y, table):
    """
    x, y = Точки
    table = Таблица Ньютона
    Выход - Интерполирующая функция Ньютоном назад"""
    
    P = ''
    h = check_step(x)
    
    for i in range(1, len(y)+1):
        skobochki = ''
        for j in range(1, i):
            skobochki += f'*(x - {table[len(y)-j][0]})'
        
        helpi = f'{table[len(y) - i][i]}/({factor(i-1)}*{h**(i-1)})'
        if i == 1:
            P += f'({helpi}){skobochki}'
        else:
            P += f'+({helpi}){skobochki}'

    return P
                  
def newton(path, deli = ',', FoS = 'F', flag = True):
    """
    Объединение всех Ньютонов
    path = путь к точкам 
    delimiter = разделитель в файле
    FoS = Путь метода
    flag = вывести или выдать интерполяционную функцию """
    x, y = csv(path = path, delimiter = deli)
    if FoS == 'F':
        F = newton_F(x, y, table(x, y))
    if FoS == 'S':
        F = newton_S(x, y, table(x, y))
    x1 = []
    y1 = []
    for i in range(int(min(x))*10, int(max(x))*10+1):
        x1.append(i/10)
        y1.append(eval(F.replace("x", str(i/10))))

    if flag:
        for ind, el in enumerate(x):
            print('x: '+ str(el)+' '+'y: '+ str(y[ind])+ ' ' + 'fi: ' + str(eval(F.replace("x", str(el)))))
        print(F)
        
        plt.plot(x1,y1,color = 'black', label = 'interpol')
        plt.scatter(x, y, color = 'orange', label = 'dots')
        plt.legend()
    else:
        return [x1, y1]
    
def linear(path, deli = ',', flag = True):
    """
    path = путь к точкам 
    delimiter = разделитель в файле
    flag = вывести или выдать интерполяционную функцию
    """
    x, y = csv(path, delimiter=deli)
    kx = solve([[sum([i**2 for i in x]),sum(x)],[sum(x), len(x)]],[[sum([y[i]*x[i] for i in range(len(x))])],[sum(y)]])
    x1 = []
    y1 = []
    for i in range(int(min(x)), int(max(x))+1):
        x1.append(i)
        y1.append(kx[0][0]*i+kx[1][0])
    if flag: 
        print(f'y = {kx[0][0]}*x+{kx[1][0]}')
        print()
        for ind, el in enumerate(x):
            print(f'x:{el} y:{y[ind]} fi:{round(kx[0][0]*el+kx[1][0])}')
        plt.plot(x1,y1,color = 'black', label = 'approx')
        plt.scatter(x, y, color = 'orange', label = 'dots')
        plt.legend()
    else:
        return [x1, y1]

def quadr(path, deli = ',', flag = True):
    """
    path = путь к точкам 
    delimiter = разделитель в файле
    flag = вывести или выдать интерполяционную функцию
    """    

    x, y = csv(path, delimiter=deli)
    kx = solve([[sum([i**4 for i in x]),sum([i**3 for i in x]),sum([i**2 for i in x])],
                [sum([i**3 for i in x]),sum([i**2 for i in x]),sum([i**1 for i in x])],
                [sum([i**2 for i in x]),sum([i**1 for i in x]),len(x)]],
                [[sum([(x[i]**2)*y[i] for i in range(len(x))])],
                 [sum([x[i]*y[i] for i in range(len(x))])],
                 [sum([y[i] for i in range(len(x))])]])
    x1 = []
    y1 = []
    for i in range(int(min(x))*10, (int(max(x))+1)*10):
        x1.append(i/10)
        y1.append(kx[0][0]*(i/10)**2+kx[1][0]*i/10+kx[2][0])
    if flag:
        print(f'y = {kx[0][0]}x**2+{kx[1][0]}x + {kx[2][0]}')
        print()
        for ind, el in enumerate(x):
            print(f'x:{el} y:{y[ind]} fi:{kx[0][0]*(el)**2+kx[1][0]*el+kx[2][0]}')
        plt.plot(x1,y1,color = 'black', label = 'approx')
        plt.scatter(x, y, color = 'orange', label = 'dots')
        plt.legend()
    else:
        return [x1, y1]
    
def norm(path, deli = ',', flag = True):
    x, y = csv(path, delimiter=deli)
    X = []
    Y = []
    for ind, el in enumerate(x):
        X.append([1, el, el**2])
        Y.append([log(y[ind])])
    B = dot(dot(linalg.inv(dot(array(X).T,array(X))), array(X).T), array(Y))
    a = exp(B[0][0] - (B[1][0]**2)/(4*B[2][0]))
    b = - (1/(B[2][0]))
    c = - (B[1][0]/(B[2][0]*2))
    f = lambda x: a * exp(-((x - c)**2/b))
    x1 = []
    y1 = []
    for el in arange(min(x), max(x), 0.1):
        x1.append(el)
        y1.append(f(el))
    if flag:
        print(f'y = {a}*e**(-((x - {c})**2/{b}')
        print()
        for ind, el in enumerate(x):
            print(f'x:{el} y:{y[ind]} fi:{f(el)}')
        plt.plot(x1,y1,color = 'black', label = 'approx')
        plt.scatter(x, y, color = 'orange', label = 'dots')
        plt.legend()
    else:
        return [x1, y1]
    
def numpy_poli(path, flag, deli = ','):
    """Вспомогательная функция для построения графиков сравнения"""
    x, y = csv(path = path, delimiter=deli)
    x1 = []
    y1 = []
    formula = ''
    m = polyfit(x, y, 2)
    for ind, el in enumerate(m):
        formula += f'+{el}*x**{len(m)-ind-1}'
    for x_el in range((int(min(x)))*10, (int(max(x))+1)*10):
        x1.append(x_el/10)
        y1.append(eval(formula.replace("x", str(eval(f'{x_el}/10')))))
    return [x1, y1]
    
def scipy_interpol(path, flag,  deli = ','):
    """Вспомогательная функция для построения графиков сравнения"""
    x, y = csv(path = path, delimiter=deli)
    f = interpolate.interp1d(x, y, kind = 'cubic')
    xnew = arange(min(x), max(x), 0.1)
    
    ynew = f(xnew)
    return [xnew, ynew]

def comparison(path, deli = ','):
    """Построение графиков сравнения"""
    x, y = csv(path = path, delimiter = deli)
    fig, axes = plt.subplots(2, 4)
    
    function_list = [[lagrange, scipy_interpol], [linear,quadr , numpy_poli]]
    for ind, i_a in enumerate(function_list):
        for i, func in enumerate(i_a):
            m = func(path = path, flag = False, deli = deli)
            axes[ind][i].set_title(func.__name__)
            axes[ind][i].plot(m[0], m[1],color = 'black', label = 'approx')
            axes[ind][i].scatter(x, y, color = 'orange', label = 'dots')
            axes[ind][i].legend()

    m = newton(path, FoS = 'F', flag =False)
    axes[0][2].set_title('Newton F')
    axes[0][2].plot(m[0], m[1],color = 'black', label = 'approx')
    axes[0][2].scatter(x, y, color = 'orange', label = 'dots')
    axes[0][2].legend()

    m = newton(path, FoS = 'S', flag =False)
    axes[0][3].set_title('Newton S')
    axes[0][3].plot(m[0], m[1],color = 'black', label = 'approx')
    axes[0][3].scatter(x, y, color = 'orange', label = 'dots')
    axes[0][3].legend()

    m = norm(path,deli = deli,  flag = False)
    axes[1][3].set_title('Нормальное')
    axes[1][3].plot(m[0], m[1],color = 'black', label = 'approx')
    axes[1][3].scatter(x, y, color = 'orange', label = 'dots')
    axes[1][3].legend()
    
    fig.set_figheight(10)
    fig.set_figwidth(16)
    
