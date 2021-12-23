"""
Пакет pract1 решает задачи практической 1.

FUNCTIONS:
    main() - основная работа программы

Latipov, Amelin, Zueva, Babayan
# License: BSD
"""

import math
import numexpr as ne
from datetime import datetime 
import plotly
import plotly.graph_objs as go
import plotly.express as px
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
from scipy import integrate
from scipy.misc import derivative


def main():
    """
    Основная функция сравнения вычисений.
    
    Определяет интеграл, на основе уже имеющейся функции.
    Определяет производную на основе уже имеющейся функции.
    Определеяет производную, посчитанную нашей фукнцией.
    Определяет интеграл, посчитанный нашей функцией
    """

    #Определенный интеграл, вычисленный по формуле трапеций
    def integr(f, l, r, step):
        """Возвращает определенный интеграл, вычисленный по формуле трапеций"""
        ans = 0
        while(l +step< r):
            k1 = f(l)
            k2 = f(l + step)
            if(k1 == np.nan):
                return np.nan
            if(k2 == np.nan):
                return np.nan
            ans += (k1 + k2) / 2 * step
            l += step
        k1 = f(l)
        k2 = f(r)
        if(k1 == np.nan):
            return np.nan
        if(k2 == np.nan):
            return np.nan
        ans += (k1 + k2) / 2 * step
        l += step
        return ans

    #Производная в точке, вычисленная по определению
    def diff(f, x, dx):
        """Производная в точке по определению"""
        k1 = f(x)
        k2 = f(x + dx)
        if(k1 == np.nan):
            return np.nan
        if(k2 == np.nan):
            return np.nan
        return (k2 - k1)/dx

    s = input("f(x) = ")
    #Корректировка ввода
    while(s.find("ln") != -1):
        s = s[0 : s.find("ln")] + "log" + s[s.find("ln") + 2 : len(s)]
    while(s.find("tg") != -1):
        s = s[0 : s.find("tg")] + "tan" + s[s.find("tg") + 2 : len(s)]
    while(s.find("^") != -1):
        s = s[0 : s.find("^")] + "**" + s[s.find("^") + 1 : len(s)]
    def g(x):
        return ne.evaluate(s)
    l = float(input("Левая граница: "))
    r = float(input("Правая граница: "))
    n = float(input("Точность: "))
    print("Количество трапеций: ", (r - l) / n)
    x = []
    y1 = [] #Список значений для определенного интеграла, посчитанного нами
    y2 = [] #Список значений для нашей функции
    y3 = [] #Список значений для производной в точке, посчитанной нами
    y4 = [] #Список значений для определенного интеграла, посчитанного системой
    y5 = [] #Список значений для производной в точке, посчитанной системой
    #Пронумерованы также, но без тройки
    test1 = []
    test2 = []
    test3 = []
    test4 = []
    #Рассчитаем значения
    while(l <= r):
        #Добавляем новый икс
        x.append(l)
        #Рассчитываем интеграл
        st = datetime.now() #Время начала расчёта
        if(x[-1] < 0):
            y1.append(integr(g, x[-1], 0.0, n))
        else:
            y1.append(integr(g, 0.0, x[-1], n))    
        test1.append(datetime.now() - st) #Фиксация времени завершения расчёта
        #Рассчитываем значение функции
        y2.append(g(x[-1]))
        #Рассчитываем значение производной
        st = datetime.now() #Время начала расчёта
        y3.append(diff(g, x[-1], n))
        test2.append(datetime.now() - st) #Фиксация времени завершения расчёта
        #Рассчитываем значение интеграла по библиотеке
        st = datetime.now() #Время начала расчёта
        if(x[-1] < 0):
            v, err = integrate.quad(g, x[-1], 0.0)
        else:
            v, err = integrate.quad(g, 0.0, x[-1])
        y4.append(v)
        test3.append(datetime.now() - st) #Фиксация времени завершения расчёта
        #Расчитываем значение производной по библиотеке
        st = datetime.now() #Время начала расчёта
        y5.append(derivative(g, x[-1], n))
        test4.append(datetime.now() - st) #Фиксация времени завершения расчёта
        l += n
    data=pd.DataFrame({
        'x' : list(map(('{:.' + str(len(str(n))) + 'f}').format, x)), 
        'f(x)' : list(map('{:.6f}'.format, y2)),
        'f\'(x)' : list(map('{:.6f}'.format, y3)),
        'F(x)' : list(map('{:.6f}'.format, y1))
        }) #создание датафрейма
    print(data)
    data1=pd.DataFrame({
        'F(x)' : test1,
        'F(x) система' : test3
        }) #создание датафрейма
    print(data1)
    data2=pd.DataFrame({
        'f\'(x)' : test2,
        'f\'(x) система' : test4
        })
    print(data2)
    #Нанесем графики, заданные списками аргументов и значений функции, на фигуру и подпишем их
    fig = go.Figure()
    fig.add_trace(go.Scatter(x = x, y = y1, name = "$$F(x)$$"))
    fig.add_trace(go.Scatter(x = x, y = y2, name = "$$f(x) = " + s + "$$"))
    fig.add_trace(go.Scatter(x = x, y = y3, name = "$$f'(x)$$"))
    fig.add_trace(go.Scatter(x = x, y = y4, name = "$$F(x) - система$$"))
    fig.add_trace(go.Scatter(x = x, y = y5, name = "$$f'(x) - система$$"))
    fig.show()#Покажем графики