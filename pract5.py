"""
Пакет pract5 решает задачи практической 5.

FUNCTIONS:
    func_ahc(csv_file) - получение АЧХ функции

    wavelets_4lvl_D(csv_file, name) - получение вейвлета
        csv_file - путь к файлу с массивом точек
        name - имя вида вейвлета

    wavelets_4lvl_I(csv_file, name) - получение вейвлета (другое представление)
        csv_file - путь к файлу с массивом точек
        name - имя вида вейвлета

    all_wavelets(csv_file) - Сравнение всех вейвлетов

Latipov, Amelin, Zueva, Babayan
# License: BSD
"""

import pywt
import matplotlib.pyplot as plt
from scipy.fft import rfft, rfftfreq, irfft
from pandas import read_csv
from numpy import std, arange, var

def func_ahc(path):
    """Получение АЧХ"""
    X, Y = list(read_csv(path, names = [0,1])[0]), list(read_csv(path, names = [0,1])[1])
    plt.plot(X, Y, color = 'black')
    plt.title('Исходный')
    plt.show()
    period = (len(X))/(X[-1] - X[0])
    X = rfftfreq(len(Y), 1 / period) 
    Y = [abs(el) for el in rfft(Y)] 
    step = X[1] - X[0]
    S0 = ((Y[Y.index(max(Y)) + 1] + Y[Y.index(max(Y)) + 1])/2 + Y[Y.index(max(Y))])*step
    Si = 0
    for i in range(1, len(Y)):
        Si += (Y[i-1] + Y[i])* step/2
    print(f'Коэффициент несинусоидальности = {(Si-S0)/S0}')
    y = irfft(Y)
    x = list(read_csv(path, names = [0,1])[0])
    plt.plot(X, Y, color = 'black')
    plt.grid()
    plt.xlabel('Гц')
    plt.title('АХЧ')
    plt.show()
    plt.plot(x, y, color = 'black')
    plt.title('Полученный')

def wavelets_4lvl_D(path, name):
    """Вейвлеты четырех уровней"""
    X, Y = list(read_csv(path, names = [0,1])[0]), list(read_csv(path, names = [0,1])[1])
    coeffs = pywt.wavedec(Y, name, level=4)
    fig, axes = plt.subplots(3, 2)
    
    axes[0][0].plot(X, Y, color = 'black')
    axes[0][0].set(title = 'Исходные данные')
    axes[0][0].grid()
    
    axes[0][1].plot(coeffs[0], color = 'black')
    axes[0][1].set(title = 'Апроксимация')
    axes[0][1].grid()
        
    axes[1][0].plot(coeffs[1], color = 'black')
    axes[1][0].set(title = 'Уровень 4')
    axes[1][0].grid()
    
    axes[1][1].plot(coeffs[2], color = 'black')
    axes[1][1].set(title = 'Уровень 3')
    axes[1][1].grid()

    axes[2][0].plot(coeffs[3], color = 'black')
    axes[2][0].set(title = 'Уровень 2')
    axes[2][0].grid()
    
    axes[2][1].plot(coeffs[4], color = 'black')
    axes[2][1].set(title = 'Уровень 1')
    axes[2][1].grid()

    fig.set_figwidth(16)
    fig.set_figheight(16)
    
    summ = 0
    for i in range(1, len(coeffs) - 1):
        summ += std(coeffs[i])**2
    summ = summ/std(coeffs[0])
    
    print(f'Коэффициент несинусиоидальности {summ}')

def wavelets_4lvl_I(path, name):
    """Вейвлеты четырех уровней другой вид"""
    X, Y = list(read_csv(path, names = [0,1])[0]), list(read_csv(path, names = [0,1])[1])
    fig, axes = plt.subplots(2,2)
    fig.suptitle(f'{name}')
    wavelet = pywt.ContinuousWavelet(name, level = 4)
    coef, freqs = pywt.cwt(Y, arange(1, 257), wavelet)
    axes[0][0].plot(coef[0])
    axes[0][0].set(title= 'level = 1')
    coef, freqs = pywt.cwt(coef[0], arange(1, 257), wavelet)
    axes[0][1].plot(coef[0])
    axes[0][1].set(title= 'level = 2')
    coef, freqs = pywt.cwt(coef[0], arange(1, 257), wavelet)
    axes[1][0].plot(coef[0])
    axes[1][0].set(title= 'level = 3')
    coef, freqs = pywt.cwt(coef[0], arange(1, 257), wavelet)
    axes[1][1].plot(coef[0])
    axes[1][1].set(title= 'level = 4')
    fig.set_figwidth(16)
    fig.set_figheight(6)
    
def all_wavelets(path):
    """Сравнение четырех вейвлетов (hhar, db2, mexh, gaus1)"""
    X, Y = list(read_csv(path, names = [0,1])[0]), list(read_csv(path, names = [0,1])[1])
    coeffs_haar = pywt.wavedec(Y, 'haar', level=4)
    coeffs_db2 = pywt.wavedec(Y, 'db2', level=4)
    wavelet_mexh = pywt.ContinuousWavelet('mexh', level = 4)
    wavelet_gaus1 = pywt.ContinuousWavelet('gaus1', level = 4)
    coef_mexh, freqs = pywt.cwt(Y, arange(1, 257), wavelet_mexh)
    coef_gaus1, freqs = pywt.cwt(Y, arange(1, 257), wavelet_gaus1)
    fig, axes = plt.subplots(3,2)
    
    axes[0][0].plot(X, Y)
    axes[0][0].set(title = 'Исходный')
    
    axes[0][1].plot(coeffs_haar[0])
    axes[0][1].set(title = f'Хаар, Дисперсия = {var(coeffs_haar[0])}')
    
    axes[1][0].plot(coeffs_db2[0])
    axes[1][0].set(title = f'ДБ2, Дисперсия = {var(coeffs_db2[0])}')
    
    axes[1][1].plot(coef_mexh[0])
    axes[1][1].set(title = f'Мексик, Дисперсия = {var(coef_mexh[0])}')
    
    axes[2][0].plot(coef_gaus1[0])
    axes[2][0].set(title = f'Гаус, Дисперсия = {var(coef_gaus1[0])}')
    
    fig.set_figwidth(16)
    fig.set_figheight(16)   

    