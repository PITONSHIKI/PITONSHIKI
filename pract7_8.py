"""
Пакет pract7_8 решает задачи практической 7,8.

Решение задачи коммивояжера. Алгоритмы оптимизации.
FUNCTIONS:
    main() - основная работа программы, сравнение методов
    get_matrix() - ввод с клавиатуры матрицы
    randomize_matrix() - генерация случайной матрицы
    get_matrix_from_csv() - матрица из .csv файла
    branch_and_bound() - решение методом ветвей и границ
    simulated_annealing() - решение методом отжига
    ACO() - решение методом муравьиной колонии
    greedy() - решение жадным алгоритмом

Latipov, Amelin, Zueva, Babayan
# License: BSD
"""

import copy
import random
import csv
import math
import time
import matplotlib.pyplot as plt

INF = float('inf')

def get_matrix():
    '''Ввод матрицы с клавиатуры'''
    def print_matrix(matrix):
        """Форматированный вывод матрицы"""
        for string in matrix:
            print('|',end='')
            for el in string:
                print('{:^5}'.format(el),end='|')
            print()
        
    print('''Ввод данных производится в виде матрицы расстояний.
Для 1 и 2 городов она не имеет смысла.
Минимальная размерность матрицы - 3 (3 пункта назначения)''')
    while True:
        length=input('Введите размерность: ')
        try:
            length=int(length)
            if length>=3:
                break
            else:
                print('Некорректный ввод. Введите число большее или равное 3') 
        except ValueError:
            print('Некорректный ввод. Введите число большее или равное 3')
    matrix=[list(el)for el in [[None]*length]*length]
    for i in range(length):
        for j in range(length):
            if i==j:
                matrix[i][j]=INF
            else:
                while True:
                    el=input(f'Введите неотрицательное \
число(вещественные числа вводятся через точку)\nЭлемент {i+1},{j+1}: ')
                    try:
                        el=float(el)
                        if el>=0:
                            break
                        else:
                            print('Некорректный ввод. Ознакомтесь с условиями еще раз')
                    except ValueError:
                        print('Некорректный ввод. Ознакомтесь с условиями еще раз')
                matrix[i][j]=el
    print('--------------------\nВведенная матрица расстояний:')            
    print_matrix(matrix)
    print('--------------------')
    return matrix

def randomize_matrix(rank = 5, min_ = 10, max_ = 99):
    """Генерация случайной матрицы стоимостей"""
    INF = float('inf')
    random_matrix = [[random.randint(min_, max_) if i != j else INF for i in range(rank)] for j in range(rank)]
    for string in random_matrix:
        print('|',end='')
        for el in string:
            print('{:^5}'.format(el),end='|')
        print()
    return random_matrix

def get_matrix_from_csv(filename = 'matrix', modify = int):
    """Чтение из csv файла"""
    matrix = []
    col = 0
    with open(f'{filename}.csv', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=';')
        for row in spamreader:
            if col:
                if len(row) != col:
                    print("Read csv error")
                    return None
            elif len(row)>1:
                col = len(row)
            else:
                print("Read matrix error")
                return None
            matrix.append(row)
            for item in row:
                if "j" in item:
                    modify = complex
                    
    if len(matrix)!=0 and (len(matrix) == len(matrix[0])):
        try:
            for i in range(len(matrix)):
                for j in range(len(matrix[i])):
                    if i!=j:
                        matrix[i][j] = modify(matrix[i][j])
                    else:
                        matrix[i][j] = float('inf')
        except:
            print("Read type error")
            return None 
    else:
        print("rank error")
        return None
    for string in matrix:
        print('|',end='')
        for el in string:
            print('{:^5}'.format(el),end='|')
        print()
    return matrix

def branch_and_bound(matrix):
    """Решение методом ветвей и границ"""
    INF = float('inf')
    def string_reduce(matrix):
        """редуция по строкам"""
        for i,el in enumerate(matrix):
            min_=min(el)
            if min_ == INF:
                min_ = 0
            matrix[i].append(min_)
            matrix[i]=[obj-matrix[i][-1] for obj in matrix[i][:-1]]+matrix[i][-1:]   
        return matrix

    def column_reduce(matrix):
        """редуция по столбцам"""
        matrix.append([])
        for i in range(len(matrix)-1):
            temp=[]
            for el in matrix[:-1]:
                temp.append(el[i])
            min_=min(temp)
            if min_ == INF:
                min_ = 0
            matrix[-1].append(min_)
            for j in range(len(matrix[i])-1):
                matrix[j][i]-=matrix[-1][i] 
        return matrix

    def lower(matrix):
        """находим нижнюю границу приведения"""
        lower=0
        for i in range(len(matrix)-1):
            lower+=(matrix[i][-1]+matrix[-1][i])
        return lower

    def reduce_border(matrix):
        matrix.pop()
        for i in range(len(matrix)):
            matrix[i].pop()
        return matrix

    def temp_matrix_create(matrix,i,j,type_=0):
        temp_matrix= copy.deepcopy(matrix)
        if type_:
            for l in range(len(temp_matrix)):
                temp_matrix[i][l]= INF
                temp_matrix[l][j]= INF
            temp_matrix[j][i]=INF
            temp_matrix=column_reduce(string_reduce(temp_matrix))
        else:
            temp_matrix[i][j]= INF
            temp_matrix=column_reduce(string_reduce(temp_matrix))
        return temp_matrix

    def max_const_summ(matrix):
        max_element=[-1,-1,-1]
        for i,string in enumerate(matrix):
            for j,el in enumerate(string):
                if el==0:
                    temp_matrix=temp_matrix_create(matrix,i,j) 
                    element_summ=lower(temp_matrix)
                    if element_summ>max_element[0]:
                        max_element=[element_summ,i,j]

        return max_element

    def max_edge_summ(matrix,i,j):
        matrix=temp_matrix_create(matrix,i,j,1)
        return lower(matrix)

    def step_finder(matrix,steps):
        asterisk_summ=max_const_summ(matrix)
        no_asterisk_summ=max_edge_summ(matrix,asterisk_summ[1],asterisk_summ[2])
        if no_asterisk_summ<=asterisk_summ[0]:
            flag=1
        else:
            flag=0
        step = asterisk_summ[1:]
        start = step[1]
        temp_start = step[0]
        count_rep = 0
        for i in range(len(matrix)):
            for x in steps:
                if x[0] == start:
                    start = x[1]
                    count_rep += 1
                    break
            if (temp_start == start) and (count_rep != 0):
                matrix[step[0]][step[1]] = INF
                return step_finder(matrix,steps)
        matrix=reduce_border(temp_matrix_create(matrix,asterisk_summ[1],asterisk_summ[2],flag))
        return matrix,asterisk_summ[1:],min(no_asterisk_summ,asterisk_summ[0])

    def end_finder(matrix,steps):
        return len(matrix)-2==len(steps)

    def end_step_finder(matrix,steps):
        temp_start,temp_end,temp_steps,temp_total=[],[],[],0
        for el in steps:
            temp_start.append(el[0])
            temp_end.append(el[1])
        temp_start={i for i in range(len(matrix))}-set(temp_start)
        temp_end={i for i in range(len(matrix))}-set(temp_end)
        temp_start=list(temp_start)
        temp_end=list(temp_end)
        finder=temp_end[0]
        while True:
            for el in steps:
                if el[0]==finder:
                    finder=el[1]
                    break
            else:
                break    
        temp_steps.append([finder,temp_end[1]])
        temp_total+=matrix[finder][temp_end[1]]
        temp_start=list(set(temp_start)-{finder})[0]
        temp_steps.append([temp_start,temp_end[0]])
        temp_total+=matrix[temp_start][temp_end[0]]
        return temp_steps,temp_total

    def steps_sort(steps):
        string_steps=[]
        started=steps[0][0]
        while steps:
            for i,el in enumerate(steps):
                if el[0]==started:
                    started=el[1]
                    string_steps.append([g+1 for g in steps.pop(i)])
        return string_steps

    def get_way(steps):
        sorted_steps = steps_sort(steps)
        way = [item[0] for item in sorted_steps]
        way.append(way[0])
        way = [str(i) for i in way]
        way = "-".join(way)
        return way
    
    matrix = copy.deepcopy(matrix)
    steps=[]
    total=0
    matrix=column_reduce(string_reduce(matrix))
    summ=lower(matrix)
    matrix=reduce_border(matrix)
    total+=summ
    iteration = 0
    while not end_finder(matrix,steps):
        matrix,step,summ=step_finder(matrix, steps)
        iteration += 1
        steps.append(step)
        total+=summ
    step,summ=end_step_finder(matrix,steps)
    steps.extend(step)
    total+=summ
    way = get_way(steps)
    return way, total, iteration

def simulated_annealing(
        matrix,
        printed = False,
        initial_T = 10,
        end_T = 0.0000000000001,
        precision = 10000
        ):
    """Решение методом имитации отжига"""
    matrix = copy.deepcopy(matrix)
    RANK = len(matrix)
    def way_len(matrix, way):
        way_length = 0
        for i in range(len(matrix)-1):
            way_length += matrix[way[i]][way[i+1]]
        way_length += matrix[way[-1]][way[0]]
        return way_length
    
    way = list(range(RANK))
    random.shuffle(way)
    T = initial_T
    current_way = way_len(matrix, way)
    i = 0
    iteration = 0
    container = []
    for i in range(precision):
        iteration += 1
        i += 1
        point1 = random.randint(0, RANK-1)
        while True:
            point2 = random.randint(0, RANK-1)
            if point1 != point2:
                break  
        candidate = copy.deepcopy(way)
        candidate[point1], candidate[point2] = candidate[point2], candidate[point1]
        candidate_way = way_len(matrix, candidate)
        if (candidate_way < current_way):
            current_way = candidate_way
            way = candidate
        else:
            p = math.exp(-(candidate_way - current_way) / T)
            if not(0<= p <= 1):
                raise ValueError 
            rand_point = random.random()
            if (rand_point <= p):
                current_way = candidate_way
                way = candidate
        T = initial_T  / math.log(1+i)
        if printed:
            print(way, current_way)
        container.append(current_way)
        if (T < end_T):
            break     
    way.append(way[0])
    way = [str(i+1) for i in way]
    way = '-'.join(way)
    return way, current_way, iteration, container

def ACO(matrix,
        printed = False,
        alpha = 0.6,
        betta = 0.65,
        ro = 0.3,
        start_pheromone = 0.1,
        precision = 999):
    """Решение алгоритмом муравьиной колонии"""
    RANK = len(matrix)
    pheromone = [[start_pheromone if i!=j else 0 for i in range(RANK)] for j in range(RANK)]
    matrix = copy.deepcopy(matrix)
    
    def way_len(matrix, way):
        way_length = 0
        for i in range(len(matrix)):
            way_length += matrix[way[i]][way[i+1]]
        return way_length
    
    def greedy(matrix, first):
        RANK = len(matrix)
        way = []
        matrix_copy = copy.deepcopy(matrix)
        way.append(first)
        for i in range(RANK-1):
            for j in range(RANK):
                matrix_copy[j][first] = INF
            ind, min_value = min(enumerate(matrix_copy[first]), key=lambda i_v: i_v[1])
            matrix_copy[first][ind] = matrix_copy[ind][first] = INF
            first = ind
            way.append(ind)
        way.append(way[0])
        way_length = way_len(matrix, way)
        return way, way_length
    
    def probability(pheromone, matrix, j_points, start_point, target_point):
        summ = 0
        for next_point in j_points:
            summ += (pheromone[start_point][next_point]**alpha)*((1/matrix[start_point][next_point])**betta)
        target = (pheromone[target_point][next_point]**alpha)*((1/matrix[target_point][next_point])**betta)
        return target/summ
    
    def find_way_by_pheromone(pheromone, matrix):
        first = 0
        way = [first]
        pheromone_temp = copy.deepcopy(pheromone)
        for i in range(RANK-1):
            for j in range(RANK):
                pheromone_temp[j][first] = 0
            ind, max_value = max(enumerate(pheromone_temp[first]), key=lambda i_v: i_v[1])
            pheromone_temp[first][ind] = pheromone_temp[ind][first] = 0
            first = ind
            way.append(ind)
        way.append(way[0])
        way_length = way_len(matrix, way)
        way = [str(i+1) for i in way]
        way = '-'.join(way)
        return way, way_length 
    
    iteration = 0
    container = []
    for _ in range(precision):
        iteration += 1
        pheromone_temp = [[0 for _ in range(RANK)] for _ in range(RANK)]
        for start in range(RANK):
            way = [start]
            j_points = list(range(RANK))
            j_points.pop(way[0])
            for point in range(len(j_points)):

                probabilities =[]
                for next_point in j_points:
                    proby = probability(pheromone, matrix, j_points, way[point], next_point)
                    probabilities.append(proby)
                probabilities = [sum(probabilities[:item+1]) for item in range(len(probabilities))]
                random_point = random.random()
                for i,prob in enumerate(probabilities):
                    if random_point< prob:
                        break
                way.append(j_points[i])
                j_points.remove(way[-1])
            way.append(way[0])
            l_min = greedy(matrix, way[0])[1]
            l_current = way_len(matrix, way)
            for i in range(len(way)-1):
                pheromone_temp[way[i]][way[i+1]] += l_min/l_current
        for i in range(RANK):
            for j in range(RANK):
                pheromone[i][j] = pheromone[i][j]*(1-ro) + pheromone_temp[i][j]
        results = find_way_by_pheromone(pheromone, matrix)
        container.append(results[1])        
        if printed:
            print(results)
    return *find_way_by_pheromone(pheromone, matrix), iteration, container

def greedy(matrix):
    """Решение жадным алгоритмом"""
    RANK = len(matrix)
    def way_len(matrix, way):
        way_length = 0
        for i in range(len(matrix)):
            way_length += matrix[way[i]][way[i+1]]
        return way_length
    way = []
    matrix_copy = copy.deepcopy(matrix)
    first = random.randint(0,RANK-1)
    way.append(first)
    for i in range(RANK-1):
        for j in range(RANK):
            matrix_copy[j][first] = INF
        ind, min_value = min(enumerate(matrix_copy[first]), key=lambda i_v: i_v[1])
        matrix_copy[first][ind] = matrix_copy[ind][first] = INF
        first = ind
        way.append(ind)
    way.append(way[0])
    way_length = way_len(matrix, way)
    way = [str(i+1) for i in way]
    way = '-'.join(way)
    return way, way_length

def main(matrix = None):
    """Решение основными методами и их сранение"""
    if matrix is None:
        matrix = randomize_matrix()
    start_time = time.time()
    result = branch_and_bound(matrix)
    end_time = time.time()
    print(f"""Алгоритм: Ветвей и границ
Количество итераций: {result[2]}
Время выполнения: {end_time-start_time}
Результат: {result[0]}
Стоимость: {result[1]}""")
    start_time = time.time()
    result = simulated_annealing(matrix, precision = 10000)
    end_time = time.time()
    print(f"""Алгоритм: имитации отжига
Количество итераций: {result[2]}
Время выполнения: {end_time-start_time}
Результат: {result[0]}
Стоимость: {result[1]}""")
    plt.figure(figsize=(20, 6))
    plt.plot(list(range(len(result[3]))), result[3])
    start_time = time.time()
    result = ACO(matrix, precision = 5000)
    end_time = time.time()
    print(f"""Алгоритм: Муравьиной колонии
Количество итераций: {result[2]}
Время выполнения: {end_time-start_time}
Результат: {result[0]}
Стоимость: {result[1]}""")
    plt.figure(figsize=(20, 6))
    plt.plot(list(range(len(result[3]))), result[3])