import math
import random
import numpy as np

def powmod(a, p, N):
    '''Возведение в степень и сравнение по модулю, т.е. a^(p-1) эквивалентно 1 по модулю p'''
    pow2 = 1
    result = 1
    while (pow2 <= p):
        pow2 *= 2
    pow2 //= 2
    while (pow2 > 0):
        if ((p // pow2) % 2 == 1):
            result *= a
        if (pow2 > 1):
            result **= 2
        pow2 //= 2
        result %= N
    return result

def matrix_powmod(a, p, N):
    '''Возведение в степень матрицы по модулю'''
    pow2 = 1
    result = np.eye(a.shape[0], dtype=int)
    while (pow2 <= p):
        pow2 *= 2
    pow2 //= 2
    while (pow2 > 0):
        if ((p // pow2) % 2 == 1):
            result = np.dot(result, a)
        if (pow2 > 1):
            result = np.dot(result, result)
        pow2 //= 2
        result %= N
    return result


def jacobi(P, Q):
    '''Проверка на число Якоби'''
    P %= Q
    R = 0
    result = 1
    if (Q % 2 == 0):
        return 0
    while (P > 1):
        while (P % 2 == 0):
            if ((Q ** 2 - 1) % 16 != 0):
                result *= -1
            P //= 2
        R = Q % P
        if ((Q - 1) % 2 == 1 and (P - 1) % 2 == 1):
            result *= -1
        Q = P
        P = R
    return result

def is_prime(N):
    '''Определяем - простое число или нет'''

    prime = True
    # Так как простые числа натуральные и больше 1
    if (N < 2):
        return False

    # Ищем все делители до 10^6, если число меньше, то выполняем обычный алгоритм для поиска простых чисел
    for i in range(2, min([N, 1000000])):
        prime = prime and (N % i != 0)
        if not prime:
            return prime
    # Если ничего не найдено, а исходное число меньше 10^12, то оно простое (ищем выше до половины)
    if (N < 10 ** 12):
        return prime
    divisors_f1 = []
    f1 = 1
    r1 = N - 1
    # Ищем простые делители меньше 1 млн у числа N-1, 
    # затем для каждого простого делителя p случайно выбираем число m от 2 до N-2 такое,
    # чтобы m^(N-1) - 1 делилось на N, а m^((N-1)/p) - 1 не делилось
    for i in range(2, 1000000):
        while (r1 % i == 0):
            r1 //= i
            f1 *= i
            if (i != 2):
                divisors_f1.append(i)
    random.seed()
    ##
    a_example = random.randrange(2, N - 1)
    prime = prime and (powmod(a_example, N - 1, N) == 1) 
    if not prime:
        return prime
    suits_a = False
    exists_a = False
    is_primitive_root = True
    for p in divisors_f1:
        for i in range(7):
            a1_example = random.randrange(2, N - 1)
            suits_a = (powmod(a1_example, N-1, N) == 1) and (powmod(a1_example, (N - 1)//p, N) != 1)
            exists_a = exists_a or suits_a
        is_primitive_root = is_primitive_root and exists_a
        exists_a = False
    prime = prime and is_primitive_root
    for i in range(12):
        # берём неразложенный остаток R1 числа N-1 и проверяем для m^((N-1)/R1) - 1
        a1_example = random.randrange(2, N - 1)
        suits_a = (powmod(a1_example, N-1, N) == 1) and (powmod(a1_example, (N - 1)//r1, N) != 1)
        exists_a = exists_a or suits_a
    prime = prime and exists_a
    if not prime:
        return prime
    # Если у числа не делителей меньше, то очевидно, что нет смысла проверять числа меньше 10^21
    if (N < 10 ** 21):
    	return prime
    # Если больше 10^21, то идём дальше и проверяем f2 по последовательности Люка
    divisors_f2 = []
    f2 = 1
    r2 = N + 1
    for i in range(2, 10000000):
        while (r2 % i == 0):
            r2 //= i
            f2 *= i
            if (i != 2):
                divisors_f2.append(i)
    suits_a = False
    exists_a = False
    is_primitive_root = True
    iters = 0
    for q in divisors_f2:
        while (iters < 288 and not exists_a):
            # задаём последовательность x_n = ax_{n-1} + bx_{n-2}, предварительно проверив, чтобы дискриминант
            # многочлена x^2 - ax - b не был квадратичным вычетом (чтобы символ Якоби равнялся -1)
            a2_example = random.randrange(2, N - 1)
            b2_example = random.randrange(2, N - 1)
            if (jacobi(a2_example ** 2 + 4 * b2_example, N) == -1):
                # Задав начальные условия x_0 = 0, x_1 = 1, считаем x_{N+1}
                # и для каждого простого делителя N+1, меньшего 1 млн, считаем x_{(N+1)/q}.
                # чтобы x_{N+1} делился на N, а x_{(N+1)/q} не делился,
                U_example = np.array([[0, 1], [b2_example, a2_example]])
                vec = np.array([0, 1])
                vecres = np.dot(matrix_powmod(U_example, N+1, N), vec)
                vecresq = np.dot(matrix_powmod(U_example, (N+1)//q, N), vec)
                suits_a = (vecres[0] == 0) and (vecresq[0] % N != 0)
                exists_a = exists_a or suits_a
            iters += 1
        is_primitive_root = is_primitive_root and exists_a
        exists_a = False
        iters = 0
    prime = prime and is_primitive_root
    exists_a = False
    while (iters < 288 and not exists_a):
        a2_example = random.randrange(2, N - 1)
        b2_example = random.randrange(2, N - 1)
        if (jacobi(a2_example ** 2 + 4 * b2_example, N) == -1):
            # Аналогично предыдущему надо проверить для x_{(N+1)/R2}
            U_example = np.array([[0, 1], [b2_example, a2_example]])
            vec = np.array([0, 1])
            vecres = np.dot(matrix_powmod(U_example, N+1, N), vec)
            vecresq = np.dot(matrix_powmod(U_example, (N+1)//r2, N), vec)
            suits_a = (vecres[0] % N == 0 and vecresq[0] % N != 0)
            exists_a = exists_a or suits_a
        iters += 1
    iters = 0
    prime = prime and exists_a
    if not prime:
        return prime
    #дальше начинается проверка N^2 + 1, F4 и R4
    divisors_f4 = []
    f4 = 1
    r4 = N ** 2 + 1
    for i in range(2, 10000000):
        while (r4 % i == 0):
            r4 //= i
            f4 *= i
            if (i != 2):
                divisors_f4.append(i)
    suits_a = False
    exists_a = False
    is_primitive_root = True
    # сноска 1
    for q in divisors_f4:
        while (iters < 288 and not exists_a):
            c_c = random.randrange(2, N - 1)
            d_d = random.randrange(2, N - 1)
            h_h = random.randrange(2, N - 1)
            k_k = random.randrange(2, N - 1)
            if (jacobi(d_d, N) == -1 and jacobi(c_c ** 2 - 16 * d_d, N) == -1):
                p1_p1 = 4 * (h_h ** 2 * 2 + 2 * h_h * k_k * c_c + 2 * k_k * k_k * d_d)
                p2_p2 = (p1_p1 ** 2 - 16 * d_d) // 4
                q_q = (p1_p1 ** 2) // 16 + d_d - (h_h * h_h * c_c + k_k * k_k * c_c * d_d +
                                                        8 * h_h * k_k * d_d)
                # сноска 2
                U_example = np.array([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1],
                                      [-q_q ** 2, q_q * p1_p1, -p2_p2 - 2 * q_q, p1_p1]])
                vec_example = np.array([0, 1, p1_p1, p1_p1 ** 2 - p2_p2 - 3 * q_q])
                vecres = np.dot(matrix_powmod(U_example, N ** 2 + 1, N), vec_example)
                vecresq = np.dot(matrix_powmod(U_example, (N ** 2 + 1)//q, N), vec_example)
                suits_a = (vecres[0] % N == 0 and vecresq[0] % N != 0)
                exists_a = exists_a or suits_a
            iters += 1
        is_primitive_root = is_primitive_root and exists_a
        exists_a = False
        iters = 0

    exists_a = False
    prime = prime and is_primitive_root
    # сноска 3
    while (iters < 288 and not exists_a):
        c_c = random.randrange(2, N - 1)
        d_d = random.randrange(2, N - 1)
        h_h = random.randrange(2, N - 1)
        k_k = random.randrange(2, N - 1)
        if (jacobi(d_d, N) == -1 and jacobi(c_c ** 2 - 16 * d_d, N) == -1):
            p1_p1 = 4 * (h_h ** 2 * 2 + 2 * h_h * k_k * c_c + 2 * k_k * k_k * d_d)
            p2_p2 = (p1_p1 ** 2 - 16 * d_d) // 4
            q_q = (p1_p1 ** 2 + 16 * d_d - 16 * (h_h * h_h * c_c + k_k * k_k * c_c * d_d +
                                                8 * h_h * k_k * d_d)) // 16
            U_example = np.array([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1],
                                [-q_q ** 2, q_q * p1_p1, -p2_p2 - 2 * q_q, p1_p1]])
            vec_example = np.array([0, 1, p1_p1, p1_p1 ** 2 - p2_p2 - 3 * q_q])
            vecres = np.dot(matrix_powmod(U_example, N ** 2 + 1, N), vec_example)
            vecresq = np.dot(matrix_powmod(U_example, (N ** 2 + 1)//r4, N), vec_example)
            suits_a = (vecres[0] % N == 0 and vecresq[0] % N != 0)
            exists_a = exists_a or suits_a
        iters += 1
    #сноска 4
    prime = prime and exists_a
    if not prime:
        return prime
    return prime




k = [-1, # false
     11099088999991232491,  # true
     1109908899999121133, # true
     1109908899999121111, #false
     0, # false
     1, # false
     2, # true
     3, # true
     10, # false
     17, # true
     15, # false
     35980201101257391549860360923563262525974949247991832187257385201689, # true
     1889, # true
     1900, # false
     1901, # true
     1109908899999123251, # true
     359802011012573915498603609235632625259749492479918321872573852016899, # true
     35980201101257391549860360923563262525974949247991832187257385201688, #false
     35980201101257391549860360923563262525974949247991832187257385202287, #true
     3598020110125739154986036092356326252597494924799183218725738520227, # false
     3598020110125739154986036092356326252597494924799183218725738520331, #true
     12345678987654322345678765432345676573, #true
     345123123267672138238128321372718382137127321839123, #false
     3451231232676721382381283213727183821371273218391234, #false
     34512312326767213823812832137271838213712732183912341232421412421412412421412, #false
     34512312326767213823812832137271838213712732183912341232421412421412412421327, #true
     35980201101257391549860360923563262525974949247991832187257385201689000992721312321323213231231312313, #false
     35980201101257391549860360923563262525974949247991832187257385201689000992721312321323213231231311921, #true
     35980201101257391549860360923563262525974949247991832187257385202287] # true
for i in k:
    n = i #35980201101257391549860360923563262525974949247991832187257385201689 - from example
    print('---------------------------')
    print(i)
    print(is_prime(n))
print('---------------------------')