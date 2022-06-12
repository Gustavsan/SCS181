#!/usr/bin/env sage

from sage.all import *
import random
import time

def simple_case_alg(a, y, p, factor_list):
    if len(factor_list) != 1 and factor_list[0][0] !=2:
        raise Exception('Wrong value of argument factor_list')
    log_p = factor_list[0][1]
    z_i = y
    B_i = Mod(a ^ (-1), p) 
    m_i = (p-1)/2
    base = 1
    x = 0
    for i in range(0, log_p):

        q = Mod(z_i ^ m_i, p)
        if q == -1: 
            x += base
            z_i = Mod(z_i * B_i, p)
        elif q != 1:
            raise Exception('Wrong value of q (not 1 and not -1). Found: %s', str(q))
        base *= 2
        m_i /= 2
        B_i = Mod(B_i ^ 2, p)
    return x


def get_r_ij(a, j, q, p):
    return Mod(a ^ (((p-1)/q)*j), p)

def get_r_ij_dict(a, p, factor_list):
    r = dict()
    for i in range(len(factor_list)):
        factor = factor_list[i]
        q = factor[0]
        r[q] = []
        for j in range(0, q):
            r_ij = get_r_ij(a, j, q, p)
            r[factor[0]].append(r_ij)
    return r

def common_case_alg(a, b, p, factor_list):
    r = get_r_ij_dict(a, p, factor_list)

    x = []
    q_x = []
    for i in range(len(factor_list)):
        factor = factor_list[i]
        q = factor[0]
        q_a = q
        base_additional = 1
        x_i = 0
        for j in range(0, factor[1]):
            result = Mod(Mod(b * base_additional, p) ^ ((p-1)/q_a), p)
            for k in range(0, len(r[q])):
                if r[q][k] == result:
                    pre_calc = Mod(a^((-1) * k * q^j), p)
                    base_additional = Mod(base_additional * pre_calc, p)
                    x_i = Mod(x_i + k * q^j, p)
                    break
            q_a *= q
        x.append(sage.rings.integer.Integer(x_i))
        q_x.append(q^factor[1])
    return CRT_list(x, q_x)

print("Example for a=2, y=3, p=5")
p = 5
print(simple_case_alg(2, 3, p, list(factor(p-1))))


print("Example for a=3, y=4, p=7")
p = 7
print(common_case_alg(3, 5, p, list(factor(p-1))))

def check_alg(alg, sample, factorization):
    return alg(sample[0], sample[2], sample[3], factorization)

def check_answer(sample, result):
    #      result == x         or     a^result mod p == b
    return int(result) == int(sample[1]) or int(Mod(sample[0]^result, sample[3])) == int(sample[2])


def generate_samples(list_p, count):
    samples = []
    for p in list_p:
        a = GF(p, modulus="primitive").gen()
        for i in range(count):
            x = randint(2, p-2)
            b = int(Mod(a^x, p))
            samples.append((a, x, b, p))
    return samples

def run_timetest(samples, first_alg, second_alg):
    f_count = 0
    f_sum = 0
    s_count = 0
    s_sum = 0
    for sample in samples:
        factorization=list(factor(sample[3]-1))
        # first algorithm
        start_time_first = time.time()
        result = check_alg(first_alg, sample, factorization)
        end_time_first = time.time()
        # check answer
        if not check_answer(sample, result):
            print("first have ERROR", sample, result)
        else:
            f_count +=1 
            f_sum += (end_time_first - start_time_first)
            
        # second algorithm
        start_time_second = time.time()
        result = check_alg(second_alg, sample, factorization)
        end_time_second = time.time()
        # check answer
        if not check_answer(sample, result):
            print("second have ERROR", sample, result)
        else:
            s_count +=1 
            s_sum += (end_time_second - start_time_second)
    # report
    report_template = "Algorithm: %s.\n Count: %d. All time: %f. Avg time: %f"
    print(report_template % (first_alg.__name__ , f_count, f_sum, (f_sum/f_count)))
    print(report_template % (second_alg.__name__ , s_count, s_sum, (s_sum/s_count)))


print("Run timetest with simple and common case on p=[65537] and counts=10000")
samples = generate_samples([65537],10000)
run_timetest(samples, simple_case_alg, common_case_alg)


print("Run timetest with simple and common case on p=[65537] and counts=100000")
samples = generate_samples([65537],100000)
run_timetest(samples, simple_case_alg, common_case_alg)


def dummy_alg(a, b, p, factor_list):
    result = 1
    degree = 0
    while degree < p:
        if result == b:
            return degree
        result = int(Mod(result*a, p))
        degree += 1
    raise Exception("Not found result for %d^x = %d (mod %d)" % (a, b, p))


def generate_list_of_p(p_min, p_max, count):
    p = []
    i = 0
    fail = 0
    while i < count:
        random_p = randint(p_min, p_max)
        next_p = next_prime(random_p)
        if next_p <= p_max:
            p.append(next_p)
            i += 1
        else:
            fail += 1
        if fail > count and len(p) == 0:
            raise Exception("Bad attempt to generate p")
    return p

print("Run timetest with common and dummy alg on p in range(10^2, 10^4) and counts=1000")
list_p = generate_list_of_p(10^2, 10^4, 10)
samples = generate_samples(list_p, 100)
run_timetest(samples, common_case_alg, dummy_alg)

#print("Run timetest with common and dummy alg on p in range(10^4, 10^6) and counts=1000")
# too slow... around 1100 sec
#list_p = generate_list_of_p(10^4, 10^6, 10)
#samples = generate_samples(list_p,100)
#run_timetest(samples, common_case_alg, dummy_alg)





