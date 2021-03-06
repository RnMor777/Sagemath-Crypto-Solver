import time
import timeout_decorator
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import os
from rich.console import Console
from reportlab.pdfgen import canvas
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase import pdfmetrics
from reportlab.lib import colors

fileName = 'report.pdf'
documentTitle = 'RSA Decoding Report'
title = 'RSA Decoding Report'
subTitle = 'Program by Ryan'
image = 'img.png'
pdf = canvas.Canvas(fileName)
pdf.setTitle(documentTitle)
pdfmetrics.registerFont(TTFont('abc', 'UbuntuMono-R.ttf'))
console = Console()
TIME_OUT = 3

# Determines the divisors of p-1
def order(a,p):
    if is_prime(p) == False or a%p == 0:
        return False
    for k in divisors(p-1):
        if power_mod(a,k,p) == 1:
            return k

# Determines if the order of a%p = p-1
def is_gernerator(a,p):
    if is_prime(p) == False:
        return False
    if order(a,p) == p-1:
        return True

#Chinese Remainder Theorem Solver
def chinese_remainder_theorem(a_values, m_values):
    for i in range(len(a_values)):
        for j in range(i+1,len(a_values)):
            if gcd(m_values[i], m_values[j]) != 1:
                return "fail"
    prod_M = prod(m_values)
    x = sum([a_values[i]*(prod_M/m_values[i])*inverse_mod(prod_M//m_values[i],m_values[i]) for i in range(len(a_values))])
    return x%prod_M

# DLP Brute Force
def discrete_log_brute_force(g,h,p):
    h = h%p
    for k in [0..p-2]:
        if power_mod(g,k,p) == h:
            return k
    return "No Solution"

# Baby Step - Giant Step DLP Algorithm
def babystep_giantstep (g,h,p):
    N = order(g,p)
    n = 1 + floor(sqrt(N))
    list1 = [1]
    list2 = [h]
    u = power_mod(inverse_mod(g,p),n,p)
    a = 1
    b = 1
    for i in [1..n]:
        a = (a*g)%p
        list1.append(a)
        b = (b*u)%p
        list2.append((h*b)%p)
    shared = first_common_element(list1, list2)
    if shared == None:
        return "No Solution"
    i = list1.index(shared)
    j = list2.index(shared)
    return (i+j*n)%N

# Pohlig Hellman DLP Algorithm
def pholig_hellman (g,h,p):
    N = order(g,p)
    factored_N = factor(N)
    q_list = [(factored_N[i][0])^(factored_N[i][1]) for i in range(len(factored_N))]
    g_list = [power_mod(g,N//q_list[i],p) for i in range(len(factored_N))]
    h_list = [power_mod(h,N//q_list[i],p) for i in range(len(factored_N))]   
    a_list = [babystep_giantstep(ghpair[0], ghpair[1],p) for ghpair in zip(g_list,h_list)]
    return chinese_remainder_theorem(a_list, q_list)

# Set Collison
def first_common_element (X,Y):
    for i in X:
        if i in Y:
            return i

# Miller Rabin Primality Test
def miller_rabin_witness (d, n):
    a = 2 + randint(2, n-1)
    x = power_mod (a, d, n)

    if (x==1 or x == n-1):
        return True

    while (d != n-1):
        x = power_mod(x,2,n)
        d *= 2

        if (x == 1):
            return False
        if (x == n-1):
            return True

    return False

# Another Implementation of prime testing
def is_prime2(n, k):
    d = n -1
    while (d&1 == 0):
            d = d >> 1

    for i in range(k):
        if (miller_rabin_witness(d,n) == False):
            return False
    return True

# Pollard P-1 Factorization Method
def pollard_pminus1(n,a,B):
    if (is_prime(n) == True):
        return n
    if (mod(n,2) == 0):
        return 2
    j = 2
    t = time.time()
    t2 = time.time() - t
    while (t2 < 30):
        a = power_mod(a,j,n)
        d = int(gcd(a-1,n))
        if (d > 1 and d < n):
            return Integer(d)
        j += 1
        t2 = time.time() - t
    return 'Fail'

# Fermat Factorization V1
def fermat_factor_v1 (N):
    k = 0
    while True:
        if (is_square(N+k^2) == True):
            return sqrt(N+k^2) + k
        k += 1

# Fermat Factorization V2
def fermat_factor_v2 (N):
    t0 = ceil(sqrt(N))
    k = 0
    while True:
        if (is_square((t0+k)^2 - N) == True):
            return t0+k+sqrt((t0+k)^2 - N)
        k += 1

#Euler Phi
def phi(n):
    result = n
    N = n
    p = 2
    while(p * p <= n):
        if (n % p == 0):
            while (n % p == 0):
                n = int(n / p)
            result -= int(result / p)
        p += 1
    if (n > 1):
        result -= int(result / n)
    return result

def quad_res(N):
    L = ceil(e^(sqrt(ln(N)*ln(ln(N)))))
    B = ceil(L^(1/sqrt(2)))
    T = ceil(L^sqrt(2))

    t0 = ceil(sqrt(N))
    list_of_as_and_cs = []
    for k in range(T):
        a = t0 + k
        c = power_mod(t0+k,2,N)
        list_of_as_and_cs.append([a,c,c,[]])

    prime_base = prime_range(B+1)

    for i in range(len(prime_base)):
        p = prime_base[i]
        for j in range(len(list_of_as_and_cs)):
            list_of_as_and_cs[j][3].append(0)
            while (mod(list_of_as_and_cs[j][2],p) == 0):
                list_of_as_and_cs[j][2] = list_of_as_and_cs[j][2]/p
                list_of_as_and_cs[j][3][i] = list_of_as_and_cs[j][3][i] + 1

    list_of_as = []
    list_of_cs = []
    list_of_cs_exponents = []
    for j in range(len(list_of_as_and_cs)):
        if list_of_as_and_cs[j][2] == 1:
            list_of_as.append(list_of_as_and_cs[j][0])
            list_of_cs.append(list_of_as_and_cs[j][1])
            list_of_cs_exponents.append(list_of_as_and_cs[j][3])

    M = matrix(GF(2),list_of_cs_exponents)
    list_of_ws = basis(kernel(M))

    for i in range(len(list_of_ws)):
        a_prod = product([a for a,k in zip(list_of_as,list_of_ws[i]) if k==1])
        c_prod = product([c for c,k in zip(list_of_cs,list_of_ws[i]) if k==1])
        if gcd(a_prod - sqrt(c_prod),N) > 1 and gcd(a_prod - sqrt(c_prod),N) < N:
            return gcd(a_prod - sqrt(c_prod),N)

def driver_rsa ():
    pfin = 0
    tot = 0
    # e = int(input("Enter the e value: "))
    # c = int(input("Enter the c value: "))
    N = int(input("Enter the N value: "))
    performance = []
    
    print ("\nRSA Analysis\n")

    with console.status("[bold green] Pollard P-1") as status:
        @timeout_decorator.timeout(TIME_OUT)
        def pm1():
            t = time.time()
            p = pollard_pminus1(N,2,30)
            t = time.time() - t
            performance.append(t)
            return p
        try:
            pfin = pm1()
        except Exception:
            performance.append(TIME_OUT)
        console.log("Pollard P-1 Complete")

    with console.status("[bold green] Fermat V1") as status:
        @timeout_decorator.timeout(TIME_OUT)
        def fem1():
            t = time.time()
            p = fermat_factor_v1(N)
            t = time.time() - t
            performance.append(t)
            return p
        try:
            pfin = fem1()
        except Exception:
            performance.append(TIME_OUT)
        console.log("Fermat V1 Complete")

    with console.status("[bold green] Fermat V2") as status:
        @timeout_decorator.timeout(TIME_OUT)
        def fem2():
            t = time.time()
            p = fermat_factor_v2(N)
            t = time.time() - t
            performance.append(t)
            return p
        try:
            pfin = fem2()
        except Exception:
            performance.append(TIME_OUT)
        console.log("Fermat V2 Completed")

    with console.status("[bold green] Euler's Totient") as status:
        @timeout_decorator.timeout(TIME_OUT)
        def eultot():
            t = time.time()
            eulphi = euler_phi(N)
            t = time.time() - t
            b = (1+N-eulphi)
            p = (b + sqrt(b^2-4*1*N))/2
            performance.append(t)
            return p, eulphi
        try:
            pfin, tot = eultot()
        except Exception:
            performance.append(TIME_OUT)
        console.log("Euler's Totient Complete")

    with console.status("[bold green] Elliptic Curve") as status:
        @timeout_decorator.timeout(TIME_OUT)
        def qr():
            t = time.time()
            x = ecm.factor(N)
            t = time.time() - t
            performance.append(t)
            x.sort()
            return int(x[0])
        try:
            pfin = qr()
        except Exception:
            performance.append(TIME_OUT)
        console.log("Elliptic Curve Complete")

    with console.status("[bold green] Built-in") as status:
        @timeout_decorator.timeout(TIME_OUT)
        def built_solve():
            t = time.time()
            fact = factor(N)
            t = time.time() - t
            performance.append(t)
            return int((str(fact).split("*"))[0])
        try:
            pfin = built_solve()
        except Exception:
            performance.append(TIME_OUT)
        console.log("Built-in Complete")

    if (pfin == 0):
        print ("No Results - Try a higher time out value")
        exit(1)

    objects = ("Pollard", "Fermat V1", "Fermat V2", "Totient", "E-Curve", "Built-in")
    y_pos = np.arange(len(objects))
    plt.bar(y_pos, performance, align='center', alpha=0.5)
    for i in range(len(objects)):
        plt.text(i,performance[i],round(performance[i],4), ha="center")
    plt.xticks(y_pos, objects)
    plt.ylabel("Time Taken (seconds)")
    plt.xlabel("Algorithm")
    plt.title("Time of factorization algorithm")
    plt.savefig("img.png")

    pdf.setFont('abc', 40)
    pdf.drawCentredString(300, 770, title)
    pdf.setFont('abc', 20)
    pdf.drawCentredString(290, 720, subTitle)
    pdf.line(30, 700, 550, 700)
    pdf.drawInlineImage(image, 70, 350, width=440, height=330)
    pdf.line(30, 340, 550, 340)
    pdf.drawString (40, 310, "N = "+str(N))
    pdf.drawString (40, 290, "P = "+str(pfin))
    pdf.drawString (40, 270, "Q = "+str(N//pfin))
    pdf.drawString (40, 250, "?? = "+str(tot))
    pdf.drawString (40, 210, "P-1  = "+str(factor(pfin-1)))
    pdf.drawString (40, 190, "Q-1  = "+str(factor((N//pfin)-1)))

    pdf.save()

    os.remove("img.png")    
    os.system("okular report.pdf")

driver_rsa()
# vim:ft=python
