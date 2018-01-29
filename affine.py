#!/usr/bin/python
# coding: UTF-8
import numpy as np
import matplotlib.pyplot as plt
import math

def init():
    f = open('kondate.txt')
    temp = f.read()
    f.close()
    temp = temp.split()
    temp.insert(284, 1)
    temp.insert(285, 1)
    temp.append(1)
    temp.append(1)
    data = np.array(temp)
    data = np.reshape(data, (23, 13))
    data = data.astype(np.float64)

    Aori = data[0:21, 0:12] / 100.0
    bori = data[0:21, 12] / 100.0
    bori = bori.reshape((21, 1))
    cori = data[21, 0:11]
    cori = cori.reshape((11, 1))

    return [Aori, bori, cori]

def make_graph(obj_list, gamma):
    X = np.arange(len(obj_list) - 1)
    Y = np.array(obj_list[:-1])
    Y = np.log10(np.absolute(Y - obj_list[len(obj_list) - 1]))
    Y = Y.reshape(Y.size)
    ax = plt.subplot()
    ax.plot(X, Y)
    ax.set_xlabel("updata_num")
    ax.set_ylabel("log10(error)")
    ax.set_title("gamma:{}".format(gamma))
    plt.show()


def affine(A, b, c, y0, gamma, epsilon):
    y_list = []
    x_list = []
    obj_list = []
    m = len(A); n = len(A[0])
    y = y0
    while (True):
        y_list.append(y)
        obj_list.append(np.dot(b.T, y)[0][0])
        S2 = np.zeros((n, n))
        for i in range(n):
            a = A[0:m, i]
            a = a.reshape((m, 1))
            S2[i][i] = c[i] - np.dot(a.T, y)
            S2[i][i] = 1 / S2[i][i]**2
        AS2A = np.dot(np.dot(A, S2), A.T)
        x = np.dot(S2, np.dot(A.T, np.dot(np.linalg.inv(AS2A), b)))
        for i in range(len(x)):
            if (abs(x[i]) < 0.000000001):
                x[i] = 0
        x_list.append(x)
        L = np.linalg.cholesky(AS2A)
        u = np.linalg.solve(L, b)
        dy = np.linalg.solve(L.T, u)

        tmin = 0; tmax = 1000000000.0
        while (tmax - tmin > 0.0000001):
            t = (tmin + tmax) / 2.0
            if (np.min(c - np.dot(A.T, y + t*dy)) >= 0):
                tmin = t
            else:
                tmax = t

        t = tmin
        nexty = y + gamma * t * dy
        delta = np.dot(b.T, nexty) - np.dot(b.T, y)
        if (np.max(delta) < epsilon):
            y_list.append(nexty)
            obj_list.append(np.dot(b.T, nexty))
            return [y_list, obj_list, x_list]
        else:
            y = nexty
            
def problem3():
    res = init()

    Aori = res[0]
    bori = res[1]
    cori = res[2]
    
    A = np.hstack((-Aori[::,:-1], -np.eye(len(Aori))))
    b = -bori
    c = np.vstack((-cori, np.zeros((len(Aori), 1))))
  
    gamma = 0.95
    y0 = 1000*np.ones((len(A), 1))
    
    [y_list, obj_list, x_list] = affine(A, b, c, y0, gamma, 0.00001)

    print "双対推定"
    for i in range(len(x_list)):
        if (i % 5 == 0):
            print x_list[i]
            print

    print x_list[-1]
        
    make_graph(obj_list, gamma)
    print "グラム数"
    print y_list[-1].astype(np.int32)
    print "合計金額"
    print -obj_list[-1]

def problem4():
    res = init()
    Aori = res[0]
    bori = res[1]
    cori = res[2]
    
    m = len(Aori)
    n = len(Aori[0])

    Aori[:, -1] *= -1
    A = np.hstack((-Aori, -np.eye(m)))
    b = -bori
    cori = np.vstack((cori, -200))
    c = np.vstack((-cori, np.zeros((m, 1))))

    y0 = 0.1 * np.ones(b.shape)
    y0[0] = y0[1] = y0[2] = y0[3] = y0[10] = y0[14] = y0[15] = y0[17] = y0[18] = y0[20] = 1000
    y0[13] = 20

    gamma = 0.95
    
    [y_list, obj_list, x_list] = affine(A, b, c, y0, gamma, 0.00001)


    make_graph(obj_list, gamma)
        
    print "グラム数"
    print y_list[-1].astype(np.int32)
    print "合計金額"
    print -obj_list[-1]
        
def problem5():
    res = init()

    Aori = res[0]
    bori = res[1]
    cori = res[2]
    
    Aori[4] = Aori[5] = Aori[6] = Aori[7] = Aori[11] = Aori[12] = Aori[13] =  np.zeros(Aori[0].shape)
    
    A = np.hstack((-Aori[::,:-1], -np.eye(len(Aori))))
    b = -bori
    c = np.vstack((-cori, np.zeros((len(Aori), 1))))
  
    gamma = 0.66
    y0 = 1000*np.ones((len(A), 1))
    
    [y_list, obj_list, x_list] = affine(A, b, c, y0, gamma, 0.00001)
    
    make_graph(obj_list, gamma)
    print "グラム数"
    print y_list[-1].astype(np.int32)
    print "合計金額"
    print -obj_list[-1]
