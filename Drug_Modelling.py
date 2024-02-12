import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def SIR(p1,p2,p3,S0,I0,R0,i):

    vect_initial = np.array([S0,I0,R0]).T
    
    output = pd.DataFrame(columns = ['day', 'S', 'I', 'R'])
    j = 0
    vect = vect_initial.copy()
    
    while j < i:
        vectplus1 = np.zeros([3])
        S = vect[0]
        I = vect[1]
        R = vect[2]
        vectplus1[0] = S - p1 *S*I
        vectplus1[1] = I +p1 *S*I + p2*R*I - p3*I
        vectplus1[2] = R - p2*R*I + p3*I
        output.loc[j] = np.concatenate((np.array([j]),vectplus1.T))
        vect = vectplus1.copy()
        j+=1
    					
    plt.plot(output['day'],output['S'])
    plt.plot(output['day'],output['I'])
    plt.plot(output['day'],output['R'])
    plt.title('p1=' + str(p1) +
              ', p2=' + str(p2) +
              ', p3=' +str(p3))
    plt.legend(['S','I','R'])
    plt.show()

    
def SUAR(p1,p2,p3, p4, p5 ,S0,U0,A0,R0,i):

    vect_initial = np.array([S0,U0,A0,R0]).T
    
    output = pd.DataFrame(columns = ['day', 'S', 'U','A', 'R'])
    j = 0
    vect = vect_initial.copy()
    
    while j < i:
        vectplus1 = np.zeros([4])
        S = vect[0]
        U = vect[1]
        A = vect[2]
        R = vect[3]
        vectplus1[0] = S -p1*S*(A+U) + p2*U
        vectplus1[1] = U +p1*S*(A+U) - p2*U - p3*U
        vectplus1[2] = A +p3*U + p5*(A+U)*R - p4*A
        vectplus1[3] = R + p4*A -p5*(A+U)*R
        output.loc[j] = np.concatenate((np.array([j]),vectplus1.T))
        vect = vectplus1.copy()
        j+=1
    					
    plt.plot(output['day'],output['S'])
    plt.plot(output['day'],output['U'])
    plt.plot(output['day'],output['A'])
    plt.plot(output['day'],output['R'])
    plt.title('p1=' + str(p1) +
              ', p2=' + str(p2) +
              ', p3=' +str(p3) +
              ', p4=' +str(p4) +
              ', p5=' +str(p5))
    plt.legend(['S','U','A','R'])
    plt.show()


def SUAR_20(aging_rate,birth_rate,death_rate, I, P, X0, steps):

    C = pd.DataFrame(columns=['s', 't', 'a', 'r'], dtype='float128')
    T = pd.DataFrame(columns=['s', 't', 'a', 'r'], dtype='float128')
    Y = pd.DataFrame(columns=['s', 't', 'a', 'r'], dtype='float128')
    A = pd.DataFrame(columns=['s', 't', 'a', 'r'], dtype='float128')
    S = pd.DataFrame(columns=['s', 't', 'a', 'r'], dtype='float128')
    
    (n_age_groups, n_status) = np.shape(X0)
    itterations  = 0
    while itterations < steps:
        C.loc[itterations] = X0[0]
        T.loc[itterations] = X0[1]
        Y.loc[itterations] = X0[2]
        A.loc[itterations] = X0[3]
        S.loc[itterations] = X0[4]
    
        i = 0
    
        IE = np.matmul(I, X0)
    
        Xn = np.zeros(np.shape(X0))
    
        while i < n_age_groups:
                group_vect = X0[i].copy()
                s = group_vect[0]
                t = group_vect[1]
                a = group_vect[2]
                r = group_vect[3]
                interaction_neg = np.sum(IE[i][[1, 2]])
                interaction_pos = np.sum(IE[i][[0, 3]])
    
                Xn[i][0] = s - (s * P[i][0] * interaction_neg) + t * P[i][1] * interaction_pos
                Xn[i][1] = t + (s * P[i][0] * interaction_neg) - t * P[i][1] * interaction_pos - t * P[i][2]
                Xn[i][2] = a + t * P[i][2] + r * P[i][4] * interaction_neg - a * P[i][3] * interaction_pos
                Xn[i][3] = r - r * P[i][4] * interaction_neg + a * P[i][3] * interaction_pos
    
                i += 1
    
        Xn_aged = np.zeros(np.shape(X0))
    
        i = 0
        j = 0
    
        while i < n_age_groups:
                j = 0
                while j < n_status:
                     if i == 0:
                            if j == 0:
                                Xn_aged[i][j] = Xn[i][j] * (1-aging_rate-death_rate[i][j])+birth_rate*sum(Xn[2]+ Xn[3])
                     j += 1           
                i+=1

        itterations += 1

def SUAR_5(p1,p2,p3, p4, p5,b,d,S0,U0,A0,R0,i):

    vect_initial = np.array([S0,U0,A0,R0]).T
    
    output = pd.DataFrame(columns = ['day', 'S', 'U','A', 'R'])
    j = 0
    vect = vect_initial.copy()
    
    while j < i:
        vectplus1 = np.zeros([4])
        S = vect[0]
        U = vect[1]
        A = vect[2]
        R = vect[3]
        vectplus1[0] = S - p1*S*(A+U) + p2*U*(S+R)  - d[0]*S + b*(S+U+A+R) 
        vectplus1[1] = U + p1*S*(A+U) - p2*U*(S+R) - p3*U - d[1]*U
        vectplus1[2] = A + p3*U + p5*(A+U)*R - p4*A*(R+S) - d[2]*A
        vectplus1[3] = R + p4*A*(S+R) - p5*(A+U)*R - d[3]*R

        output.loc[j] = np.concatenate((np.array([j]),vectplus1.T))
        vect = vectplus1.copy()
        j+=1
        
    output['pop'] = output['S'] + output['U'] + output['A'] + output['R']
    plt.plot(output['day'],output['S']/output['pop'])
    plt.plot(output['day'],output['U']/output['pop'])
    plt.plot(output['day'],output['A']/output['pop'])
    plt.plot(output['day'],output['R']/output['pop'])
    plt.title('Portion: p1=' + str(p1) +
              ', p2=' + str(p2) +
              ', p3=' +str(p3) +
              ', p4=' +str(p4) +
              ', p5=' +str(p5))
    plt.legend(['S','U','A','R'])
    plt.show()

def SUAR_norm(p1,p2,p3, p4, p5,d,S0,U0,A0,R0,i,h, plot = True):

    vect_initial = np.array([S0,U0,A0,R0]).T

    output = pd.DataFrame(columns = ['day', 's', 'u','a', 'r'])
    j = 0
    vect = vect_initial.copy()
    
    ds = d[0]
    du = d[1]
    da = d[2]
    dr = d[3]
    
    while j < i:
        vectplus1 = np.zeros([4])
        s = vect[0]
        u = vect[1]
        a = vect[2]
        r = vect[3]
        vectplus1[0] = s+h*((-s*p1*(u+a)+u*p2*(s+r) -ds*s + 1)-s *(1 - ds*s - du*u - da*a - dr*r))
        vectplus1[1] = u+h*((s*p1*(u+a)-u*p2*(s+r) - u*p3 -du*u )-u*(1 - ds*s - du*u - da*a - dr*r))
        vectplus1[2] = a+h*((u*p3 -a*p4*(s+r) + r*p5*(u+a) -da*a )-a*(1 - ds*s - du*u - da*a - dr*r))
        vectplus1[3] = r+h*((a*p4*(s+r) - r*p5*(u+a) -dr*r )-r*(1 - ds*s -du*u - da*a - dr*r))
        output.loc[j] = np.concatenate((np.array([j*h]),vectplus1.T))
        vect = vectplus1.copy()
        j+=1
    if plot:    
        plt.plot(output['day'],output['s'])
        plt.plot(output['day'],output['u'])
        plt.plot(output['day'],output['a'])
        plt.plot(output['day'],output['r'])
        plt.title('Portion: p1=' + str(p1) +
                  ', p2=' + str(p2) +
                  ', p3=' +str(p3) +
                  ', p4=' +str(p4) +
                  ', p5=' +str(p5))
        plt.legend(['S','U','A','R'])
        plt.show()
        return
    else:
        return [s,u,a,r]
