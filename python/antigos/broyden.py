from numpy import array, transpose, matmul
from numpy.linalg import solve, inv

class Broyden:

    def __init__ (self, h=10e-5, metodo_sistema_linear=solve):
        self.h = h
        self.metodo_sistema_linear = metodo_sistema_linear

    def derivada_parcial (self, f, ph, p):
        return (f(*ph)-f(*p))/self.h

    def gradiente (self, f, p):
        grad = []

        for i in range(len(p)):
            p0 = [*p]
            p0[i] += self.h
            der_par = self.derivada_parcial(f, p0, p)
            grad.append(der_par)

        return grad

    def jacobiana (self, F, p):
        Jac = []

        for f in F:
            grad = self.gradiente(f, p)
            Jac.append(grad)
        
        return Jac

    def resolver_sistema (self, A, b):
        return self.metodo_sistema_linear(A,b)

    def norma_2 (self, x):
        return sum(xi**2 for xi in x)**.5

    def norma_infinito (self, y):
        return max(abs(x) for x in y)

    def norma_infinito_matriz (self, M):
        return max(sum(abs(x) for x in Mi) for Mi in M)

    def numero_condicao (self, M, norma):
        return norma(M)

    def funcao (self, F, x):
        return array([f(*x) for f in F])

    def passo_broyden (self, F, p0, Jac=[], erro=10e-5):
        # matriz jacobiana
        if len(Jac) == 0: J = array(self.jacobiana(F, p0))
        else: J = array([[jij(*p0) for jij in ji] for ji in Jac])

        invA = inv(J)       
        print(invA) 
        x0 = p0
        F0 = self.funcao(F, x0)
        s = -matmul(invA,F0)
        x1 = x0 + s

        qntdPassos = 0

        while True:
            qntdPassos += 1

            F0 = self.funcao(F, x0)
            F1 = self.funcao(F, x1)
            
            y = F1 - F0
            
            st = transpose(s)
            
            q = matmul(st, invA)
            k = matmul(q,y)
            
            interno = (1/k)*(s - matmul(invA,y))
            invA = invA + matmul(interno,st)*invA
            
            s = -matmul(invA, F1)
            x0 = x1
            x1 = x1 + s
            
            if self.norma_infinito(s) < erro:
                return x1, self.norma_infinito(s)
        
# from math import sin, cos, exp, pi

# f1 = lambda x1, x2, x3: 3*x1-cos(x2*x3)-0.5
# f2 = lambda x1, x2, x3: x1**2 - 81*(x2+.1)**2+sin(x3)+1.06
# f3 = lambda x1, x2, x3: exp(-x1*x2)+20*x3+(10*pi-3)/3

# F = [f1, f2, f3]

# x0 = [0.1, 0.1, -0.1]

# broyden = Broyden()

# x2 = broyden.passo_broyden(F, x0, erro=10e-6)
# print(x2)