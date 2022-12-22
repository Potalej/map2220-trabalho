from numpy import array
from numpy.linalg import solve

class Newton:

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

    def norma_infinito (self, y):
        return max(abs(x) for x in y)

    def norma_infinito_matriz (self, M):
        return max(sum(abs(x) for x in Mi) for Mi in M)

    def numero_condicao (self, M, norma):
        return norma(M)

    def passo (self, F, p0, Jac=[]):
        # matriz jacobiana
        if len(Jac) == 0: J = array(self.jacobiana(F, p0))
        else: J = array([[jij(*p0) for jij in ji] for ji in Jac])
        
        # - F(x0)
        Fs = - array([f(*p0) for f in F])

        # J(x0)*y = -F(x0)
        y = self.resolver_sistema(J, Fs)

        # aplica o passo
        p = [p0[i]+y[i] for i in range(len(p0))]

        return p, y, J