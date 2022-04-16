from multiprocessing.sharedctypes import Value
import sympy as sp


a = sp.symbols('alpha')


def _phase2golden(alpha_l, alpha_u, functionOfAlpha, error, k, debug): #implementação do SLIDE 13
    
    currentError = error + 1 # força o algorítimo a passar para o loop while
    k += 1
    
    while abs(currentError) > error: # Compara o erro atual com a tolerância desejada
        alpha = (alpha_u+alpha_l) / 2
        fk_current = functionOfAlpha.subs(a, alpha)
               
        tau = 0.618
        alpha_a = alpha_l + (1-tau) * (alpha_u-alpha_l)
        alpha_b = alpha_l + tau * (alpha_u-alpha_l)
        
        fka = functionOfAlpha.subs(a, alpha_a)
        fkb = functionOfAlpha.subs(a, alpha_b)
              
        if fka < fkb:
            alpha_u = alpha_b
        elif fka > fkb:
            alpha_l = alpha_a
        
        alpha = (alpha_u+alpha_l) / 2
        fk_new = functionOfAlpha.subs(a, alpha)
        
        currentError = abs(fk_current - fk_new)
        k += 1
        
        if k > 1000:
            return ValueError('Convergence not reached until 1000 iterations')

        if debug:
            print(f'Iteration {k}')
            print(f'Number of iterations {k}')
            print(f'optimum alpha = {alpha}')
            print(f'f(alpha*) = {fk_new}')

    #print(f'Number of iterations {k}')
    #print(f'optimum alpha = {alpha}')
    #print(f'f(alpha*) = {fk_new}')
    return alpha


def goldenSearch(functionOfAlpha, delta, debug, error):
    '''
    Parameters
    ----------
    functionOfAlpha : sympy.core.add.Add
        A sympy made equation.
        It should be a equation parametrized in function of alpha
    debug : bool
        To see printed partial results
    error : float, default 10e-4
        Error allowed in convergence
    delta : float , default 0.1
        Iteration step

    Returns
    -------
    alpha : float
        returns the optimum value of alpha
    '''


    currentError = 100 # Número grande qualquer
    k = 1 # número de steps
    q = 0 # usado no método golden search
    alpha = 0 # Delta acumulado
        
    while abs(currentError) > error: # Compara o erro atual com a tolerância desejada
        alpha0 = alpha + delta * 1.618 ** q
        alpha1 = alpha0 + delta * 1.618 ** (q+1)
        alpha2 = alpha1 + delta * 1.618 ** (q+2)
        
        fk0 = functionOfAlpha.subs(a, alpha0) #step k
        fk1 = functionOfAlpha.subs(a, alpha1) # step k+1
        fk2 = functionOfAlpha.subs(a, alpha2) # step k+1

        if fk1 < fk2 :
            alpha_l = alpha0
            alpha_u = alpha2
            alpha = _phase2golden(alpha_l, alpha_u, functionOfAlpha, error, k, debug) #final alpha
            break
        else:
            q += 1
            alpha = alpha0 # garante a continuação do algoritmo
       
        k += 1
        
        if debug:
            print(f'alpha0 = {alpha0}')
            print(f'fk0 = {fk0}')
            print(f'alpha1 = {alpha1}')
            print(f'fk1 = {fk1}')
            print(f'alpha2 = {alpha2}')
            print(f'fk2 = {fk2}')


        if k > 1000:
            return ValueError('Convergence not reached until 1000 iterations')
    return alpha