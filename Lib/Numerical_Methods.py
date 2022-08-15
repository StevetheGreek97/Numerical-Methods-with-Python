
from math import sqrt
import numpy as np

class Optimizers():

    def minimumsearch(f, a, b, steps):
        """lo, hi = minimumsearch(f, a, b, steps).
        Searches the interval (a,b) in a number of steps for
        the bounds (lo,hi) of the minima of f(x).
        """
        if steps < 2:
            steps = 2
        h = (b- a) / steps
        lo = a
        md = a + h 
        f_lo = f(lo)
        f_md = f(md)
        for i in range(2, steps +1):
            hi = a + h * i
            f_hi = f(hi)
            
            if f_md < f_lo and f_md <= f_hi:
                yield lo, hi
                
            lo, f_lo = md, f_md    
            md, f_md= hi, f_hi

    def maximumsearch(f, a, b, steps):
        """lo, hi = minimumsearch(f, a, b, steps).
        Searches the interval (a,b) in a number of steps for
        the bounds (lo,hi) of the minima of f(x).
        """
        if steps < 2:
            steps = 2
        h = (b- a) / steps
        lo = a
        md = a + h 
        f_lo = f(lo)
        f_md = f(md)
        for i in range(2, steps +1):
            hi = a + h * i
            f_hi = f(hi)
            
            if f_md > f_lo and f_md >= f_hi:  # if operant '>' algo finds maxima else finds minimum
                yield lo, hi
                
            lo, f_lo = md, f_md    
            md, f_md= hi, f_hi

    def min_golden_search(f, a, b, tol=1e-8):
        """minimum = golden_search(f, a, b, tol=...).
        Finds a minimum of the function f using golden
        section search, starting from bracketing values
        a and b and ending when |b-a| < tol.
        """
        factor = (3.0 - sqrt(5.0)) / 2.0

        x1 = a + factor * (b - a)
        x2 = b - factor * (b - a)
        f1 = f(x1)
        f2 = f(x2)
        while abs(b-a) > tol:
            if f1 < f2: # if operant '>' algo finds maxima else finds minimum
                b = x2
                x2 = x1
                f2 = f1
                x1 = a + factor * (b - a)
                f1 = f(x1)
            else:
                a = x1
                x1 = x2
                f1 = f2
                x2 =  b - factor * (b - a)
                f2 = f(x2)         
        return (a + b) / 2.0

    def max_golden_search(f, a, b, tol=1e-8):
        """maximum = golden_search(f, a, b, tol=...).
        Finds a minimum of the function f using golden
        section search, starting from bracketing values
        a and b and ending when |b-a| < tol.
        """
        factor = (3.0 - sqrt(5.0)) / 2.0

        x1 = a + factor * (b - a)
        x2 = b - factor * (b - a)
        f1 = f(x1)
        f2 = f(x2)
        while abs(b-a) > tol:
            if f1 > f2: # if operant '>' algo finds maxima else finds minimum
                b = x2
                x2 = x1
                f2 = f1
                x1 = a + factor * (b - a)
                f1 = f(x1)
            else:
                a = x1
                x1 = x2
                f1 = f2
                x2 =  b - factor * (b - a)
                f2 = f(x2)
        return (a + b) / 2.0

class Integrators():

    def trapezoid(f, a, b, n=1000):
        """df = trapezoid(f, a, b, n=...).
        Calculates the definite integral of the function f(x)
        from a to b using the composite trapezoidal rule with
        n subdivisions (with default n=...).
        """
        h = (b - a ) / n
        I  = f(a) + f(b) / 2
        for i in range(1, n-1):
            xi = a + i * h 
            I += f(xi)
        I *= h 
        return I 

    def rec_trapezoid(f, a, b, tol=1e-8):
        """df = recurcize trapezoid(f, a, b, tol=...).
        Calculates the definite integral of the function f(x)
        from a to b using the recursive trapezoidal rule with
        an absolute tolerance tol (with default 1e-8).
        """
        h = (b - a) # Interval size
        panels = 1 # No. of intervals
        I_old = (f(a)+f(b)) * h/2 
        
        while True:
            h /= 2
            panels *= 2

            I_new = 0.5 * I_old + sum(f(a +i*h) for i in range(1, panels, 2)) * h
            
            
            if abs(I_new - I_old) < tol:
                return I_new
            else:
                I_old = I_new

    def simpson(f, a, b, n=100):
        """df = simpson(f, a, b, n=...).
        Calculates the definite integral of the function f(x)
        from a to b using the composite Simpson's
        rule with n subdivisions (with default n=...).
        """
        
        n += n % 2 # force to be even
        
        h = (b -a) / n 
        
        I = f(a) + f(b)
        
        for i in range(1, n, 2):
            xi = a + i*h
            I += 4*f(xi)
            
        for i in range(2, n, 2):
            xi = a + i*h
            I += 2*f(xi)
        I *= h/3
        
        return I

    def romberg(f, a, b, tol = 1e-8):
        """df = simpson(f, a, b, tol=...).
        Calculates the definite integral of the function f(x)
        from a to b using Romberg integration based on the
        trapezoidal rule until a specified tolerance tol is
        reached (with default tol=...).
        """
        h = (b - a) # Interval size
        n = 1 # No. of intervals
        Rold = [ (f(a)+f(b)) * h/2 ]
        
        while True:
            h /= 2
            n *= 2 
            Rnew = [ 0.5 * Rold[0] + sum(f(a +o*h) for o in range(1, n, 2)) * h  ]
            factor = 1
            for R in Rold:
                factor *= 4
                Rnew.append( (factor*Rnew[-1] - R) / (factor-1) )
                
            if abs(Rnew[-1] - Rold[-1]) < tol:
                return Rnew[-1]
            Rold = Rnew

class Rootfinders():

    def rootsearch(f, a, b, steps):
        """lo, hi = rootsearch(f, a, b, steps).
        Searches the interval (a,b) in a number of steps for
        the bounds (lo,hi) of the roots of f(x).
        """
        h = (b - a) / steps
        f_lo = f(a)
        for step in range(steps):
            lo = a + step * h
            hi = lo + h
            f_hi = f(hi)
            if f_lo * f_hi <= 0.0:
                yield lo, hi
            f_lo = f_hi

    def regula_falsi(f, a, b, tol=1e-8, maxiter = 50):
        """
        """
        fa = f(a)
        if fa == 0:
            return a, fa
        
        fb = f(b)
        if fb == 0:
            return b, fb
        
        if fa * fb > 0.0:
            raise ValueError('Root is not bracketed')
        # as long as difference between my x values is bigger than the tolerance:
        
        while maxiter > 0:
            maxiter -= 1
            #calculate an new xn value
            xn = (a * fb - b * fa) / (fb - fa)
            fxn = f(xn)
            #print(fa, fxn, fb)

            if fxn == 0:
                return xn, fxn
            
            elif abs(fxn) < tol:
                return xn, fxn
            
            else:
                if (fxn * fa) > 0:
                    a, fa = xn, fxn

                elif (fxn * fb) > 0:
                    b, fb = xn, fxn

                else:
                    print('bad')
        print('no root found with given maximal number of iterations')

    def ridder(f,a,b,tol=1.0e-9):
        
        lo, f_lo = a, f(a)
        if f_lo == 0.0:
            return lo
        hi, f_hi = b, f(b)
        if f_hi == 0.0:
            return hi
        if f_lo * f_hi > 0.0:
            raise ValueError('Root is not bracketed')
            
        while abs(hi - lo) > tol:
            mid = (hi + lo) / 2.0
            f_mid = f(mid)
            
            s= np.sqrt(f_mid**2 - f_lo*f_hi)
            if s == 0.0:
                return None
            dx = (mid - lo)*f_mid/s
            if (f_lo - f_hi) < 0.0:
                dx = -dx
            x = mid + dx
            fx = f(x)
            if (f_mid * fx > 0):
                if (f_lo * fx < 0):
                    hi = x
                    f_hi = f(hi)
                else:
                    lo = x
                    f_lo = f(lo)
            else:
                lo, hi, f_lo, f_hi = mid, x, f_mid, fx
        return x
    def bisection(f, a, b, tol=1e-4):
        """root = bisection(f, a, b, tol=...).
        Finds a root of f(x) = 0 by bisection.
        The root must be bracketed in (a,b).
        """
        count = 0
        lo, f_lo = a, f(a)
        if f_lo == 0.0:
            return lo, f_lo
        hi, f_hi = b, f(b)
        if f_hi == 0.0:
            return hi, f_hi
        if f_lo * f_hi > 0.0:
            raise ValueError('Root is not bracketed')
        while abs(hi - lo) > tol:
    #         count += 1
            mid = (hi + lo) / 2.0
            f_mid = f(mid)
            if f_mid == 0.0:
                return mid, f_mid
            if (f_mid * f_hi > 0):
                hi = mid
                f_hi = f_mid
            else:
                lo = mid
                f_lo = f_mid

    def secant(f, a, b, tol= 1e-8):
        """root = secant(f, a, b, tol=...).
        Finds a root of f(x) = 0 by the secant method.
        """
        x1 = a
        f1 = f(x1)
            
        if f1 == 0:
            return x1
            
        x2 = b
        f2 = f(x2) 
        
        if f2 == 0:
            return x2
        while abs(x2 - x1) > tol:
            x3 = (f1 *x2 - f2 *x1) / (f1 - f2)
            f3= f(x3)
            if f3 == 0:
                return x3
            
            x1 = x2
            x2 = x3
            f1 = f2
            f2 = f3
        return x2

    def newton_raphson(f, df, a, b, tol=...):
        """root = newton_raphson(f, df, a, b, tol=....).
        Finds a root of f(x) = 0 by combining the Newton-Raphson
        method with bisection. The root must be bracketed in (a,b).
        Calls user-supplied functions f(x) and its derivative df(x).
        """

        x0 = (a + b) / 2
        f0 = f(x0)
        df0 = df(x0)
        while True:
            delta = f0/df0
            x1 = x0 - delta
            if abs(x1- x0) < tol:
                return x1
            x0 = x1
            f0 = f(x0)
            df0 = df(x0)


    def newton_raphson_der(f, a, b, tol=1e-4):
        """root = newton_raphson(f, df, a, b, tol=....).
        Finds a root of f(x) = 0 by combining the Newton-Raphson
        method with bisection. The root must be bracketed in (a,b).
        Calls user-supplied functions f(x) and its derivative df(x).
        """

        x0 = (a + b) / 2
        f0 = f(x0)
        df0 = central_derivative(f, x0, h=1e-6)
        while True:
            delta = f0/df0
            x1 = x0 - delta
            if abs(x1- x0) < tol:
                return x1
            x0 = x1
            f0 = f(x0)
            df0 = central_derivative(f, x0, h=1e-6)

class RungeKutta():
    def euler(f, y0, t, h):
        """xs, ys = euler(f, y0, x0, x1, steps).
        Euler's method for solving the
        initial value problem {y}' = {f(x,{y})},
        where {y} = {y[0],y[1],...,y[n-1]}.
        x0, y0 = initial conditions
        x1     = terminal value of x
        steps  = number of integration steps
        f      = user-supplied function that returns the
                array f(x,y) = {y’[0],y’[1],...,y’[n-1]}.
        """
        steps = int(t / h) + 1
        h = t / steps

        xs = np.linspace(0, t, steps + 1)
        y = y0
        for x in xs[:-1]:
            y = y + h * f(x, y)
  
        return y

    def heun(f, y0, x0, x1, steps):
        """xs, ys = heun(f, y0, x0, x1, steps).
        Heun's method for solving the
        initial value problem {y}' = {f(x,{y})},
        where {y} = {y[0],y[1],...,y[n-1]}.
        x0, y0 = initial conditions
        x1     = terminal value of x
        steps  = number of integration steps
        f      = user-supplied function that returns the
                array f(x,y) = {y’[0],y’[1],...,y’[n-1]}.
        """
        h = (x1 - x0) / steps
        xs = np.linspace(x0, x1, steps + 1)
        y = y0
        ys =[y]
        for x in xs[:-1]:
            k1 = h *  f(x, y)
            k2 = h * f(x + h, y + k1 )
            y = y +  0.5 * (k1 + k2)
            ys.append(y)
        return  xs, ys

    def Kutta(f, y0, x0, x1, steps):
        
        h = (x1 - x0) / steps
        xs = np.linspace(x0, x1, steps + 1)
        y = y0
        ys =[y]
        for x in xs[:-1]:
            k0 = f(x, y)
            k1 = f(x + (h/2) , y + (h/2)*k0)
            k2 = f(x + (h/2) , y + (-h)*k0 + 2*h*k1 )
            y = y + (h*(k0 + 4*k1 + k2))/ 6
            ys.append(y)
        return  xs, ys 

    def runge_kutta(f, y0, x0, x1, steps):
        """xs, ys = runge_kutta(f, y0, x0, x1, steps).
        4th-order Runge-Kutta method for solving the
        initial value problem {y}' = {f(x,{y})},
        where {y} = {y[0],y[1],...,y[n-1]}.
        x0, y0 = initial conditions
        x1     = terminal value of x
        steps  = number of integration steps
        f      = user-supplied function that returns the
                array f(x,y) = {y’[0],y’[1],...,y’[n-1]}.
        """
        h = (x1 - x0) / steps
        xs = np.linspace(x0, x1, steps + 1)
        y = y0
        ys =[y]
        for x in xs[:-1]:
            #Initial calculation
            k0 = h * f(x, y)
            # Middle calculations
            k1 = h * f(x + 0.5 * h, y + 0.5 *k0)
            k2 = h * f(x + 0.5 * h, y + 0.5 *k1)
            # End calculation
            k3 = h * f(x +h, y + k2)
            
            
            y = y + (k0 + 2.0 * k1 + 2.0 * k2 +k3) / 6
            ys.append(y)
        return  xs, ys

    def Butcher(f, y0, x0, x1, steps):
        h = (x1 - x0) / steps
        xs = np.linspace(x0, x1, steps + 1)
        y = y0
        ys =[y]
        for x in xs[:-1]:
            k1 = f(x, y)
            k2 = f(x + (h/4), y + (h/4)*k1        )
            k3 = f(x + (h/4), y + (h/8)*(k1 + k2) ) 
            k4 = f(x + (h/2), y + (-h/2)*k2 + k3*h  )
            k5 = f(x + ((12*h)/16), y + ((3*h)/16)*k1 + ((9*h)/16)*k4      )
            k6 = f(x + h, y +h*((-3/7)*k1  + (2/7)*k2 + (12/7)*k3 + (-12/7)*k4 + (8/7)*k5))
            
            y = y + (h*(7*k1 + 32*k3 + 12*k4 + 32*k5 + 7*k6))/90
            ys.append(y)
        return  xs, ys 

    def midpoint(f, y0, x0, x1, steps):
        h = (x1 - x0) / steps
        xs = np.linspace(x0, x1, steps + 1)
        y = y0
        ys =[y]
        for x in xs[:-1]:
            k1 = f(x, y)
            k2 = f(x + (h/2), y + (h/2)*k1)
            
            y = y + h*(k2)
            ys.append(y)
        return  xs, ys

    def SSPRK3(f, y0, x0, x1, steps):
        
        h = (x1 - x0) / steps
        xs = np.linspace(x0, x1, steps + 1)
        y = y0
        ys =[y]
        for x in xs[:-1]:
            k0 = f(x, y)
            k1 = f(x +h, y + k0*h)
            k2 = f(x +(h/2), y + (h/4)*h)
            y = y + (h*(k0 + k1 + 4*k2))/ 6
            ys.append(y)
        return  xs, ys 

    def SSPRK4(f, y0, t, steps):
        
        h = t / steps
        xs = np.linspace(0, t, steps + 1)
        y = y0
        ys =[y]
        for x in xs[:-1]:

            k1 = f(x, y)
            k2 = f(x + (h/2), y + (h/2)*k1             )
            k3 = f(x +  h,    y + (h/2)*(k1 + k2)      )
            k4 = f(x + (h/2), y + (h/6)*(k1 + k2 + k3) )
            
            
            y = y + (h*(k1 + k2 + k3 +3.0*k4))/ 6

        return  y

    def three_eight_rule(f, y0, x0, x1, steps):
        
        h = (x1 - x0) / steps
        xs = np.linspace(x0, x1, steps + 1)
        y = y0
        ys =[y]
        for x in xs[:-1]:
            #Initial calculation
            k1 = f(x, y)
            # Middle calculations
            k2 = f(x + (h/3)    , y + (h/3)*k1)
            k3 = f(x + ((2*h)/3), y -h*(k1/3 - k2))
            # End calculation
            k4 = f(x +h, y + h*(k1 - k2 + k3) )
            
            
            y = y + (h*(k1 + 3.0*(k2 + k3) +k4))/ 8
            ys.append(y)
        return  xs, ys

    def ralston(f, y0, t, steps):
        """y = ralston(f, y0, t, h)"""
        
        h = t / steps
        xs = np.linspace(0, t, steps + 1)
        y = y0
        ys =[y]
        for x in xs[:-1]:
            #Initial calculation
            k1 = f(x, y)
            # Middle calculations
            k2 = f(x + ((3*h)/4), y + ((3*h)/4)*k1)

            
            
            y = y + (k1*(1/3) + k2*(2/3))*h

        return  y