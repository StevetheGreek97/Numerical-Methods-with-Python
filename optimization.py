
from math import sqrt

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
            

    def golden_search(f, a, b, tol=1e-8):
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
    
    def maximumsearch(f, a, b, steps):
        """lo, hi = maximumsearch(f, a, b, steps).
        Searches the interval (a,b) in a number of steps for
        the bounds (lo,hi) of the maxima of f(x).
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
            
            if f_md > f_lo and f_md >= f_hi:  
                yield lo, hi
                
            lo, f_lo = md, f_md    
            md, f_md= hi, f_hi

    def golden_search_max(f, a, b, tol=1e-8):
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
    #     print(count) 
        return (lo + hi) / 2.0