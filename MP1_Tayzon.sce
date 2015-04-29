//Keith Tayzon

a = -2              // Initial Guess
b = 0               // Initial Guess
tol = 1 * 10 ^ -10  // Maximum Error
i = 0               // Iteration

// The equation f(x) = x^3 + 2*tan(x) + 4
function[y] = f(x)
   y = x^3 + 2*tan(x) + 4
endfunction

// The derivative of f(x)
function[y] = f1(x)
   y = 3*x^2 + 2*(sec(x))^2
endfunction

// Implements the Interval Bisection method
function[a,b,i] = IBM(a,b,tol,i)
    mprintf('Interval Bisection method\n')
    while (b - a) > tol
        m = a + (b - a) / 2 // Midpoint between a and b
        mprintf('%i\ta = %.11f\tf(a) = %.11f\tb = %.11f\tf(b) = %.11f\tm = %.11f\tb-a = %.11f\n',i,a,f(a),b,f(b),m,b-a)
        if sign(f(a)) == sign(f(m)) then
            a = m
        else
            b = m
        end
        i = i + 1
    end
    mprintf('%i\ta = %.11f\tf(a) = %.11f\tb = %.11f\tf(b) = %.11f\tm = %.11f\tb-a = %.11f\n\n',i,a,f(a),b,f(b),m,b-a)
endfunction

// Implements Newton's method
function[xk,i] = NM(xk,tol,i)
    mprintf('Newton''s method\n')
    while i < 100
        h = - f(xk) / f1(xk)    
        xkp1 = xk + h           
        mprintf('%i\tx = %.11f\tf(x) = %.11f\tf''(x) = %.11f\th = %.11f\n',i,xk,f(xk),f1(xk),h)
        if abs(xkp1 - xk) < tol then
            break
        end
        xk = xkp1
        i = i + 1
    end
    mprintf('\n')
endfunction

// Implements the Secant method
function[xk,i] = SM(x0,x1,tol,i)
    mprintf('Secant method\n')
    xkm1 = x0
    fkm1 = f(x0)
    xk = x1
    fk = f(x1)
    while i < 100
        xkp1 = xk - (fk * (xk-xkm1)) / (fk-fkm1)
        mprintf('%i\tx = %.11f\tf(x) = %.11f\n',i,xk,f(xk))
        if abs(xkp1 - xk) < tol then
            break
        end
        xkm1 = xk
        fkm1 = fk
        xk = xkp1
        fk = f(xkp1)
        i = i + 1
    end
endfunction

IBM(a,b,tol,i)
NM(a,tol,i)
SM(a,b,tol,i)
