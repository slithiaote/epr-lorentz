# Statement

Compute Hessian of the max-likelihood estimator in the Lorentz model.
Use Python Sympy library

# Assets and Artefacts : not applicable

# Method

1. Define the log likelihood of the Lorentz Model
1. Compute derivative with symbolic library
1. Convert expression to matlab code

## First use of the sympy library


```python
import sympy as sym
sym.init_printing();

x = sym.Symbol('x')
sym.diff(x ** 5)

```

```text
   4
5⋅x 
```


## Mathematical representation of the EPR signal

In this document, we focus on single line EPR spectra and present the mathematical model.
EPR spectra be represented mathematically as a mother shape (a Gaussian or Lorentzian function) after translation and dilation. 
Correspondingly, we introduce the dma and scale operators to represent these operations and derive their properties with regard to differentiation.

We first introduce required variable names.


```python
from sympy.abc import e # generic name for an expression

# variables names for EPR model parameters
from sympy.abc import v, C, B
Br = sym.Symbol('Br');
gamma = sym.Symbol("gamma") # FWHM = sym.Symbol('FWHM');
MA = sym.Symbol('MA');

from sympy.abc import n, Y  # variables for the observed EPR

```

```text
```

### Testing function

The following function is used throughout the document to test equality of two expressions, and check that a formula is correct.


```python
def test_diff(a,e1,e2):
    if (e1.expand() - e2.expand() == 0):
        print(a,"ok")
    else:
        print(e1)
        print(e2)

test_diff("test",sym.diff(x**5),5*x**4)

e = sym.Function("e")(B)

```

```text
test ok
```

### The scale operator

The scale operator applies a dilation to a given expression, while preserving overall Area-Under-the-Curve. 
In EPR, this means that the overall intensity or quantity of matter is the same.


```python
def scale(e):
    # apply scaling transform to expression e while preserving area
    return e.subs(B,B/gamma)/gamma

```

```text
```

The scale operator has the following properties with respect to differentiation. The scale operator commutes with
differetiation by v, C, Br and MA because these variables are not involved. We provide two useful formulas when 
differentiating by B and gamma.


```python
scale(e).diff(v) == scale(e.diff(v))
test_diff("formula", scale(e).diff(B), 1/gamma * scale(e.diff(B)))
test_diff("formula", scale(e).diff(gamma),  - 1/gamma * ( scale(e) + scale(B*e.diff(B)) ) )

```

```text
True
formula ok
formula ok
```


### The dma operator

In EPR, the signal is not measured directly as intensity as a function of the magnetic field B. Instead, a finite difference
scheme is used, with parameter MA. We introduce the dma operator to represent both translation and this finite difference scheme.


```python
def dma(e):
    # applies the dma operator
    return e.subs(B,B-Br+MA/2)-e.subs(B,B-Br-MA/2)

```

```text
```

Let us examine the properties of the dma operator with respect to differentiation.
The dma operator involves B, MA and Br, so differentiation with other variables commute.

```python
dma(e).diff(C) == dma(e.diff(C))
dma(e).diff(B) == dma(e.diff(B))
dma(e).diff(Br) == - dma(e.diff(B))

```

```text
True
True
True
```

To differentiate by MA, we introduce another operator.

```python
def dmam(e):
    # companion operator to dma wher differentiating by MA
    return e.subs(B,B-Br+MA/2)/2+e.subs(B,B-Br-MA/2)/2

dma(e).diff(MA) == dmam(e.diff(B))

```

```text
True
```

As we will perform double differentiation, we need to understand how to differentiate the dmam operator. This operator
has similar properties to dma.


```python
dmam(e).diff(Br) == - dmam(e.diff(B))
dmam(e).diff(MA) == dma(e.diff(B))/4

```

```text
True
True
```


### EPR signal and likelihood function

We are now ready to define the model for the EPR signal. Let Abs denote the normalized shape of a peak in the model. 
Then a single line EPR spectrum with AUC C, location Br, scale gamma is modelled by the following expression.
The Modulation Amplitude parameter MA is implicitly defined in the dma operator.


```python
Abs = sym.Function("Abs")(B)
EPR = C * dma(scale(Abs))
EPR

```

```text
  ⎛     ⎛         MA⎞      ⎛         MA⎞⎞
  ⎜     ⎜B - Br - ──⎟      ⎜B - Br + ──⎟⎟
  ⎜     ⎜         2 ⎟      ⎜         2 ⎟⎟
  ⎜  Abs⎜───────────⎟   Abs⎜───────────⎟⎟
  ⎜     ⎝     γ     ⎠      ⎝     γ     ⎠⎟
C⋅⎜- ──────────────── + ────────────────⎟
  ⎝         γ                  γ        ⎠
```

The Abs function represents any model such as the Lorentzian, Gaussian or Voigt lineshape models. The following expressions
are provided in terms of Abs, which makes it possible to obtain generic formulae that are valid for generic location-scale 
families of shapes. The Voigt lineshape model requires more work as it has 2 parameters for its shape.

As an example, this is the Lorentzian model:

```python
from sympy.abc import pi
FWHM = sym.Symbol("FWHM")
# spec = (B.^2+1).^(-1) ./ pi;
#  mol$C / (pi * (mol$fwhm/2) * (1 + ((mac$B - mol$Br)/ (mol$fwhm/2) )^2))
Lorentz = 1 / pi / (1 + B**2)
(C * scale(Lorentz)).subs(gamma, FWHM/2)
(C * dma(scale(Lorentz))).subs(gamma,FWHM/2)

```

```text
       2⋅C        
──────────────────
       ⎛    2    ⎞
       ⎜ 4⋅B     ⎟
FWHM⋅π⋅⎜───── + 1⎟
       ⎜    2    ⎟
       ⎝FWHM     ⎠
  ⎛              2                               2              ⎞
C⋅⎜───────────────────────────── - ─────────────────────────────⎟
  ⎜       ⎛                   2⎞          ⎛                   2⎞⎟
  ⎜       ⎜      ⎛         MA⎞ ⎟          ⎜      ⎛         MA⎞ ⎟⎟
  ⎜       ⎜    4⋅⎜B - Br + ──⎟ ⎟          ⎜    4⋅⎜B - Br - ──⎟ ⎟⎟
  ⎜       ⎜      ⎝         2 ⎠ ⎟          ⎜      ⎝         2 ⎠ ⎟⎟
  ⎜FWHM⋅π⋅⎜1 + ────────────────⎟   FWHM⋅π⋅⎜1 + ────────────────⎟⎟
  ⎜       ⎜             2      ⎟          ⎜             2      ⎟⎟
  ⎝       ⎝         FWHM       ⎠          ⎝         FWHM       ⎠⎠
```

We assume that the observed EPR spectrum Y is obtained in the instrument according to the EPR model with additive white
Gaussian noise (standard deviation v). Consequently, the log-likelihood function ll is defined in the following chunk.
In this expression, B are the input values of magnetic field and Y are the measured spectrum values at B. 
In reality, we measure Y_i = Spectrum(B_i) and the log-likelihood is the sum of the vectorized expression.


```python
ll = - n * sym.log(v) - 1 / (2* v ** 2) * (Y - EPR)**2; # complete definition
ll # complete expression

```

```text
                                                             2
            ⎛    ⎛     ⎛         MA⎞      ⎛         MA⎞⎞    ⎞ 
            ⎜    ⎜     ⎜B - Br - ──⎟      ⎜B - Br + ──⎟⎟    ⎟ 
            ⎜    ⎜     ⎜         2 ⎟      ⎜         2 ⎟⎟    ⎟ 
            ⎜    ⎜  Abs⎜───────────⎟   Abs⎜───────────⎟⎟    ⎟ 
            ⎜    ⎜     ⎝     γ     ⎠      ⎝     γ     ⎠⎟    ⎟ 
            ⎜- C⋅⎜- ──────────────── + ────────────────⎟ + Y⎟ 
            ⎝    ⎝         γ                  γ        ⎠    ⎠ 
-n⋅log(v) - ──────────────────────────────────────────────────
                                      2                       
                                   2⋅v                        
```

## Computing the Hessian matrix

### Aliases

We define the following aliases. These will be defined as functions in the Matlab implementation.


```python
AbsB = Abs.diff(B)
AbsBB = Abs.diff(B).diff(B)
BAbsB = B*AbsB
BAbsBB = B*AbsBB
BBAbsBB = B*B*AbsBB

```

```text
```

In the case of the Lorentz lineshape model:

```python
test_diff("Lorentz", Lorentz.diff(B), -2 / pi * B * (B**2+1)**(-2))
test_diff("Lorentz", Lorentz.diff(B).diff(B), 8/pi * B**2 * (B**2+1)**(-3) - 2/pi * (B**2+1)**(-2))

```

```text
Lorentz ok
Lorentz ok
```

The final expressions only depend on the following two "condensed" operators.


```python
def ds(e):
    # dma(scale(e))
    return e.subs(B,(B-Br+MA/2)/gamma)/gamma - e.subs(B,(B-Br-MA/2)/gamma)/gamma

def dms(e):
    # dmam(scale(e))
    return e.subs(B,(B-Br+MA/2)/gamma)/gamma/2 + e.subs(B,(B-Br-MA/2)/gamma)/gamma/2

```

```text
```

### Differentiate by v

Let's perform differentiation starting with v.

```python
dv = -n / v + 1 / (v ** 3) * (Y - EPR)**2;

dvv =  n / v**2 - 3 / v**4 * (Y - EPR)**2;
dvC = -2 / v**3 / C * (Y - EPR) * EPR
dvB =  2 / gamma * C / v**3 * (Y - EPR) * ds(AbsB)
dvM = -2 / gamma * C / v**3 * (Y - EPR) * dms(AbsB)
dvg =  2 / gamma * C / v**3 * (Y - EPR) * (ds(Abs)+ds(BAbsB))

test_diff("dv", ll.diff(v), dv)
test_diff("dvv", dv.diff(v), dvv)
test_diff("dvC", dv.diff(C), dvC)
test_diff("dvB", dv.diff(Br), dvB)
test_diff("dvM", dv.diff(MA), dvM)
test_diff("dvg", dv.diff(gamma), dvg)

```

```text
dv ok
dvv ok
dvC ok
dvB ok
dvM ok
dvg ok
```

### Differentiate by C


```python
dC = 1 / v ** 2 * (Y - EPR) * ds(Abs)

dCC = -1 / v ** 2 * ds(Abs) ** 2
dCB = -1 / v ** 2 / gamma * (Y - 2*EPR) * ds(AbsB)
dCM =  1 / v ** 2 / gamma * (Y - 2*EPR) * dms(AbsB)
dCg = -1 / v ** 2 / gamma * (Y - 2*EPR) * (ds(Abs)+ds(BAbsB))

test_diff("dC", ll.diff(C), dC)
test_diff("dCv is dvC", dC.diff(v), dvC)
test_diff("dCC", dC.diff(C), dCC)
test_diff("dCB", dC.diff(Br), dCB)
test_diff("dCM", dC.diff(MA), dCM)
test_diff("dCg", dC.diff(gamma), dCg)

```

```text
dC ok
dCv is dvC ok
dCC ok
dCB ok
dCM ok
dCg ok
```

### Differentiate by Br


```python
dB = -C / v ** 2 / gamma * (Y - EPR) * ds(AbsB)
dBB = -C / v**2 / gamma**2 * ( C * ds(AbsB)**2 - (Y-EPR) * ds(AbsBB) )
dBM = -C / v**2 / gamma**2 * ( -C * dms(AbsB) * ds(AbsB) + (Y-EPR) * dms(AbsBB) )
dBg = -C / v**2 / gamma**2 * ( (3*EPR-2*Y + C*ds(BAbsB)) * ds(AbsB) - (Y-EPR) * ds(BAbsBB)  )

test_diff("dB", ll.diff(Br), dB)
test_diff("dBv is dvB", dB.diff(v), dvB)
test_diff("dBC is dCB", dB.diff(C), dCB)
test_diff("dBB", dB.diff(Br), dBB)
test_diff("dBM", dB.diff(MA), dBM)
test_diff("dBg", dB.diff(gamma), dBg)

```

```text
dB ok
dBv is dvB ok
dBC is dCB ok
dBB ok
dBM ok
dBg ok
```

### Differentiate by MA


```python
dM = C / v**2 / gamma * (Y-EPR) * dms(AbsB)
dMM = C / v**2 / gamma**2 * (-C * dms(AbsB)**2 + (Y-EPR) * ds(AbsBB)/4 )
dMg = C / v**2 / gamma**2 * ( C * (ds(Abs)+ds(BAbsB)) * dms(AbsB) - (Y-EPR) * (2*dms(AbsB) + dms(BAbsBB)) )

test_diff("dM",ll.diff(MA),dM)
test_diff("dMv is dvM", dM.diff(v),dvM)
test_diff("dMC is dCM", dM.diff(C),dCM)
test_diff("dMB is dBM", dM.diff(Br),dBM)
test_diff("dMM", dM.diff(MA), dMM)
test_diff("dMg", dM.diff(gamma), dMg)

```

```text
dM ok
dMv is dvM ok
dMC is dCM ok
dMB is dBM ok
dMM ok
dMg ok
```

### Differentiate by gamma


```python
dg = C / v**2 * (Y-EPR)*dma(scale(Abs).diff(gamma))
dgg = C / v**2 / gamma**2 * ( -C*(ds(Abs)+ds(BAbsB))**2 + (Y-EPR)*(2*ds(Abs) + 4*ds(BAbsB) + ds(BBAbsBB)))

test_diff("dg",ll.diff(gamma),dg)
test_diff("dgv is dvg", dg.diff(v),dvg)
test_diff("dgC is dCg", dg.diff(C),dCg)
test_diff("dgB is dBg", dg.diff(Br),dBg)
test_diff("dgM is dMg", dg.diff(MA),dMg)
test_diff("dgg", dg.diff(gamma), dgg)

```

```text
dg ok
dgv is dvg ok
dgC is dCg ok
dgB is dBg ok
dgM is dMg ok
dgg ok
```




## Application to the Lorentz lineshape model

The Lorentz lineshape model corresponds to 

```python
Lorentz

```

```text
    1     
──────────
  ⎛ 2    ⎞
π⋅⎝B  + 1⎠
```

We compute here using sympy the expressions for the aliases in the case of the Lorentz lineshape model.

```python
Lorentz.diff(B)
Lorentz.diff(B).diff(B)

```

```text
   -2⋅B    
───────────
          2
  ⎛ 2    ⎞ 
π⋅⎝B  + 1⎠ 
       2                 
    8⋅B            2     
─────────── - ───────────
          3             2
  ⎛ 2    ⎞      ⎛ 2    ⎞ 
π⋅⎝B  + 1⎠    π⋅⎝B  + 1⎠ 
```

We apply / prepare a check for matlab.


```python
AbsB = Lorentz.diff(B); AbsB.subs(B,3)
AbsBB = Lorentz.diff(B).diff(B); AbsBB.subs(B,3)
ds(AbsB).subs([(B,3),(Br,10),(MA,2),(gamma,6)])
dms(AbsB).subs([(B,3),(Br,10),(MA,2),(gamma,6)])

```

```text
-3  
────
50⋅π
  13 
─────
250⋅π
 193  
──────
7500⋅π
  1057 
───────
15000⋅π
```


```python
AbsB = Lorentz.diff(B); AbsB.subs(B,5)
AbsBB = Lorentz.diff(B).diff(B); AbsBB.subs(B,5)
ds(AbsB).subs([(B,5),(Br,10),(MA,2),(gamma,6)])
dms(AbsB).subs([(B,5),(Br,10),(MA,2),(gamma,6)])

```

```text
 -5  
─────
338⋅π
  37  
──────
4394⋅π
  47  
──────
2028⋅π
 385  
──────
4056⋅π
```





