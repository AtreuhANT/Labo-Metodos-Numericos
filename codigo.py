#Para ejecutar el codigo escribir en la terminal: pip install numpy sympy

y se instalara la libreria requerida

"""
Ejercicio 1 - Clase de Laboratorio
Integración Numérica: Trapecio, Simpson 1/3, Simpson 3/8 y Gauss-Legendre
"""

import numpy as np
import sympy as sp
import re

# ─────────────────────────────────────────────
# UTILIDADES GENERALES
# ─────────────────────────────────────────────

def parse_function(expr_str):
    """Convierte string de función a callable numpy."""
    x = sp.Symbol('x')
    expr = sp.sympify(expr_str)
    f_sympy = sp.lambdify(x, expr, modules=['numpy'])
    return f_sympy, expr

def is_polynomial(expr_str):
    """Devuelve True y el grado si la expresión es un polinomio en x."""
    x = sp.Symbol('x')
    try:
        expr = sp.sympify(expr_str)
        poly = sp.Poly(expr, x)
        return True, poly.degree()
    except Exception:
        return False, None

# ─────────────────────────────────────────────
# MÉTODOS DE NEWTON-COTES
# ─────────────────────────────────────────────

def trapecio_compuesto(f, a, b, N=12):
    """Regla del Trapecio Compuesta con N subintervalos."""
    h = (b - a) / N
    x = np.linspace(a, b, N + 1)
    y = f(x)
    result = h * (y[0]/2 + np.sum(y[1:-1]) + y[-1]/2)
    puntos = N + 1
    return result, puntos

def simpson_1_3_compuesto(f, a, b, N=12):
    """Regla de Simpson 1/3 Compuesta. N debe ser par."""
    if N % 2 != 0:
        raise ValueError("N debe ser par para Simpson 1/3")
    h = (b - a) / N
    x = np.linspace(a, b, N + 1)
    y = f(x)
    result = h/3 * (y[0] + 4*np.sum(y[1:-1:2]) + 2*np.sum(y[2:-2:2]) + y[-1])
    puntos = N + 1
    return result, puntos

def simpson_3_8_compuesto(f, a, b, N=12):
    """Regla de Simpson 3/8 Compuesta. N debe ser múltiplo de 3."""
    if N % 3 != 0:
        raise ValueError("N debe ser múltiplo de 3 para Simpson 3/8")
    h = (b - a) / N
    x = np.linspace(a, b, N + 1)
    y = f(x)
    S = y[0] + y[-1]
    for i in range(1, N):
        if i % 3 == 0:
            S += 2 * y[i]
        else:
            S += 3 * y[i]
    result = 3*h/8 * S
    puntos = N + 1
    return result, puntos

# ─────────────────────────────────────────────
# GAUSS-LEGENDRE (raíces y pesos calculados algorítmicamente)
# ─────────────────────────────────────────────

def gauss_legendre_nodes_weights(n):
    """
    Calcula algorítmicamente las raíces y pesos de Gauss-Legendre de n puntos
    usando el método de Newton sobre los polinomios de Legendre.
    """
    # Estimación inicial de raíces (Golub-Welsch / aproximación trigonométrica)
    nodes = np.zeros(n)
    weights = np.zeros(n)

    for i in range(1, n // 2 + 1):
        # Estimación inicial
        xi = np.cos(np.pi * (i - 0.25) / (n + 0.5))

        # Iteración de Newton
        for _ in range(100):
            p0, p1 = 1.0, xi
            for k in range(2, n + 1):
                p2 = ((2*k - 1) * xi * p1 - (k - 1) * p0) / k
                p0, p1 = p1, p2
            # p2 = P_n(xi), dp = P_n'(xi)
            dp = n * (p0 - xi * p2) / (1 - xi**2)
            dx = p2 / dp
            xi -= dx
            if abs(dx) < 1e-15:
                break

        nodes[i - 1] = -xi
        nodes[n - i] = xi
        # Peso
        p0, p1 = 1.0, xi
        for k in range(2, n + 1):
            p2 = ((2*k - 1) * xi * p1 - (k - 1) * p0) / k
            p0, p1 = p1, p2
        dp = n * (p0 - xi * p2) / (1 - xi**2)
        weights[i - 1] = 2.0 / ((1 - xi**2) * dp**2)
        weights[n - i] = weights[i - 1]

    # Si n es impar, el nodo central es 0
    if n % 2 == 1:
        nodes[n // 2] = 0.0
        xi = 0.0
        p0, p1 = 1.0, xi
        for k in range(2, n + 1):
            p2 = ((2*k - 1) * xi * p1 - (k - 1) * p0) / k
            p0, p1 = p1, p2
        dp = n * (p0 - xi * p2) / (1 - xi**2) if xi != 1 else n * (-1)**(n+1) * n / 2
        weights[n // 2] = 2.0 / ((1 - xi**2) * dp**2) if xi != 1 else 2.0 / (n * (n + 1))

    return nodes, weights

def gauss_legendre_integrate(f, a, b, n):
    """Integra f en [a,b] con n puntos de Gauss-Legendre."""
    nodes, weights = gauss_legendre_nodes_weights(n)
    # Cambio de variable: t en [-1,1] -> x en [a,b]
    x_mapped = 0.5 * (b - a) * nodes + 0.5 * (b + a)
    result = 0.5 * (b - a) * np.dot(weights, f(x_mapped))
    return result, n  # n es el número de evaluaciones de f

def gauss_legendre_exact_n(grado):
    """
    Devuelve el n mínimo de Gauss-Legendre para integrar exactamente
    un polinomio de grado 'grado'.
    Precisión exacta si 2n-1 >= grado  =>  n >= (grado+1)/2
    """
    return int(np.ceil((grado + 1) / 2))

# ─────────────────────────────────────────────
# INTEGRAL EXACTA (SymPy)
# ─────────────────────────────────────────────

def integral_exacta(expr_str, a, b):
    x = sp.Symbol('x')
    expr = sp.sympify(expr_str)
    val = float(sp.integrate(expr, (x, a, b)))
    return val

# ─────────────────────────────────────────────
# TABLA DE RESULTADOS
# ─────────────────────────────────────────────

def imprimir_tabla(resultados, valor_exacto):
    col_metodo = 24
    col_puntos = 12
    col_aprox  = 16
    col_error  = 16
    sep = "-" * (col_metodo + col_puntos + col_aprox + col_error + 12)

    print(f"\n{'Método':<{col_metodo}} | {'Puntos F(x)':^{col_puntos}} | {'Aproximación':^{col_aprox}} | {'Error Absoluto':^{col_error}}")
    print(sep)
    for nombre, aprox, puntos in resultados:
        error = abs(aprox - valor_exacto)
        print(f"{nombre:<{col_metodo}} | {puntos:^{col_puntos}} | {aprox:^{col_aprox}.10f} | {error:^{col_error}.2e}")
    print(sep)
    print(f"  Valor exacto: {valor_exacto:.10f}")

# ─────────────────────────────────────────────
# MENÚ PRINCIPAL
# ─────────────────────────────────────────────

def main():
    print("=" * 65)
    print("   INTEGRACIÓN NUMÉRICA - Laboratorio de Métodos Numéricos")
    print("=" * 65)
    print("Ingrese la función a integrar en términos de x.")
    print("Ejemplos: (1/sqrt(2*pi))*exp(-x**2/2)  |  cos(x**2)  |  x**5-2*x**3+4")
    print()

    expr_str = input("f(x) = ").strip()
    a = float(input("Límite inferior a = "))
    b = float(input("Límite superior b = "))
    N = 12  # subintervalos fijos

    try:
        f, expr = parse_function(expr_str)
    except Exception as e:
        print(f"Error al parsear la función: {e}")
        return

    # Verificar si es polinomio
    es_poli, grado = is_polynomial(expr_str)
    if es_poli:
        print(f"\n✔ La función es POLINÓMICA de grado {grado}.")
        n_gl_exacto = gauss_legendre_exact_n(grado)
        print(f"  → Gauss-Legendre usará n = {n_gl_exacto} nodo(s) para solución EXACTA.")
        n_gl = n_gl_exacto
    else:
        print("\n✘ La función NO es polinómica. Gauss-Legendre usará n = 5 nodos.")
        n_gl = 5

    # Valor exacto
    try:
        val_exacto = integral_exacta(expr_str, a, b)
        print(f"\n  Integral exacta (SymPy): {val_exacto:.10f}")
    except Exception as e:
        print(f"  No se pudo calcular la integral exacta simbólicamente: {e}")
        val_exacto = None

    # Aplicar métodos
    resultados = []

    aprox_trap, pts_trap = trapecio_compuesto(f, a, b, N)
    resultados.append(("Trapecio Compuesto", aprox_trap, pts_trap))

    # Simpson 1/3 requiere N par; N=12 es par ✔
    aprox_s13, pts_s13 = simpson_1_3_compuesto(f, a, b, N)
    resultados.append(("Simpson 1/3 Compuesto", aprox_s13, pts_s13))

    # Simpson 3/8 requiere N múltiplo de 3; N=12 ✔
    aprox_s38, pts_s38 = simpson_3_8_compuesto(f, a, b, N)
    resultados.append(("Simpson 3/8 Compuesto", aprox_s38, pts_s38))

    aprox_gl, pts_gl = gauss_legendre_integrate(f, a, b, n_gl)
    resultados.append(("Gauss-Legendre", aprox_gl, pts_gl))

    # Mostrar tabla
    if val_exacto is not None:
        imprimir_tabla(resultados, val_exacto)

if __name__ == "__main__":
    main()
