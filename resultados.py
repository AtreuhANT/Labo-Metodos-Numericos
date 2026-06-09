=== ∫₀³ (1/√2π) e^(-x²/2) dx ===
Método                  | Puntos F(x) | Aproximación   | Error Absoluto
Trapecio Compuesto      |     13      | 0.4985812865   |   6.88e-05
Simpson 1/3 Compuesto   |     13      | 0.4986483856   |   1.72e-06
Simpson 3/8 Compuesto   |     13      | 0.4986462838   |   3.82e-06
Gauss-Legendre (n=5)    |      5      | 0.4986490706   |   1.03e-06
Valor exacto: 0.4986501020

=== ∫₀² cos(x²) dx ===
Método                  | Puntos F(x) | Aproximación   | Error Absoluto
Trapecio Compuesto      |     13      | 0.4685037984   |   7.04e-03
Simpson 1/3 Compuesto   |     13      | 0.4613260847   |   1.35e-04
Simpson 3/8 Compuesto   |     13      | 0.4611780673   |   2.83e-04
Gauss-Legendre (n=5)    |      5      | 0.4612266926   |   2.35e-04
Valor exacto: 0.4614614624

=== ∫₋₂³ (x⁵ - 2x³ + 4) dx ===  ← POLINOMIO de grado 5
Método                  | Puntos F(x) | Aproximación   | Error Absoluto
Trapecio Compuesto      |     13      | 102.588714...  |   4.26e+00
Simpson 1/3 Compuesto   |     13      | 98.383568...   |   5.02e-02
Simpson 3/8 Compuesto   |     13      | 98.446361...   |   1.13e-01
Gauss-Legendre (n=3)    |      3      | 98.333333...   |   4.26e-14  ✔ EXACTO
Valor exacto: 98.3333333333
