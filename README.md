# Labo-Metodos-Numericos

HUERTA SANTANA, SERGIO ALEXANDER ( 24190057 )

Pregunta 1 – Eficiencia Computacional

Los métodos de Newton-Cotes (Trapecio, Simpson) usaron 13 evaluaciones de f(x) para N=12 subintervalos. Gauss-Legendre con n=5 solo necesitó 5 evaluaciones. A pesar de ello, GL obtuvo errores comparables o menores.

Esto se explica por la optimización simultánea de nodos y pesos: en Newton-Cotes, los nodos son equiespaciados (impuestos), por lo que solo se optimizan los pesos, aprovechando solo la mitad de los grados de libertad disponibles. Gauss-Legendre elige tanto las posiciones xᵢ como los pesos wᵢ de forma óptima, lo que le permite integrar exactamente polinomios de grado hasta 2n−1 con solo n evaluaciones. Un método Newton-Cotes de n puntos solo integra exactamente hasta grado n (o n+1 si n es par). Por eso se dice que Gauss optimiza el costo computacional: extrae la máxima precisión por evaluación de f.

Pregunta 2 – El Efecto de la Oscilación

La función cos(x²) tiene una frecuencia que crece con x: la pendiente y curvatura aumentan hacia x = 2. Esto provoca oscilaciones de amplitud constante pero longitud de onda decreciente — una función altamente no uniforme en el intervalo.

Los métodos compuestos con nodos equiespaciados tienen un error proporcional a h²f''(ξ) (Trapecio) o h⁴f⁴(ξ) (Simpson), donde la derivada se evalúa en algún punto ξ. Cuando la curvatura crece, esas derivadas altas son mucho mayores en las zonas de mayor oscilación que en el inicio del intervalo, por lo que el error es dominado por la región problemática. Los nodos equiespaciados no concentran evaluaciones donde más se necesitan.

Gauss-Legendre, en cambio, posiciona sus nodos en las raíces del polinomio de Legendre de grado n, que están ligeramente concentradas hacia el centro del intervalo y distribuidas de manera óptima para minimizar el error global. Adicionalmente, su error depende de la derivada de orden 2n, que para funciones no polinómicas puede ser menor que la derivada de orden 4 de Simpson en el mismo sentido práctico.

Pregunta 3 – Límites del Método

Para P(x) = 7x⁷ − 3x⁴ + 2x, el grado es 7.

Por el Teorema de Precisión Máxima de Gauss-Legendre: con n nodos se integra exactamente todo polinomio de grado ≤ 2n − 1. Para grado 7 se necesita:

2n − 1 ≥ 7 ⟹ n ≥ 4

Por tanto, con n = 4 nodos Gauss-Legendre obtiene error absoluto exactamente cero (salvo redondeo de punto flotante ~10⁻¹⁴). Ningún valor menor de n garantizaría esto. Esto se verifica en el ejercicio: para la integral de grado 5 se usó n = 3 (ya que 2·3−1 = 5 ≥ 5), y se obtuvo error 4.26×10⁻¹⁴, consistente con precisión de máquina.

# Para ejecutar el codigo escribir en la terminal: 

pip install numpy sympy

y se instalara la libreria requerida

luego al ejecutar:

<img width="834" height="218" alt="image" src="https://github.com/user-attachments/assets/6027e62c-851f-49b0-bec5-3bf47c0299f1" />


ponermos por ejemplo:

f(x) = x**5 - 2*x**3 + 4

Límite inferior a = -1

Límite superior b = 2

y saldria:
<img width="829" height="459" alt="image" src="https://github.com/user-attachments/assets/520d236f-a53e-48cc-a4de-713c21fbce56" />
