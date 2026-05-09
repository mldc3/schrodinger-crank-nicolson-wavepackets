import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags, csr_matrix, lil_matrix
from scipy.sparse.linalg import spsolve
from scipy.linalg import expm
import time
import warnings
import matplotlib.animation as animation
import os
from scipy.linalg import solve_banded


########################### APARTADO 1: THOMAS PARA V=0 (esta sección saltatela) ###########################
#en este caso simplemente lo aplicamos a V=0 pero se puede generalizar definiendo otro CN_step
def thomas(a, b, c, d):
    n = len(b)                                  # Se obtiene el tamaño del sistema (n ecuaciones)
    ac, bc, cc, dc = map(np.copy, (a, b, c, d)) # definimos la subdiagonal, diagonal principal, superdiagonal y un vector de términos independientes

    # Se crean copias de los vectores de entrada para no modificar los originales
    for i in range(1, n):
        m = ac[i-1] / bc[i-1]                   # Factor multiplicador: m = a_{i-1} / b_{i-1} Este valor se usa para eliminar el término a_i
        bc[i] -= m * cc[i-1]                    # Se actualiza la diagonal principal:  b_i = b_i - m * c_{i-1}
        dc[i] -= m * dc[i-1]                    # Se actualiza el vector de términos independientes: d_i = d_i - m * d_{i-1}

    # Back substitution
    x = np.zeros(n, dtype=complex)             
    x[-1] = dc[-1] / bc[-1]                     # Se calcula el último valor de la solución: x_{n-1} = d_{n-1} / b_{n-1}
    for i in range(n-2, -1, -1):                # Se calculan los valores restantes desde atrás hacia adelante
        x[i] = (dc[i] - cc[i] * x[i+1]) / bc[i] # Fórmula de sustitución hacia atrás: x_i = (d_i - c_i * x_{i+1}) / b_i
    return x


def CN_step(psi, dx, dt):
    N = len(psi)
    psi_int = psi[1:-1]         
    n = len(psi_int)

    r = 1j * dt / (2 * dx * dx)

    # Matriz A (izquierda)
    mainA = (1 + 2*r) * np.ones(n, dtype=complex)
    offA  = -r * np.ones(n-1, dtype=complex)

    # Matriz B (derecha)
    mainB = (1 - 2*r) * np.ones(n, dtype=complex)
    offB  = r * np.ones(n-1, dtype=complex)

    # b = B * psi^n
    b = mainB * psi_int
    b[1:] += offB * psi_int[:-1]
    b[:-1] += offB * psi_int[1:]

    # Resolver A psi^{n+1} = b
    new_int = thomas(offA, mainA, offA, b)

    psi_new = np.zeros_like(psi, dtype=complex)
    psi_new[1:-1] = new_int
    return psi_new

"""
# Parámetros
L = 100
N = 800
x = np.linspace(-L, L, N)
dx = x[1] - x[0]
dt = 0.01
steps = 3000
# Paquete gaussiano inicial
x0 = -40
k0 = 1.0
sigma = 6.0
psi = np.exp(-(x-x0)**2/(2*sigma**2)) * np.exp(1j*k0*x)
psi /= np.sqrt(np.trapz(np.abs(psi)**2, x))             # Normalizar

fig, ax = plt.subplots()                                # Creamos la figura y los ejes
line, = ax.plot(x, np.abs(psi)**2)                      # Graficamos la densidad de probabilidad |ψ|²
ax.set_ylim(0, np.max(np.abs(psi)**2) * 1.3)            # Definimos límites del eje y con margen
ax.set_xlabel("x")                                      # Etiquetamos del eje x
ax.set_ylabel("|psi|^2")                                # Etiquetamos del eje y
ax.set_title("Propagación libre (V=0) con Thomas")      # Título de la gráfica

frames = []                                             # Lista para guardar los frames
skip = 10                                               # Guardamos un frame cada 10 pasos

for n in range(steps):
    psi = CN_step(psi, dx, dt)
    if n % skip == 0:                                   # Verificamos si toca guardar frame
        frames.append(np.abs(psi)**2)
        
def update(i):                                          # Función de actualización de la animación
    line.set_ydata(frames[i])                           # Actualizamos la curva con el frame i
    return line,                                        # Devuelvemos el objeto gráfico

ani = animation.FuncAnimation(fig, update, frames=len(frames), blit=True)
ani.save("schrodinger_libre.gif", writer="pillow", fps=20)

print("GIF generado como schrodinger_libre.gif")
"""




########################### APARTADO 2: CONDICIONES CONTORNO CUALQUIER POTENCIAL ###########################
#hacemos una función que genere los coeficientes de thomas para que despues la de thomas lo resuelva
def solveThomas(A, b):
    N = A.shape[0]                               # Tamaño del sistema

    # Extraemos diagonales de A
    a = np.diag(A, k=-1)                         # Subdiagonal
    d = np.diag(A, k=0)                          # Diagonal principal
    c = np.diag(A, k=1)                          # Superdiagonal

    # Resolvemos usando el método de Thomas
    return thomas(a, d, c, b)

#Ahora generamos una funcion que aplique CN independientemente del potencial
def CN_step(psi, V, dx, dt, resolver="normal", cc_tipo="dirichlet"):
    N = len(psi)
    r = 1j * dt / (2 * dx * dx)

    # CASO DIRICHLET (paredes duras)
    if cc_tipo == "dirichlet":

        # solo puntos interiores
        psi_i = psi[1:-1]
        V_i   = V[1:-1]
        M = N - 2

        # diagonales
        mainA = (1 + 2*r + 1j * dt/2 * V_i).astype(complex)
        offA  = -r * np.ones(M-1, dtype=complex)

        mainB = (1 - 2*r - 1j * dt/2 * V_i).astype(complex)
        offB  =  r * np.ones(M-1, dtype=complex)

        # matrices tridiagonales
        A = np.zeros((M, M), dtype=complex)
        B = np.zeros((M, M), dtype=complex)

        np.fill_diagonal(A, mainA)
        np.fill_diagonal(B, mainB)

        for i in range(M-1):
            A[i, i+1] = offA[i]
            A[i+1, i] = offA[i]
            B[i, i+1] = offB[i]
            B[i+1, i] = offB[i]

        # lado derecho
        b = B @ psi_i

        # resolver
        if resolver == "normal":
            psi_new_i = np.linalg.solve(A, b)
        elif resolver == "thomas":
            psi_new_i = solveThomas(A, b)

        # reconstruir solución completa
        psi_new = np.zeros_like(psi)
        psi_new[1:-1] = psi_new_i
        return psi_new


    # CASO PERIÓDICO
    elif cc_tipo == "periodic":

        mainA = (1 + 2*r + 1j * dt/2 * V).astype(complex)
        offA  = -r * np.ones(N-1, dtype=complex)
        mainB = (1 - 2*r - 1j * dt/2 * V).astype(complex)
        offB  =  r * np.ones(N-1, dtype=complex)

        A = np.zeros((N, N), dtype=complex)
        B = np.zeros((N, N), dtype=complex)

        np.fill_diagonal(A, mainA)
        np.fill_diagonal(B, mainB)

        for i in range(N-1):
            A[i,   i+1] = offA[i]
            A[i+1, i  ] = offA[i]
            B[i,   i+1] = offB[i]
            B[i+1, i  ] = offB[i]

        # acoplos periódicos
        A[0,  -1] = -r
        A[-1,  0] = -r
        B[0,  -1] =  r
        B[-1,  0] =  r

        b = B @ psi

        # periodic solo ponemos normal porque thomas no se puede aplicar 
        psi_new = np.linalg.solve(A, b)
        return psi_new


# Parámetros espaciales 
L = 10                      # Longitud total
N = 800                     # Número de puntos
x = np.linspace(-L, L, N)   # Malla espacial
dx = x[1] - x[0]

# Parámetros temporales 
dt = 0.1
steps = 2000   #subir a unoa 1500

 
# Hemos intentando aplicar los sigueintes potenciales para comprobar si funcionaba.
#Ejemplo: pozo armónico
"""
def V_func(x): #0.02*x**2  
    return 0.02 * x**2  
"""

def V_func(x): #0
    return np.zeros_like(x)

"""
def V_func(x): #V0*(x**2-a**2)**2
    V0 = 0.0005   # escala del potencial
    a = 15.0      # separación entre los pozos
    return V0 * (x**2 - a**2)**2
"""
"""
def V_func(x): #pozofinito
    V0 = 5.0     # profundidad del pozo
    a  = 8.0     # semiancho del pozo
    return np.where(np.abs(x) < a, -V0, 0.0)
"""
"""
def V_func(x):    #barrera
    V0 = 2.0      # altura de la barrera
    a = 5.0       # semiancho
    return np.where(np.abs(x) < a, V0, 0.0)
"""

V = V_func(x)

# Estado inicial: paquete gaussiano 
x0 = 5
k0 = 4 
sigma = 1.7 
psi0 = np.exp(-(x-x0)**2/(2*sigma**2)) * np.exp(1j * k0 * x)

# Normalización
psi0 /= np.sqrt(np.trapz(np.abs(psi0)**2, x))
"""
#definimos una función que grafica la función a diferentes instantes de tiempo para estudiar su comportamiento
def graficar_instantes(psi0, V, dx, dt, steps, cc_tipo="dirichlet", resolver="normal", N_instantes=10):
    N = len(psi0)
    psi = psi0.copy()                      # Copia del estado inicial
    instantes = np.linspace(0, steps, N_instantes, dtype=int)  # índices de pasos donde graficar
    frames = []                            # Guardar |\phi|^2 en esos instantes

    # Evolucionamos la onda
    for n in range(steps+1):
        if n in instantes:
            frames.append(np.abs(psi)**2)
            print("imagenes de instantes;", n, "de", steps)
        if n < steps:  # para no hacer un paso extra
            print("imagenes de instantes;", n, "de", steps)
            psi = CN_step(psi, V, dx, dt, resolver=resolver, cc_tipo=cc_tipo)

    # Graficar
    fig, ax = plt.subplots(figsize=(8,4))
    colores = plt.cm.viridis(np.linspace(0,1,N_instantes))  # Colores para cada instante
    for i, frame in enumerate(frames):
        ax.plot(x, frame, color=colores[i], lw=2, label=f"t = {instantes[i]*dt:.2f}")
    ax.plot(x, V*0.5/np.max(V+1e-12), 'r--', lw=1, label="V(x) (escalado)")
    ax.set_xlabel("x")
    ax.set_ylabel("|phi|^2")
    ax.set_title(f"Densidad de probabilidad ({cc_tipo}, {resolver})")
    ax.legend()
    plt.show()
    
# Dirichlet + normal
graficar_instantes(psi0, V, dx, dt, steps, cc_tipo="dirichlet", resolver="normal", N_instantes=10)

# Dirichlet + Thomas
graficar_instantes(psi0, V, dx, dt, steps, cc_tipo="dirichlet", resolver="thomas", N_instantes=10)

# Periódicas + normal
graficar_instantes(psi0, V, dx, dt, steps, cc_tipo="periodic", resolver="normal",N_instantes=10)

# Periódicas + Thomas: no lo añado porque da mal porque no es tridiagonal ahora 


#lo comento porque no es lo que te quiero enseñar
def mapa_calor(psi0, V, dx, dt, steps, cc_tipo="dirichlet", resolver="normal", skip=10):
    N = len(psi0)
    psi = psi0.copy()
    frames = []  # Guardamos densidad

    for n in range(steps):
        psi = CN_step(psi, V, dx, dt, resolver=resolver, cc_tipo=cc_tipo)
        if n % skip == 0:
            frames.append(np.abs(psi)**2)

    frames = np.array(frames)           # Convertimos a array 2D: time x space
    t_vals = np.arange(0, steps, skip)*dt
    fig, ax = plt.subplots(figsize=(8,4))
    c = ax.imshow(frames, extent=[x[0], x[-1], t_vals[-1], t_vals[0]],aspect='auto', origin='upper', cmap='viridis')
    ax.set_xlabel("x")
    ax.set_ylabel("t")
    ax.set_title(f"Mapa de calor |phi(x,t)|² ({cc_tipo}, {resolver})")
    fig.colorbar(c, ax=ax, label="|phi|²")
    plt.show()

# Dirichlet + normal
mapa_calor(psi0, V, dx, dt, steps, cc_tipo="dirichlet", resolver="normal")

# Dirichlet + Thomas
mapa_calor(psi0, V, dx, dt, steps, cc_tipo="dirichlet", resolver="thomas")

# Periódicas + normal
mapa_calor(psi0, V, dx, dt, steps, cc_tipo="periodic", resolver="normal")

# Periódicas + Thomas
mapa_calor(psi0, V, dx, dt, steps, cc_tipo="periodic", resolver="thomas")



# Función para crear GIF comparativo
def generar_gif(cc_tipo="dirichlet", nombre_gif="schrodinger.gif"):
    # Copias independientes
    psi_normal = psi0.copy()
    psi_thomas = psi0.copy()

    fig, ax = plt.subplots()
    # Graficamos líneas iniciales
    line_normal, = ax.plot(x, np.abs(psi_normal)**2, lw=2, label="CN + solve")
    line_thomas, = ax.plot(x, np.abs(psi_thomas)**2, lw=2, ls="--", label="CN + Thomas")

    pot_scale = np.max(np.abs(psi_normal)**2)/np.max(V+1e-12)*0.5
    ax.plot(x, V*pot_scale, 'r--', lw=2, label="Potencial V(x) (escalado)")

    ax.set_ylim(0, np.max(np.abs(psi_normal)**2)*1.5)
    ax.set_xlabel("x"); ax.set_ylabel("|ψ|²")
    ax.set_title(f"CN: solver normal vs Thomas ({cc_tipo})")
    ax.legend()

    frames_normal = []
    frames_thomas = []
    skip = 2

    # Evolución temporal
    for n in range(steps):
        psi_normal = CN_step(psi_normal, V, dx, dt, resolver="normal", cc_tipo=cc_tipo)
        psi_thomas = CN_step(psi_thomas, V, dx, dt, resolver="thomas", cc_tipo=cc_tipo)
        if n % skip == 0:
            print(f"{n} de {steps} ({cc_tipo})")
            frames_normal.append(np.abs(psi_normal)**2)
            frames_thomas.append(np.abs(psi_thomas)**2)

    # Función de actualización para la animación
    def update(i):
        line_normal.set_ydata(frames_normal[i])
        line_thomas.set_ydata(frames_thomas[i])
        return line_normal, line_thomas

    ani = animation.FuncAnimation(fig, update, frames=len(frames_normal), interval=30, blit=True)
    ani.save(nombre_gif, writer="pillow", fps=20)
    print(f"GIF generado: {nombre_gif}")

generar_gif(cc_tipo="dirichlet", nombre_gif="schrodinger_dirichlet.gif")
generar_gif(cc_tipo="periodic", nombre_gif="schrodinger_periodic.gif")

"""















########################### APARTADO 3: CONSERVACIÓN DE LA NORMA ###########################
"""
def comprobar_norma(psi0, V, dx, dt, steps):
    # Diccionario para almacenar las normas
    normas = {}

    # Lista de casos a estudiar (solo combinaciones válidas)
    casos = [
        ("Dirichlet + normal",  "dirichlet", "normal",  "-"),
        ("Dirichlet + Thomas",  "dirichlet", "thomas",  "--"),
        ("Periódicas + normal", "periodic",  "normal",  ":")
    ]

    # Bucle sobre cada caso
    for nombre, cc_tipo, resolver, estilo in casos:

        # Copia del estado inicial
        psi = psi0.copy()

        # Lista para guardar la norma en cada paso temporal
        norms = []

        # Evolución temporal
        for n in range(steps):
            if n%100==0:
                print(casos, ":", n," de ", steps)
            # Calculamos la norma en el instante actual
            norm = np.sum(np.abs(psi)**2) * dx
            norms.append(norm)

            # Avanzamos un paso temporal con CN
            psi = CN_step(psi, V, dx, dt,
                          resolver=resolver,
                          cc_tipo=cc_tipo)

        # Guardamos la norma y el estilo de línea
        normas[nombre] = (np.array(norms), estilo)

    # Array de tiempos
    t_array = np.arange(steps) * dt


    fig, ax = plt.subplots(figsize=(7,4))
    for nombre, (norms, estilo) in normas.items():
        ax.plot(t_array, norms, estilo, lw=2, label=nombre)
    
    ax.ticklabel_format(axis='y', style='plain', useOffset=False)
    ax.set_xlabel("Tiempo t")
    ax.set_ylabel(r"$N(t) = \int |\psi(x,t)|^2 dx$")
    ax.set_title("Conservación de la norma")
    ax.legend(); ax.grid(True); fig.tight_layout()
    plt.show()

    print("Errores relativos máximos:")
    for nombre, (norms, _) in normas.items():
        err = np.max(np.abs(norms - norms[0]) / norms[0])
        print(f"{nombre}: {err*100:.6f} %")

comprobar_norma(psi0, V, dx, dt, steps)
"""



















########################### APARTADO 4: COMPARACIÓN CON SOLUCIONES ANALÍTICAS ###########################
######### APARTADO 4.1
"""
def potencial_pozo_infinito(x, L_pozo, V0=1e6):
    return np.where(np.abs(x) < L_pozo/2, 0.0, V0)

def solucion_analitica_pozo(x, n, L_pozo, t, m=0.5):
    x_centrado = x + L_pozo/2
    psi = np.zeros_like(x, dtype=complex)
    dentro = (x_centrado >= 0) & (x_centrado <= L_pozo)
    psi[dentro] = np.sqrt(2.0 / L_pozo) * np.sin(n * np.pi * x_centrado[dentro] / L_pozo)
    E_n = (n**2) * (np.pi**2) / (L_pozo**2)       # Con m=1/2
    psi *= np.exp(-1j * E_n * t)
    return psi


L_pozo = 30.0
L_dom  = 50.0
N_pozo = 800
x_pozo = np.linspace(-L_dom/2, L_dom/2, N_pozo)
dx_pozo = x_pozo[1] - x_pozo[0]

dt_pozo    = 0.001
steps_pozo = 5000
n_nivel    = 1

# Estado inicial analítico del pozo infinito
V_pozo = potencial_pozo_infinito(x_pozo, L_pozo)
psi_analitica_0 = solucion_analitica_pozo(x_pozo, n_nivel, L_pozo, t=0.0)
psi_analitica_0 /= np.sqrt(np.trapz(np.abs(psi_analitica_0)**2, x_pozo))

# Copia para la evolución numérica
psi_num = psi_analitica_0.copy()

# Listas para almacenar resultados
frames_pozo_num   = []
frames_pozo_exac  = []
errores_absolutos = []
tiempos_pozo      = []
skip_pozo = 10

for step in range(steps_pozo):
    t_actual = step * dt_pozo
    psi_num = CN_step(psi_num, V_pozo, dx_pozo, dt_pozo, resolver="thomas")
    psi_exac = solucion_analitica_pozo(x_pozo, n_nivel, L_pozo, t_actual)
    psi_exac /= np.sqrt(np.trapz(np.abs(psi_exac)**2, x_pozo))
    error_abs = np.max(np.abs(psi_num - psi_exac))
    errores_absolutos.append(error_abs)
    tiempos_pozo.append(t_actual)
    if step % skip_pozo == 0:
        frames_pozo_num.append(np.abs(psi_num)**2)
        frames_pozo_exac.append(np.abs(psi_exac)**2)

# Animación pozo infinito
fig_comp, ax_comp = plt.subplots()
line_num, = ax_comp.plot(x_pozo, np.abs(psi_analitica_0)**2, 'b-', lw=2, label='Numérica')
line_exa, = ax_comp.plot(x_pozo, np.abs(psi_analitica_0)**2, 'g--', lw=2, label='Analítica')
escala_pot = np.max(np.abs(psi_analitica_0)**2) / (np.max(V_pozo)+1e-12) * 0.5
ax_comp.plot(x_pozo, V_pozo*escala_pot, 'r-.', lw=2, label='V(x) escalado')
ax_comp.set_ylim(0, np.max(np.abs(psi_analitica_0)**2)*1.5)
ax_comp.set_xlabel("x")
ax_comp.set_ylabel("|ψ(x,t)|²")
ax_comp.set_title("Pozo infinito: CN vs solución analítica")
ax_comp.legend()
texto_info = ax_comp.text(0.02, 0.95, "", transform=ax_comp.transAxes)

def actualizar_pozo(i):
    line_num.set_ydata(frames_pozo_num[i])
    line_exa.set_ydata(frames_pozo_exac[i])
    texto_info.set_text(f"Frame {i+1}/{len(frames_pozo_num)}")
    return line_num, line_exa, texto_info

ani_pozo = animation.FuncAnimation(fig_comp, actualizar_pozo, frames=len(frames_pozo_num), interval=30, blit=True)
ani_pozo.save("pozo_infinito_CN_vs_analitica.gif", writer="pillow", fps=20)

plt.figure()
plt.plot(tiempos_pozo, errores_absolutos)
plt.xlabel("t")
plt.ylabel("Error absoluto |ψ_num - ψ_exac|")
plt.title("Error máximo en el pozo infinito")
plt.grid(True)
plt.show()
"""

######### APARTADO 4.2: no me ha salido
"""
print("Comienzo oscilador")
# Potencial armónico: V(x) = 1/2 * \omega^2 * x^2
def potencial_armonico(x, omega):
    return 0.5 * (omega**2) * x**2

# Paquete gaussiano/coherente aproximado en el oscilador armónico
# (forma gaussiana centrada en x0 con impulso p0, width sigma0 ≈ ancho del gs)
def psi_armonico_gaussiano(x, x0, p0, sigma0):
    psi = np.exp(-(x - x0)**2 / (4 * sigma0**2)) * np.exp(1j * p0 * x)
    psi /= np.sqrt(np.trapz(np.abs(psi)**2, x))
    return psi

# Evolución analítica aproximada del centro para el oscilador armónico:
# <x>(t) = x0 cos(ω t) + (p0/çomega) sin(\pmega t)   (Ehrenfest / coherente) 
def x_centro_analitico(t, x0, p0, omega):
    return x0 * np.cos(omega * t) + (p0/omega) * np.sin(omega * t)

# Energía analítica aproximada (invariante): H = p0^2/2 + 1/2 * omega^2 * x0^2
# (ignorando la energía de punto cero)
def energia_analitica(x0, p0, omega):
    return 0.5 * p0**2 + 0.5 * (omega**2) * x0**2


# Parámetros espaciales y del paquete
L_ho   = 10.0            # tamaño del dominio
N_ho   = 2000
x_ho   = np.linspace(-L_ho, L_ho, N_ho)
dx_ho  = x_ho[1] - x_ho[0]

omega  = 1.0             # frecuencia del oscilador
V_ho   = potencial_armonico(x_ho, omega)

x0_ho    = 2.0           # centro inicial del paquete
p0_ho    = 0.0           # impulso inicial
sigma0_ho = 1.0          # ancho inicial (elige cercano al gs si quieres paquete coherente)

psi_ho_0 = psi_armonico_gaussiano(x_ho, x0_ho, p0_ho, sigma0_ho)

# Copia para evolución numérica
psi_ho = psi_ho_0.copy()


# Parámetros temporales
dt_ho    = 0.001
steps_ho = 3000          # t_max = steps_ho * dt_ho

# Almacenamos centro y energía
t_list        = []
x_num_list    = []
E_num_list    = []
x_ana_list    = []
E_ana_const   = energia_analitica(x0_ho, p0_ho, omega)

# Para animación de |psi|^2

# Bucle temporal

frames_ho_num  = []
frames_ho_ana  = []
frames_ho_pot  = []
skip_ho        = 20

for n in range(steps_ho+1):
    t = n * dt_ho

    # Densidad y valores esperados numéricos
    dens = np.abs(psi_ho)**2
    x_num = np.trapz(x_ho * dens, x_ho)

    dpsi_dx = np.gradient(psi_ho, dx_ho)
    T_num = np.trapz(np.abs(dpsi_dx)**2, x_ho)
    V_num = np.trapz(V_ho * dens, x_ho)
    E_num = T_num + V_num

    # Centro analítico (gaussiano coherente)
    x_ana = x_centro_analitico(t, x0_ho, p0_ho, omega)

    t_list.append(t)
    x_num_list.append(x_num)
    E_num_list.append(E_num)
    x_ana_list.append(x_ana)
    print(n, " de ", steps_ho)

    # Estado analítico aproximado en x (para el gif)
    psi_ana_t = psi_armonico_gaussiano(x_ho, x_ana, p0_ho*np.cos(omega*t), sigma0_ho)
    dens_ana  = np.abs(psi_ana_t)**2

    if n % skip_ho == 0:
        frames_ho_num.append(dens.copy())
        frames_ho_ana.append(dens_ana.copy())
        frames_ho_pot.append(V_ho.copy())

    # Evolución CN (no en el último paso)
    if n < steps_ho:
        psi_ho = CN_step(psi_ho, V_ho, dx_ho, dt_ho,
                         resolver="thomas", cc_tipo="dirichlet")
        psi_ho /= np.sqrt(np.trapz(np.abs(psi_ho)**2, x_ho))


fig_ho, ax_ho = plt.subplots(figsize=(10,4))
line_num, = ax_ho.plot(x_ho, frames_ho_num[0], 'b-', lw=2, label='|ψ|² numérico')
line_ana, = ax_ho.plot(x_ho, frames_ho_ana[0], 'g--', lw=2, label='|ψ|² analítico')
escala_pot_ho = np.max(frames_ho_num[0]) / (np.max(V_ho)+1e-12) * 0.5
line_V,  = ax_ho.plot(x_ho, V_ho*escala_pot_ho, 'r-.', lw=1, label='V(x) escalado')

ax_ho.set_xlabel("x")
ax_ho.set_ylabel("|ψ|²")
ax_ho.set_title("Paquete gaussiano en potencial armónico")
ax_ho.legend()
txt_ho = ax_ho.text(0.02, 0.9, "", transform=ax_ho.transAxes)

def actualizar_ho(i):
    line_num.set_ydata(frames_ho_num[i])
    line_ana.set_ydata(frames_ho_ana[i])
    txt_ho.set_text(f"t = {t_list[i*skip_ho]:.3f}")
    return line_num, line_ana, line_V, txt_ho

ani_ho = animation.FuncAnimation(fig_ho, actualizar_ho,frames=len(frames_ho_num),interval=30, blit=True)
ani_ho.save("paquete_armonico_CN_vs_analitica.gif", writer="pillow", fps=20)

plt.show()
"""





















########################### APARTADO 5: ESTUDIO DE CONVERGENCIA para ej V=0 ###########################
"""
def estudio_convergencia_CN(V_func, x0, k0, sigma, T_final,dx_list=[0.5, 0.25, 0.15, 0.10],dt_list=[0.02, 0.015, 0.010, 0.005]):

    print("Estudio de convergencia Crank–Nicolson\n")
    resultados = []

    L = 50.0  # Semilongitud del dominio espacial (dominio = [-L, L])

    # Casos que vamos a estudiar
    casos = [("normal", "dirichlet"),("thomas", "dirichlet"),("normal", "periodic")]

    # Recorremos cada combinación de solver y CC
    for resolver, cc_tipo in casos:
        print(f"\nSolver = {resolver}, CC = {cc_tipo}")

        # Recorremos cada dx objetivo de la lista
        for dx_target in dx_list:
            # Número de puntos espaciales para aproximar ese dx
            N = int(2 * L / dx_target)
            # Malla espacial equiespaciada
            x = np.linspace(-L, L, N)
            # dx real resultante de la malla
            dx_real = x[1] - x[0]

            # Potencial sobre la malla (puede ser libre u otro)
            V = V_func(x)

            # Recorremos cada paso temporal dt
            for dt in dt_list:
                # Número de pasos de tiempo para llegar hasta T_final
                n_steps = int(T_final / dt)

                # Copia del estado inicial para evolucionarlo
                psi = psi0.copy()

                # Norma inicial (debería mantenerse si el esquema es unitario)
                norm_ini = np.trapz(np.abs(psi)**2, x)

                # Medimos tiempo de CPU para este par (dx, dt)
                t0 = time.perf_counter()
                for step in range(n_steps):
                    psi = CN_step(psi, V, dx_real, dt,resolver=resolver, cc_tipo=cc_tipo)
                t1 = time.perf_counter()

                # Norma final tras la evolución
                norm_fin = np.trapz(np.abs(psi)**2, x)
                # Error absoluto en la norma
                error = abs(norm_fin - norm_ini)
                # Tiempo que tarda
                tiempo = t1 - t0

                # Almacenamos resultados en una lista de diccionarios
                resultados.append({"dx_target": dx_target,"dx_real": dx_real,"dt": dt,"N": N,"solver": resolver,"cc_tipo": cc_tipo,"error_norma": error,"tiempo": tiempo})

                # Imprimimos resumen de este experimento
                print(f" dx≈{dx_real:.4f}, dt={dt:.4f}, "f"Error={error:.2e}, tiempo={tiempo:.2f}s")
    return resultados



def graficar_error(resultados, dx_list):
    for (solver, cc) in {
        (r["solver"], r["cc_tipo"]) for r in resultados
    }:
        plt.figure(figsize=(8,6))

        for dx in dx_list:
            dts = [r["dt"] for r in resultados
                   if r["solver"]==solver and r["cc_tipo"]==cc
                   and abs(r["dx_target"]-dx)<1e-12]

            errs = [r["error_norma"] for r in resultados
                    if r["solver"]==solver and r["cc_tipo"]==cc
                    and abs(r["dx_target"]-dx)<1e-12]

            if dts:
                plt.loglog(dts, errs, "o-", label=f"dx={dx}")

        plt.xlabel("dt")
        plt.ylabel("Error absoluto en la norma")
        plt.title(f"Convergencia CN — {solver}, {cc}")
        plt.grid(True, which="both", ls="--")
        plt.legend()
        plt.tight_layout()
        plt.show()


def graficar_tiempo(resultados, dx_list):
    for (solver, cc) in {
        (r["solver"], r["cc_tipo"]) for r in resultados
    }:
        plt.figure(figsize=(8,6))

        for dx in dx_list:
            dts = [r["dt"] for r in resultados
                   if r["solver"]==solver and r["cc_tipo"]==cc
                   and abs(r["dx_target"]-dx)<1e-12]

            times = [r["tiempo"] for r in resultados
                     if r["solver"]==solver and r["cc_tipo"]==cc
                     and abs(r["dx_target"]-dx)<1e-12]

            if dts:
                plt.loglog(dts, times, "o-", label=f"dx={dx}")

        plt.xlabel("dt")
        plt.ylabel("Tiempo de computación (s)")
        plt.title(f"Coste computacional CN — {solver}, {cc}")
        plt.grid(True, which="both", ls="--")
        plt.legend()
        plt.tight_layout()
        plt.show()


def V_libre(x):
    return np.zeros_like(x)

resultados = estudio_convergencia_CN(V_libre,x0=-20.0,k0=1.0,sigma=4.0,T_final=5.0)

graficar_error(resultados, dx_list=[0.5,0.25,0.15,0.10])
graficar_tiempo(resultados, dx_list=[0.5,0.25,0.15,0.10])
"""


########################### APARTADO 6: EVOLUCIÓN ###########################
# CASO 1: Gaussiano libre
"""

def psi_libre_analitica(x, t, x0, sigma0, k0):                      # Define la solución analítica del paquete gaussiano libre
    sigma_t  = sigma0 * np.sqrt(1 + (t**2) / (sigma0**4))           # Calcula el ancho σ(t) del paquete libre (dispersión temporal)
    x_mean_t = x0 + k0 * t                                          # Posición media del paquete: se desplaza con velocidad k0


    prefactor = 1.0 / ((2*np.pi)**0.25 * sigma_t)                   # Factor de normalización del gaussiano en tiempo t
    gauss = np.exp(-(x - x_mean_t)**2 / (4 * sigma_t**2))           # Parte gaussiana real centrada en x_mean_t con ancho sigma(t)
    phase = np.exp(1j * k0 * (x - x_mean_t - 0.5 * k0 * t))         # Fase compleja que incorpora el momento medio k0 y evolución temporal
    psi_t = prefactor * gauss * phase                               # Construye la función de onda completa phi(x,t) = prefactor·gauss·fase

    psi_t /= np.sqrt(np.trapz(np.abs(psi_t)**2, x))                 # Renormaliza phi(x,t) en la malla para asegurar que esta normalizado
    return psi_t                                                    # Devuelve la función de onda analítica en el tiempo t



# Parámetros espaciales y estado inicial
L_libre = 50.0               # Semilongitud del dominio espacial [-L_libre, L_libre], grande para evitar interacción con paredes
N_libre = 2000               # Número de puntos espaciales
x_libre = np.linspace(-L_libre, L_libre, N_libre)   # Malla espacial uniforme
dx_libre = x_libre[1] - x_libre[0]                  # Paso espacial real de la malla


x0_libre    = -15.0          # Posición inicial del centro del paquete
k0_libre    = 1.0            # Momento medio inicial (velocidad del paquete)
sigma_libre = 1.0            # Ancho inicial σ(0) del paquete gaussiano


psi0_libre = np.exp(-(x_libre - x0_libre)**2/(2*sigma_libre**2)) \  # Parte gaussiana real del estado inicial
             * np.exp(1j * k0_libre * x_libre)                      # Fase inicial con momento medio k0_libre
psi0_libre /= np.sqrt(np.trapz(np.abs(psi0_libre)**2, x_libre))     # Normaliza el estado inicial en la malla

V_libre = np.zeros_like(x_libre)

dt    = 0.01       # paso de tiempo
steps = 1000       # t_max = steps * dt


# Guardar estados y anchos
psi_cn_list  = []
psi_a_list   = []
sigma_cn_list = []
sigma_a_list  = []

psi_cn = psi0_libre.copy()
psi0_analitica = psi0_libre.copy()   # mismo estado inicial

for n in range(steps+1):
    t = n * dt

    # Guardar estado CN actual
    psi_cn_list.append(psi_cn.copy())

    # Estado analítico en t
    psi_a = psi_libre_analitica(x_libre, t, x0_libre, psi0_analitica, dx_libre)

    psi_a_list.append(psi_a.copy())

    # Ancho sigma_CN(t)
    dens_cn = np.abs(psi_cn)**2
    x_mean_cn  = np.trapz(x_libre * dens_cn, x_libre)
    x2_mean_cn = np.trapz(x_libre**2 * dens_cn, x_libre)
    sigma_cn_list.append(np.sqrt(x2_mean_cn - x_mean_cn**2))

    # Ancho sigma_a(t)
    dens_a = np.abs(psi_a)**2
    x_mean_a  = np.trapz(x_libre * dens_a, x_libre)
    x2_mean_a = np.trapz(x_libre**2 * dens_a, x_libre)
    sigma_a_list.append(np.sqrt(x2_mean_a - x_mean_a**2))

    # Evolución numérica CN (si no es el último paso)
    if n < steps:
        psi_cn = CN_step(psi_cn, V_libre, dx_libre, dt,
                         resolver="thomas", cc_tipo="dirichlet")
        psi_cn /= np.sqrt(np.trapz(np.abs(psi_cn)**2, x_libre))
        print(n, "de", steps)


# Animación
fig, ax = plt.subplots(figsize=(12,5))
line_cn, = ax.plot([], [], lw=2, label="CN")
line_a,  = ax.plot([], [], lw=2, ls=':', label="Analítica")

ax.set_xlim(x_libre[0], x_libre[-1])
ax.set_ylim(0, np.max(np.abs(psi0_libre)**2)*1.2)
ax.set_xlabel("x")
ax.set_ylabel("|ψ|²")
ax.set_title("Evolución del paquete gaussiano libre")
ax.legend()

def update(frame):
    line_cn.set_data(x_libre, np.abs(psi_cn_list[frame])**2)
    line_a.set_data(x_libre,  np.abs(psi_a_list[frame])**2)
    ax.set_title(f"t = {frame*dt:.3f}")
    return line_cn, line_a

ani = animation.FuncAnimation(fig, update,frames=range(0, steps+1, 10),blit=True)
plt.show()

ani.save("paquete_libre_dirichlet.gif", writer='pillow', fps=30)

# Evolución del ancho
t_array = np.arange(steps+1) * dt
plt.figure(figsize=(10,4))
plt.plot(t_array, sigma_cn_list, 'o-',  label="CN")
plt.plot(t_array, sigma_a_list,  '-',   label="Analítica")
plt.xlabel("t")
plt.ylabel("σ(t)")
plt.title("Evolución del ancho del paquete gaussiano libre")
plt.legend()
plt.grid()
plt.show()
"""






"""
# CASO 2: Gaussiano en pozo infinito

def psi_pozo_infinito_analitica(x, t, psi0, L, N_max=200):                  # Solución analítica en pozo infinito simétrico de [-L, L]
    psi_t = np.zeros_like(x, dtype=complex)                                 # Inicializa ψ(x,t) = 0 en toda la malla (tipo complejo)
    dentro = (x >= -L) & (x <= L)                                           # Máscara booleana: puntos dentro del pozo

    for n in range(1, N_max+1):
        phi_n = np.zeros_like(x, dtype=complex)
        phi_n[dentro] = np.sqrt(1/L) * np.sin(n*np.pi*(x[dentro] + L)/(2*L))
        E_n = (n*np.pi/(2*L))**2
        c_n = np.trapz(np.conj(phi_n) * psi0, x)
        psi_t += c_n * phi_n * np.exp(-1j * E_n * t)
        
    for n in range(1, N_max+1):                                             # Recorre los modos estacionarios n=1,...,N_max
        phi_n = np.zeros_like(x, dtype=complex)                             # Inicializa el modo propio \psi_n(x) a cero
        phi_n[dentro] = np.sqrt(1/L) * np.sin(                              
            n*np.pi*(x[dentro] + L)/(2*L)                                   # Define \psi_n(x) dentro del pozo: seno con nodos en ±L
        )
        E_n = (n*np.pi/(2*L))**2                                            # Energía del nivel n (m=1/2 en tus unidades)
        c_n = np.trapz(np.conj(phi_n) * psi0, x)                            # Coeficiente de expansión c_n = ⟨\psi_n|\phi(0)⟩
        psi_t += c_n * phi_n * np.exp(-1j * E_n * t)                        # Suma la contribución c_n \psi_n(x) e^{-iE_n t} al paquete total


    psi_t /= np.sqrt(np.trapz(np.abs(psi_t)**2, x))                         # Renormaliza ψ(x,t) para asegurar normalizacion
    return psi_t                                                            # Devuelve la solución analítica en el tiempo t



L_pozo = 10                                                                 # Semilongitud del pozo infinito (pozo en [-L_pozo, L_pozo])
N_pozo = 800                                                                # Número de puntos espaciales
x_pozo = np.linspace(-L_pozo, L_pozo, N_pozo)                               # Malla espacial uniforme en el pozo
dx_pozo = x_pozo[1] - x_pozo[0]                                            # Paso espacial real

x0_pozo = -5.0                                                              # Posición inicial del centro del paquete
k0_pozo = 2.0                                                               # Momento medio inicial
sigma_pozo = 1.0                                                            # Anchura inicial del paquete
psi0_pozo = np.exp(-(x_pozo - x0_pozo)**2/(2*sigma_pozo**2)) *             \
            np.exp(1j*k0_pozo*x_pozo)                                       # Paquete gaussiano inicial con fase e^{ikx}
psi0_pozo /= np.sqrt(np.trapz(np.abs(psi0_pozo)**2, x_pozo))               # Normaliza el estado inicial


V_pozo = np.zeros_like(x_pozo)                                             # Potencial V(x)=0 dentro del pozo (paredes se imponen numéricamente)
dt = 0.005                                                                 # Paso de tiempo
steps = 500                                                                # Número de pasos de tiempo (t_max = steps·dt)


psi_thomas = psi0_pozo.copy()
psi0_analitica = psi0_pozo.copy()

psi_thomas_list = []
psi_a_list = []

for n in range(steps+1):
    t = n * dt
    if n < steps:
        psi_thomas = CN_step(psi_thomas, V_pozo, dx_pozo, dt, resolver="thomas", cc_tipo="dirichlet")
        psi_thomas /= np.sqrt(np.trapz(np.abs(psi_thomas)**2, x_pozo))
    
    psi_a = psi_pozo_infinito_analitica(x_pozo, t, psi0_analitica, L_pozo)
    print(n, " de ", steps)
    
    psi_thomas_list.append(psi_thomas.copy())
    psi_a_list.append(psi_a.copy())


fig, ax = plt.subplots(figsize=(12,6))  # figura más grande
line_thomas, = ax.plot([], [], lw=2, ls='--', label="CN Thomas")
line_a, = ax.plot([], [], lw=2, ls=':', label="Analítica")

ax.set_xlim(x_pozo[0]*1.5, x_pozo[-1]*1.5)
ax.set_ylim(0, 1.5*np.max(np.abs(psi0_pozo)**2))  # eje Y más alto para ver paredes

# Representar paredes del pozo
V_visual = np.zeros_like(x_pozo)
V_visual[x_pozo <= -L_pozo] = 1.5*np.max(np.abs(psi0_pozo)**2)
V_visual[x_pozo >= L_pozo] = 1.5*np.max(np.abs(psi0_pozo)**2)
ax.fill_between(x_pozo, 0, V_visual, color='red', alpha=0.2, label="Paredes pozo")

ax.set_xlabel("x")
ax.set_ylabel("|ψ|²")
ax.legend()

def update(frame):
    line_thomas.set_data(x_pozo, np.abs(psi_thomas_list[frame])**2)
    line_a.set_data(x_pozo, np.abs(psi_a_list[frame])**2)
    ax.set_title(f"Paquete en pozo infinito t = {frame*dt:.3f}")
    return line_thomas, line_a

ani = animation.FuncAnimation(fig, update, frames=range(0, steps+1, 2), blit=True)
plt.show()

# Guardar animación como GIF
ani.save("pozo_infinito_con_paredes.gif", writer='pillow', fps=30)
"""








########################### APARTADO 7: DINÁMICA ###########################
"""
def psi_pozo_analitica(x, t, psi0, V, dx):                         # Solución “analítica” por expansión espectral para pozo/barrera genérico
    N = len(x)                                                      # Número de puntos de la malla espacial
    # Hamiltoniano H = -1/2 d^2/dx^2 + V(x)
    diag = 1/dx**2 + V                                              # Diagonal principal: término cinético + potencial
    off  = -0.5*np.ones(N-1)/dx**2                                  # Diagonal sub/superior: acoplos cinéticos entre vecinos
    H = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)          # Matriz tridiagonal completa del Hamiltoniano discretizado

    # Autovalores y autofunciones completas
    E, phi = np.linalg.eigh(H)                                      # Diagonaliza H: E[n] autovalores, phi[:,n] autofunciones

    psi_t = np.zeros_like(x, dtype=complex)                         # Inicializa ψ(x,t)=0
    for n in range(N):                                              # Recorre todos los modos propios n
        phi_n = phi[:, n]                                           # Autofunción φ_n(x)
        c_n = np.trapz(np.conj(phi_n) * psi0, x)                    # Coeficiente c_n = ⟨φ_n|ψ(0)⟩
        psi_t += c_n * phi_n * np.exp(-1j * E[n] * t)               # Suma c_n φ_n(x) e^{-i E_n t} a la superposición

    psi_t /= np.sqrt(np.trapz(np.abs(psi_t)**2, x))                 # Renormaliza ψ(x,t) para evitar errores numéricos
    return psi_t                                                    # Devuelve la función de onda en el tiempo t




def calcular_reflexion_transmision(x, psi_t,x_left=-2.0, x_right=+2.0):
    izquierda = x < x_left      # región de reflexión (muy a la izquierda del pozo)
    derecha   = x > x_right     # región de transmisión (muy a la derecha del pozo)

    R = np.trapz(np.abs(psi_t[izquierda])**2, x[izquierda])         # Probabilidad reflejada: \int_{x<x_left} |\ohi|^2 dx
    T = np.trapz(np.abs(psi_t[derecha])**2,   x[derecha])           # Probabilidad transmitida: \int_{x>x_right} |\phi|^2 dx
    return R, T                                                     # Devuelve R y T en ese instante


def graficar_instantes_pozo(x, psi0, V, dx, dt, steps,cc_tipo="dirichlet", N_instantes=3):
    psi_thomas = psi0.copy()                                       # Copia de çphi(0) para evolucionar numéricamente con CN+Thomas
    instantes = np.linspace(0, steps, N_instantes, dtype=int)      # Pasos de tiempo en los que se quiere comparar con la solución espectral

    R_list = []                                                    # Lista para guardar R(t) en cada paso
    T_list = []                                                    # Lista para guardar T(t) en cada paso

    for n in range(steps+1):                                       # Bucle temporal desde n=0 a n=steps
        t = n * dt                                                 # Tiempo físico actual
        print(n, " de ", steps)
        if n < steps:                                              # Evoluciona con CN mientras no sea el último paso
            psi_thomas = CN_step(psi_thomas, V, dx, dt,            # Un paso de CN con solver Thomas y CC elegidas
                                 resolver="thomas", cc_tipo=cc_tipo)
            psi_thomas /= np.sqrt(np.trapz(np.abs(psi_thomas)**2, x))  # Renormaliza para mantener la norma ≈1

        # Probabilidades reflejada / transmitida 
        R, T = calcular_reflexion_transmision(                     # Calcula R(t) y T(t) usando regiones asintóticas x_left/x_right
            x, psi_thomas, x_left=-2.0, x_right=+2.0)
        R_list.append(R)                                           # Guarda R(t)
        T_list.append(T)                                           # Guarda T(t)

        if n in instantes:                                         # Si este paso es uno de los instantes de comparación:
            psi_a = psi_pozo_analitica(x, t, psi0, V, dx)          # Calcula solución espectral phi_a(x,t)

            fig, ax = plt.subplots(figsize=(8,4))                  # Nueva figura para comparar perfiles en este t
            ax.plot(x, np.abs(psi_thomas)**2, lw=2, ls="--",       # Perfil numérico |ψ_CN|^2 (CN Thomas)
                    label="CN Thomas")
            ax.plot(x, np.abs(psi_a)**2,      lw=2, ls=":",        # Perfil “analítico” |\phi_a|^2 (expansión espectral)
                    label="Analítica pozo")
            # Pozo escalado (V < 0)
            ax.plot(
                x,
                V/np.max(np.abs(V)+1e-12)                          # Normaliza V(x) a [−1,0] aproximadamente
                  *np.max(np.abs(psi_thomas)**2)*0.5,              # Reescala para que quepa en la misma escala vertical que |\phi|^2
                'r--', lw=1, label="V(x) (escalado)"
            )

            ax.set_xlabel("x")
            ax.set_ylabel("|ψ|²")
            ax.set_title(f"Paquete sobre pozo rectangular, t = {t:.2f}")
            ax.legend()
            plt.show()

    # Evolución de R y T
    plt.figure(figsize=(10,4))
    plt.plot(np.arange(steps+1)*dt, R_list, 'b-', lw=2, label="Reflexión")
    plt.plot(np.arange(steps+1)*dt, T_list, 'g-', lw=2, label="Transmisión")
    plt.xlabel("t")
    plt.ylabel("Probabilidad")
    plt.title("Evolución de la reflexión y transmisión (pozo)")
    plt.grid()
    plt.legend()
    plt.show()


L = 100
N = 1000
x = np.linspace(-L, L, N)
dx = x[1] - x[0]

x0 = -5.0
k0 = 2.0
sigma = 1.0

psi0 = np.exp(-(x - x0)**2/(2*sigma**2)) * np.exp(1j*k0*x)
psi0 /= np.sqrt(np.trapz(np.abs(psi0)**2, x))

V0 = -10.0    # profundidad del pozo
a  = 2.0     # ancho
V = np.zeros_like(x)
V[np.abs(x) < a/2] = V0

dt = 0.0001
steps = 100000

graficar_instantes_pozo(x, psi0, V, dx, dt, steps,cc_tipo="dirichlet", N_instantes=3)
"""

"""

# CASO 2: Barrera rectangular o gaussiana

def psi_barrera_analitica(x, t, psi0, V, dx):                      # Análogo a psi_pozo_analitica pero aplicado a una barrera
    N = len(x)
    # Hamiltoniano H = -1/2 d^2/dx^2 + V(x)
    diag = 1/dx**2 + V                                             # Diagonal principal del Hamiltoniano
    off  = -0.5*np.ones(N-1)/dx**2                                 # Off-diagonales cinéticas
    H = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)         # Matriz tridiagonal H

    # Autovalores y autofunciones completas
    E, phi = np.linalg.eigh(H)                                     # Diagonalización: espectro y modos propios

    psi_t = np.zeros_like(x, dtype=complex)                        # Inicializa \phi(x,t)=0
    for n in range(N):                                             # Superpone todos los modos
        phi_n = phi[:, n]                                          # Autofunción \psi_n(x)
        c_n = np.trapz(np.conj(phi_n) * psi0, x)                   # Coeficiente de proyección c_n
        psi_t += c_n * phi_n * np.exp(-1j * E[n] * t)              # Contribución temporal c_n çpsi_n(x) e^{-iE_n t}

    psi_t /= np.sqrt(np.trapz(np.abs(psi_t)**2, x))                # Renormaliza
    return psi_t                                                   # Devuelve solución espectral en t



def calcular_reflexion_transmision(x, psi_t, x_barrera=0):

    izquierda = x < x_barrera
    derecha   = x > x_barrera

    R = np.trapz(np.abs(psi_t[izquierda])**2, x[izquierda])
    T = np.trapz(np.abs(psi_t[derecha])**2, x[derecha])
    return R, T


def graficar_instantes_barrera(x, psi0, V, dx, dt, steps, cc_tipo="dirichlet", N_instantes=3):
    psi_thomas = psi0.copy()
    instantes = np.linspace(0, steps, N_instantes, dtype=int)

    R_list = []
    T_list = []

    for n in range(steps+1):
        t = n * dt
        print(n, " de ", steps)
        if n < steps:
            psi_thomas = CN_step(
                psi_thomas, V, dx, dt,
                resolver="thomas", cc_tipo=cc_tipo
            )
            psi_thomas /= np.sqrt(np.trapz(np.abs(psi_thomas)**2, x))  # mantener norma

        # Probabilidades reflejada/transmitida
        R, T = calcular_reflexion_transmision(x, psi_thomas, x_barrera=0)
        R_list.append(R)
        T_list.append(T)
        
        
        if n in instantes:
            psi_a = psi_barrera_analitica(x, t, psi0, V, dx)

            fig, ax = plt.subplots(figsize=(8,4))
            ax.plot(x, np.abs(psi_thomas)**2, lw=2, ls="--", label="CN Thomas")
            ax.plot(x, np.abs(psi_a)**2, lw=2, ls=":", label="Analítica barrera")

            # Barrera escalada
            ax.plot(
                x,
                V/np.max(V+1e-12)*np.max(np.abs(psi_thomas)**2)*0.5,
                'r--', lw=1, label="V(x) (escalado)"
            )

            ax.set_xlabel("x")
            ax.set_ylabel("|ψ|²")
            ax.set_title(f"Paquete sobre barrera rectangular, t = {t:.2f}")
            ax.legend()
            plt.show()
            

    plt.figure(figsize=(10,4))
    plt.plot(np.arange(steps+1)*dt, R_list, 'b-', lw=2, label="Reflexión")
    plt.plot(np.arange(steps+1)*dt, T_list, 'g-', lw=2, label="Transmisión")
    plt.xlabel("t")
    plt.ylabel("Probabilidad")
    plt.title("Evolución de la reflexión y transmisión")
    plt.grid()
    plt.legend()
    plt.show()


L = 100
N = 1000
x = np.linspace(-L, L, N)
dx = x[1] - x[0]


x0 = -5.0
k0 = 2.0
sigma = 1.0

psi0 = np.exp(-(x - x0)**2/(2*sigma**2)) * np.exp(1j*k0*x)
psi0 /= np.sqrt(np.trapz(np.abs(psi0)**2, x))


V0 = 1.0     # altura de la barrera
a  = 2.0      # ancho
V = np.zeros_like(x)
V[np.abs(x) < a/2] = V0


dt = 0.0001
steps = 30000

graficar_instantes_barrera(x, psi0, V, dx, dt, steps, cc_tipo="dirichlet", N_instantes=3)

"""
































































