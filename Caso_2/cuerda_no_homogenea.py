#!/usr/bin/env python3
"""Cuerda no homogenea: diferencias finitas + Jacobi.

Modelo implementado segun el documento:

    T X''(x) + omega^2 mu(x) X(x) = 0
    X(0) = X(L) = 0
    mu(x) = mu_1 + (mu_2 - mu_1) x / L

La discretizacion produce A X = lambda M X, lambda = omega^2.
Luego se transforma a C Y = lambda Y con C = D A D y D = M^(-1/2).
Los valores y vectores propios de C se obtienen con el metodo de Jacobi.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
from dataclasses import asdict, dataclass
from typing import Iterable


Matrix = list[list[float]]
Vector = list[float]


@dataclass
class StringResult:
    nodes: int
    length: float
    tension: float
    mu_1: float
    mu_2: float
    epsilon: float
    max_iterations: int
    dx: float
    x_full: Vector
    x_inner: Vector
    mu_inner: Vector
    A: Matrix
    M: Matrix
    D: Matrix
    C_initial: Matrix
    eigenvalues: Vector
    y_vectors: Matrix
    x_vectors: Matrix
    x_modes_full: Matrix
    angular_frequencies: Vector
    ordinary_frequencies: Vector
    jacobi_iterations: int
    jacobi_converged: bool
    jacobi_max_offdiag: float


def zeros(n: int, m: int) -> Matrix:
    return [[0.0 for _ in range(m)] for _ in range(n)]


def identity(n: int) -> Matrix:
    out = zeros(n, n)
    for i in range(n):
        out[i][i] = 1.0
    return out


def diagonal(values: Iterable[float]) -> Matrix:
    values = list(values)
    out = zeros(len(values), len(values))
    for i, value in enumerate(values):
        out[i][i] = value
    return out


def matmul(a: Matrix, b: Matrix) -> Matrix:
    rows = len(a)
    inner = len(b)
    cols = len(b[0]) if b else 0
    out = zeros(rows, cols)
    for i in range(rows):
        for k in range(inner):
            aik = a[i][k]
            if aik == 0.0:
                continue
            for j in range(cols):
                out[i][j] += aik * b[k][j]
    return out


def column(a: Matrix, j: int) -> Vector:
    return [row[j] for row in a]


def set_column(a: Matrix, j: int, values: Vector) -> None:
    for i, value in enumerate(values):
        a[i][j] = value


def linear_density(x: float, length: float, mu_1: float, mu_2: float) -> float:
    return mu_1 + (mu_2 - mu_1) * x / length


# ---------------------------------------------------------------------------
# Algoritmo 1: Construccion de las matrices A y M
# ---------------------------------------------------------------------------


def algoritmo_1_construccion_matrices(
    nodes: int,
    length: float,
    tension: float,
    mu_1: float,
    mu_2: float,
) -> tuple[float, Vector, Vector, Vector, Matrix, Matrix]:
    """Entrada: L, T, mu(x) definida por mu_1/mu_2 y N.

    Salida: matrices A y M de tamano N x N.

    Nota: se usa T/(Delta x)^2, consistente con la ecuacion diferencial
    y con el ejemplo del documento. En el pseudocodigo aparece T/Delta x,
    pero eso contradice el ejemplo numerico.
    """
    validate_physical_inputs(nodes, length, tension, mu_1, mu_2)

    # 1. Calcular Delta x.
    dx = length / (nodes + 1)

    # 2. Definir posiciones de los nodos interiores.
    x_full = [j * dx for j in range(nodes + 2)]
    x_inner = x_full[1:-1]

    # 3. Inicializar A y M.
    a = zeros(nodes, nodes)
    m = zeros(nodes, nodes)

    # 4. Para i = 1, ..., N.
    coef = tension / (dx * dx)
    for i in range(nodes):
        # 5. Evaluar la densidad mu_i = mu(x_i).
        mu_i = linear_density(x_inner[i], length, mu_1, mu_2)

        # 6. Asignar M_ii <- mu_i.
        m[i][i] = mu_i

        # 7. Asignar A_ii y los vecinos de la matriz tridiagonal.
        a[i][i] = 2.0 * coef
        if i > 0:
            a[i][i - 1] = -coef
        if i < nodes - 1:
            a[i][i + 1] = -coef

    mu_inner = [m[i][i] for i in range(nodes)]
    return dx, x_full, x_inner, mu_inner, a, m


# ---------------------------------------------------------------------------
# Algoritmo 2: Transformacion del problema generalizado
# ---------------------------------------------------------------------------


def algoritmo_2_transformacion_problema_generalizado(
    a: Matrix,
    m: Matrix,
) -> tuple[Matrix, Matrix]:
    """Entrada: matrices A y M.

    Salida: matriz simetrica C y matriz D = M^(-1/2).
    """
    n = len(m)

    # 1. Inicializar D.
    d = zeros(n, n)

    # 2. Para i = 1, ..., N, calcular D_ii <- 1 / sqrt(M_ii).
    for i in range(n):
        if m[i][i] <= 0:
            raise ValueError("La matriz M debe tener diagonal positiva.")
        d[i][i] = 1.0 / math.sqrt(m[i][i])

    # 3. Construir C <- D A D.
    c = matmul(matmul(d, a), d)
    symmetrize_in_place(c)
    return d, c


def symmetrize_in_place(a: Matrix) -> None:
    n = len(a)
    for i in range(n):
        for j in range(i + 1, n):
            value = 0.5 * (a[i][j] + a[j][i])
            a[i][j] = value
            a[j][i] = value


def largest_offdiag(a: Matrix) -> tuple[float, int, int]:
    n = len(a)
    best_value = 0.0
    best_p = 0
    best_q = 1 if n > 1 else 0
    for p in range(n - 1):
        for q in range(p + 1, n):
            value = abs(a[p][q])
            if value > best_value:
                best_value = value
                best_p = p
                best_q = q
    return best_value, best_p, best_q


# ---------------------------------------------------------------------------
# Algoritmo 3: Metodo de Jacobi
# ---------------------------------------------------------------------------


def algoritmo_3_metodo_jacobi(
    c_input: Matrix,
    tolerance: float = 1.0e-10,
    max_iterations: int = 100_000,
) -> tuple[Vector, Matrix, Matrix, int, bool, float]:
    """Metodo de Jacobi para C Y = lambda Y.

    Devuelve valores propios, vectores propios Y en columnas, la matriz casi
    diagonal final, iteraciones usadas, bandera de convergencia y mayor termino
    fuera de la diagonal restante.
    """
    if tolerance <= 0:
        raise ValueError("La tolerancia debe ser positiva.")
    if max_iterations < 1:
        raise ValueError("El numero maximo de iteraciones debe ser mayor que cero.")

    n = len(c_input)
    c = [row[:] for row in c_input]

    # 1. Inicializar V <- I_N.
    v = identity(n)

    if n == 1:
        return [c[0][0]], v, c, 0, True, 0.0

    last_offdiag = math.inf

    # 2. Para k = 1, ..., N_max.
    for iteration in range(1, max_iterations + 1):
        # 3-8. Buscar p < q correspondiente al mayor elemento fuera de la diagonal.
        last_offdiag, p, q = largest_offdiag(c)

        # 9-10. Detener si |C_pq| < epsilon.
        if last_offdiag < tolerance:
            return [c[i][i] for i in range(n)], v, c, iteration - 1, True, last_offdiag

        # 11. Calcular el angulo de rotacion theta en radianes.
        theta = 0.5 * math.atan2(2.0 * c[p][q], c[p][p] - c[q][q])
        cos_t = math.cos(theta)
        sin_t = math.sin(theta)

        # 12-16. Construir Q y actualizar C <- Q^T C Q, V <- V Q.
        apply_jacobi_rotation(c, v, p, q, cos_t, sin_t)

    # 17-18. Extraer la diagonal como lambdas y las columnas de V como Y^(n).
    last_offdiag, _, _ = largest_offdiag(c)
    return [c[i][i] for i in range(n)], v, c, max_iterations, False, last_offdiag


def apply_jacobi_rotation(
    c: Matrix,
    v: Matrix,
    p: int,
    q: int,
    cos_t: float,
    sin_t: float,
) -> None:
    """Actualiza C <- Q^T C Q y V <- V Q sin construir Q completa."""
    n = len(c)

    c_pp = c[p][p]
    c_qq = c[q][q]
    c_pq = c[p][q]

    for k in range(n):
        if k == p or k == q:
            continue
        c_kp = c[k][p]
        c_kq = c[k][q]
        new_kp = cos_t * c_kp + sin_t * c_kq
        new_kq = -sin_t * c_kp + cos_t * c_kq
        c[k][p] = new_kp
        c[p][k] = new_kp
        c[k][q] = new_kq
        c[q][k] = new_kq

    c[p][p] = cos_t * cos_t * c_pp + 2.0 * cos_t * sin_t * c_pq + sin_t * sin_t * c_qq
    c[q][q] = sin_t * sin_t * c_pp - 2.0 * cos_t * sin_t * c_pq + cos_t * cos_t * c_qq
    c[p][q] = 0.0
    c[q][p] = 0.0

    for i in range(n):
        v_ip = v[i][p]
        v_iq = v[i][q]
        v[i][p] = cos_t * v_ip + sin_t * v_iq
        v[i][q] = -sin_t * v_ip + cos_t * v_iq


def sort_eigenpairs(eigenvalues: Vector, y_vectors: Matrix) -> tuple[Vector, Matrix]:
    order = sorted(range(len(eigenvalues)), key=lambda i: eigenvalues[i])
    sorted_values = [eigenvalues[i] for i in order]
    sorted_vectors = zeros(len(y_vectors), len(order))
    for new_col, old_col in enumerate(order):
        values = column(y_vectors, old_col)
        # Fija una convencion de signo reproducible: el mayor componente positivo.
        if values:
            pivot = max(range(len(values)), key=lambda i: abs(values[i]))
            if values[pivot] < 0:
                values = [-value for value in values]
        set_column(sorted_vectors, new_col, values)
    return sorted_values, sorted_vectors


# ---------------------------------------------------------------------------
# Algoritmo 4: Recuperacion de los vectores propios originales
# ---------------------------------------------------------------------------


def algoritmo_4_recuperacion_vectores_originales(
    d: Matrix,
    y_vectors: Matrix,
) -> tuple[Matrix, Matrix]:
    """Entrada: D = M^(-1/2) y vectores propios Y^(n).

    Salida: vectores propios originales X^(n) y modos en toda la malla.
    """
    n = len(y_vectors)
    x_vectors = zeros(n, n)
    x_modes_full = zeros(n + 2, n)

    # 1. Para n = 1, ..., N.
    for mode in range(n):
        # 2. Calcular X^(n) <- D Y^(n).
        y_mode = column(y_vectors, mode)
        x_mode = matvec_diagonal(d, y_mode)

        # 3. Normalizar X^(n), si es necesario.
        # En este caso se conservan las amplitudes resultantes del cambio
        # X = D Y para coincidir con los valores del documento.
        set_column(x_vectors, mode, x_mode)

        # 4-5. Construir el vector de toda la malla e imponer fronteras fijas.
        x_modes_full[0][mode] = 0.0
        for i, value in enumerate(x_mode, start=1):
            x_modes_full[i][mode] = value
        x_modes_full[n + 1][mode] = 0.0

    return x_vectors, x_modes_full


def matvec_diagonal(d: Matrix, x: Vector) -> Vector:
    return [d[i][i] * x[i] for i in range(len(x))]


# ---------------------------------------------------------------------------
# Algoritmo 5: Calculo de las frecuencias naturales
# ---------------------------------------------------------------------------


def algoritmo_5_calculo_frecuencias(eigenvalues: Vector) -> tuple[Vector, Vector]:
    """Entrada: valores propios lambda_1, ..., lambda_N.

    Salida: frecuencias naturales angulares omega_n y ordinarias f_n.
    """
    angular_frequencies = []
    ordinary_frequencies = []

    # 1. Para n = 1, ..., N.
    for eigenvalue in eigenvalues:
        # 2. Calcular omega_n <- sqrt(lambda_n).
        omega = math.sqrt(max(eigenvalue, 0.0))
        angular_frequencies.append(omega)

        # 3. Calcular, si se requiere, f_n <- omega_n / (2 pi).
        ordinary_frequencies.append(omega / (2.0 * math.pi))

    return angular_frequencies, ordinary_frequencies


# ---------------------------------------------------------------------------
# Algoritmo 6: Solucion numerica de la cuerda no homogenea
# ---------------------------------------------------------------------------


def algoritmo_6_solucion_numerica_cuerda(
    nodes: int,
    length: float,
    tension: float,
    mu_1: float,
    mu_2: float,
    tolerance: float,
    max_iterations: int,
) -> StringResult:
    # 1. Construir las matrices A, M.
    dx, x_full, x_inner, mu_inner, a, m = algoritmo_1_construccion_matrices(
        nodes, length, tension, mu_1, mu_2
    )

    # 2-3. Construir D = M^(-1/2) y C = D A D.
    d, c = algoritmo_2_transformacion_problema_generalizado(a, m)

    # 4-13. Aplicar Jacobi a C para obtener lambdas y vectores Y.
    eigenvalues, y_vectors, _c_final, iterations, converged, max_offdiag_value = (
        algoritmo_3_metodo_jacobi(c, tolerance, max_iterations)
    )

    # 14. Ordenar los valores propios y las columnas correspondientes de V.
    eigenvalues, y_vectors = sort_eigenpairs(eigenvalues, y_vectors)

    # 15-19. Recuperar X^(n) e incorporar condiciones de frontera.
    x_vectors, x_modes_full = algoritmo_4_recuperacion_vectores_originales(d, y_vectors)

    # 20. Calcular omega_n y f_n.
    angular_frequencies, ordinary_frequencies = algoritmo_5_calculo_frecuencias(eigenvalues)

    return StringResult(
        nodes=nodes,
        length=length,
        tension=tension,
        mu_1=mu_1,
        mu_2=mu_2,
        epsilon=tolerance,
        max_iterations=max_iterations,
        dx=dx,
        x_full=x_full,
        x_inner=x_inner,
        mu_inner=mu_inner,
        A=a,
        M=m,
        D=d,
        C_initial=c,
        eigenvalues=eigenvalues,
        y_vectors=y_vectors,
        x_vectors=x_vectors,
        x_modes_full=x_modes_full,
        angular_frequencies=angular_frequencies,
        ordinary_frequencies=ordinary_frequencies,
        jacobi_iterations=iterations,
        jacobi_converged=converged,
        jacobi_max_offdiag=max_offdiag_value,
    )


def validate_physical_inputs(
    nodes: int,
    length: float,
    tension: float,
    mu_1: float,
    mu_2: float,
) -> None:
    if nodes < 1:
        raise ValueError("El numero de nodos interiores debe ser al menos 1.")
    if length <= 0:
        raise ValueError("La longitud L debe ser positiva.")
    if tension <= 0:
        raise ValueError("La tension T debe ser positiva.")
    if mu_1 <= 0 or mu_2 <= 0:
        raise ValueError("Las densidades mu_1 y mu_2 deben ser positivas.")


def prompt_int(label: str, default: int | None = None) -> int:
    suffix = f" [{default}]" if default is not None else ""
    while True:
        raw = input(f"{label}{suffix}: ").strip()
        if not raw and default is not None:
            return default
        try:
            return int(raw)
        except ValueError:
            print("Ingrese un entero valido.")


def prompt_float(label: str, default: float | None = None) -> float:
    suffix = f" [{default}]" if default is not None else ""
    while True:
        raw = input(f"{label}{suffix}: ").strip()
        if not raw and default is not None:
            return default
        try:
            return float(raw)
        except ValueError:
            print("Ingrese un numero valido.")


def format_vector(values: Vector, precision: int) -> str:
    return "[" + ", ".join(f"{value:.{precision}g}" for value in values) + "]"


def format_matrix(matrix: Matrix, precision: int) -> str:
    if not matrix:
        return "[]"
    rows = []
    for row in matrix:
        rows.append("  " + format_vector(row, precision))
    return "[\n" + "\n".join(rows) + "\n]"


def print_results(result: StringResult, precision: int, include_matrices: bool) -> None:
    print("\n=== Parametros ===")
    print(f"Nodos interiores N: {result.nodes}")
    print(f"Longitud L: {result.length:.{precision}g}")
    print(f"Tension T: {result.tension:.{precision}g}")
    print(f"mu_1: {result.mu_1:.{precision}g}")
    print(f"mu_2: {result.mu_2:.{precision}g}")
    print(f"Epsilon: {result.epsilon:.{precision}g}")
    print(f"Iteraciones maximas: {result.max_iterations}")
    print(f"Delta x: {result.dx:.{precision}g}")
    print(f"Nodos interiores x_i: {format_vector(result.x_inner, precision)}")
    print(f"Densidades mu_i: {format_vector(result.mu_inner, precision)}")

    print("\n=== Metodo de Jacobi ===")
    print(f"Convergio: {'si' if result.jacobi_converged else 'no'}")
    print(f"Iteraciones: {result.jacobi_iterations}")
    print(f"Maximo fuera de la diagonal: {result.jacobi_max_offdiag:.{precision}g}")

    if include_matrices:
        print("\n=== Matrices del problema ===")
        print("A =")
        print(format_matrix(result.A, precision))
        print("M =")
        print(format_matrix(result.M, precision))
        print("D = M^(-1/2) =")
        print(format_matrix(result.D, precision))
        print("C = D A D =")
        print(format_matrix(result.C_initial, precision))

    print("\n=== Valores propios y frecuencias ===")
    print("n, lambda_n, omega_n [rad/s], f_n [Hz]")
    for i, (lam, omega, freq) in enumerate(
        zip(result.eigenvalues, result.angular_frequencies, result.ordinary_frequencies),
        start=1,
    ):
        print(f"{i}, {lam:.{precision}g}, {omega:.{precision}g}, {freq:.{precision}g}")

    print("\n=== Vectores propios Y del problema transformado ===")
    print("Columnas: Y^(n)")
    print(format_matrix(result.y_vectors, precision))

    print("\n=== Vectores propios X del problema original ===")
    print("Columnas: X^(n), solo nodos interiores")
    print(format_matrix(result.x_vectors, precision))

    print("\n=== Modos X sobre toda la malla ===")
    print("Incluye fronteras fijas: X_0 = X_(N+1) = 0")
    print(format_matrix(result.x_modes_full, precision))


def save_json(result: StringResult, path: str) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(asdict(result), handle, indent=2)


def save_csv(result: StringResult, prefix: str) -> None:
    summary_path = f"{prefix}_frecuencias.csv"
    with open(summary_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["modo", "lambda", "omega_rad_s", "frecuencia_hz"])
        for i, (lam, omega, freq) in enumerate(
            zip(result.eigenvalues, result.angular_frequencies, result.ordinary_frequencies),
            start=1,
        ):
            writer.writerow([i, lam, omega, freq])

    y_path = f"{prefix}_vectores_Y.csv"
    write_vector_matrix_csv(y_path, result.y_vectors, "nodo_interior")

    x_path = f"{prefix}_vectores_X.csv"
    write_vector_matrix_csv(x_path, result.x_vectors, "nodo_interior")

    modes_path = f"{prefix}_modos_malla.csv"
    with open(modes_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["x"] + [f"X_modo_{i + 1}" for i in range(result.nodes)])
        for row_index, x_value in enumerate(result.x_full):
            writer.writerow([x_value] + result.x_modes_full[row_index])


def write_vector_matrix_csv(path: str, matrix: Matrix, first_column_name: str) -> None:
    rows = len(matrix)
    cols = len(matrix[0]) if rows else 0
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([first_column_name] + [f"modo_{i + 1}" for i in range(cols)])
        for i in range(rows):
            writer.writerow([i + 1] + matrix[i])


def simulate_with_vispy(
    result: StringResult,
    mode_count: int,
    amplitude: float,
    speed: float,
    backend: str,
) -> None:
    try:
        import numpy as np
        from vispy import app, scene
    except ImportError as exc:
        raise RuntimeError(
            "VisPy no esta instalado. Instale dependencias con: "
            "python3 -m pip install -r requirements.txt"
        ) from exc

    selected_backend = select_vispy_backend(app, backend)

    mode_count = max(1, min(mode_count, result.nodes))
    x = np.array(result.x_full, dtype=float)
    modes = np.array(result.x_modes_full, dtype=float)[:, :mode_count]
    omegas = np.array(result.angular_frequencies[:mode_count], dtype=float)

    for j in range(mode_count):
        scale = np.max(np.abs(modes[:, j]))
        if scale > 0:
            modes[:, j] /= scale

    weights = np.array([amplitude / (j + 1) for j in range(mode_count)], dtype=float)
    x_centered = x - 0.5 * result.length
    y0 = np.zeros_like(x)
    pos = np.column_stack([x_centered, y0])

    canvas = scene.SceneCanvas(
        keys="interactive",
        show=True,
        title=f"Cuerda no homogenea - modos normales ({selected_backend})",
        bgcolor="#111111",
        size=(1100, 650),
    )
    view = canvas.central_widget.add_view()
    view.camera = scene.PanZoomCamera(aspect=1)
    view.camera.set_range(
        x=(-0.55 * result.length, 0.55 * result.length),
        y=(-1.4 * amplitude, 1.4 * amplitude),
    )

    base_pos = np.column_stack([x_centered, np.zeros_like(x_centered)])
    scene.visuals.Line(base_pos, color=(0.35, 0.35, 0.35, 1.0), width=1.0, parent=view.scene)
    line = scene.visuals.Line(pos, color=(0.08, 0.72, 1.0, 1.0), width=3.0, parent=view.scene)
    nodes = scene.visuals.Markers(parent=view.scene)
    nodes.set_data(pos, face_color=(1.0, 0.9, 0.25, 1.0), size=7)

    title = (
        f"N={result.nodes}, L={result.length:g}, mu1={result.mu_1:g}, "
        f"mu2={result.mu_2:g}, modos={mode_count}"
    )
    scene.visuals.Text(
        title,
        color="white",
        pos=(0, 0.92 * amplitude),
        anchor_x="center",
        parent=view.scene,
        font_size=11,
    )

    start = time.perf_counter()

    def update(_event) -> None:
        elapsed = (time.perf_counter() - start) * speed
        phases = np.cos(omegas * elapsed)
        y = modes @ (weights * phases)
        pos[:, 1] = y
        line.set_data(pos)
        nodes.set_data(pos, face_color=(1.0, 0.9, 0.25, 1.0), size=7)
        canvas.update()

    timer = app.Timer(interval=1.0 / 60.0, connect=update, start=True)
    try:
        app.run()
    finally:
        timer.stop()


def select_vispy_backend(app_module, backend: str) -> str:
    candidates = ["glfw", "pyqt6"] if backend == "auto" else [backend]
    errors = []

    for candidate in candidates:
        try:
            app_module.use_app(candidate)
            return candidate
        except Exception as exc:
            errors.append(f"{candidate}: {exc}")

    details = "; ".join(errors)
    raise RuntimeError(
        "No se pudo inicializar ningun backend de VisPy. "
        f"Intentos: {details}. "
        "En Manjaro/Arch, si fuerza --backend pyqt6 y falla con xcb, instale "
        "el paquete del sistema xcb-util-cursor."
    )


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Modelo de una cuerda no homogenea con el metodo descrito en el PDF.",
    )
    parser.add_argument("-N", "--nodes", type=int, help="numero de nodos interiores")
    parser.add_argument("-L", "--length", type=float, help="longitud de la cuerda")
    parser.add_argument("--mu1", type=float, help="densidad lineal en x=0")
    parser.add_argument("--mu2", type=float, help="densidad lineal en x=L")
    parser.add_argument("-T", "--tension", type=float, default=1.0, help="tension constante")
    parser.add_argument(
        "--epsilon",
        "--tol",
        dest="tol",
        type=float,
        default=None,
        help="epsilon/tolerancia de convergencia para Jacobi",
    )
    parser.add_argument("--max-iter", type=int, default=100_000, help="iteraciones maximas")
    parser.add_argument("--precision", type=int, default=8, help="cifras significativas impresas")
    parser.add_argument("--show-matrices", action="store_true", help="imprime A, M, D y C")
    parser.add_argument("--save-json", help="guarda todos los resultados en un archivo JSON")
    parser.add_argument("--save-csv-prefix", help="guarda frecuencias y vectores en CSV con este prefijo")
    parser.add_argument("--simulate", action="store_true", help="abre una simulacion con VisPy")
    parser.add_argument("--sim-modes", type=int, default=3, help="numero de modos en la simulacion")
    parser.add_argument("--sim-amplitude", type=float, default=0.18, help="amplitud grafica de la simulacion")
    parser.add_argument("--sim-speed", type=float, default=1.0, help="factor de velocidad de la animacion")
    parser.add_argument(
        "--backend",
        choices=["auto", "glfw", "pyqt6"],
        default="auto",
        help="backend grafico de VisPy para la simulacion",
    )
    parser.add_argument(
        "--example",
        action="store_true",
        help="usa el caso del documento: N=3, L=1, T=1, mu1=1, mu2=2",
    )
    return parser.parse_args(argv)


def complete_args(args: argparse.Namespace) -> argparse.Namespace:
    if args.example:
        args.nodes = 3
        args.length = 1.0
        args.tension = 1.0
        args.mu1 = 1.0
        args.mu2 = 2.0
        if args.tol is None:
            args.tol = 1.0e-10
        return args

    interactive = any(value is None for value in (args.nodes, args.length, args.mu1, args.mu2))
    if interactive:
        print("Ingrese los parametros del modelo.")
    if args.nodes is None:
        args.nodes = prompt_int("Nodos interiores N", 3)
    if args.length is None:
        args.length = prompt_float("Longitud L", 1.0)
    if args.mu1 is None:
        args.mu1 = prompt_float("Densidad mu_1 en x=0", 1.0)
    if args.mu2 is None:
        args.mu2 = prompt_float("Densidad mu_2 en x=L", 2.0)
    if interactive and args.tol is None:
        args.tol = prompt_float("Epsilon para Jacobi", 1.0e-10)
    elif args.tol is None:
        args.tol = 1.0e-10
    return args


def main(argv: list[str] | None = None) -> int:
    args = complete_args(parse_args(argv or sys.argv[1:]))
    try:
        result = algoritmo_6_solucion_numerica_cuerda(
            nodes=args.nodes,
            length=args.length,
            tension=args.tension,
            mu_1=args.mu1,
            mu_2=args.mu2,
            tolerance=args.tol,
            max_iterations=args.max_iter,
        )
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2

    print_results(result, args.precision, args.show_matrices)

    if args.save_json:
        save_json(result, args.save_json)
        print(f"\nJSON guardado en: {os.path.abspath(args.save_json)}")

    if args.save_csv_prefix:
        save_csv(result, args.save_csv_prefix)
        print(f"\nCSV guardados con prefijo: {os.path.abspath(args.save_csv_prefix)}")

    if args.simulate:
        try:
            simulate_with_vispy(
                result,
                args.sim_modes,
                args.sim_amplitude,
                args.sim_speed,
                args.backend,
            )
        except RuntimeError as exc:
            print(f"\nNo se pudo iniciar la simulacion: {exc}", file=sys.stderr)
            return 3

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
