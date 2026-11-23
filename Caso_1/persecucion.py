import numpy as np
import matplotlib.pyplot as plt

# Parametros editables
H = 100.0          # altura del objeto [m]
v = 20.0           # rapidez horizontal del objeto [m/s]
a_motor = 20.0     # aceleracion del motor del cohete [m/s^2]
g = 9.81           # aceleracion gravitacional [m/s^2]
vc0x, vc0y = v, 0.0
dt, tmax = 0.005, 90.0
radio, tol = 1.0, 1e-3


def objeto(t):
    return np.array([v * t, H])


def f(t, Y):
    pos = Y[:2]
    vel = Y[2:]
    r_vec = objeto(t) - pos               # r_p = r_objeto - r_cohete
    r = np.linalg.norm(r_vec)
    r = r if r > 0 else 1e-30
    u = r_vec / r                          # vector unitario r_p(t)
    gravedad = np.array([0.0, -g])
    return np.concatenate([vel, a_motor * u + gravedad])


def rk4(t, Y, h):
    k1 = f(t, Y)
    k2 = f(t + h / 2, Y + h / 2 * k1)
    k3 = f(t + h / 2, Y + h / 2 * k2)
    k4 = f(t + h, Y + h * k3)
    return Y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6


def distancia(t, Y):
    return np.linalg.norm(objeto(t) - Y[:2])


def biseccion(t0, Y0, t1):
    a, b, Ya = t0, t1, Y0.copy()
    print("\nIteraciones de biseccion")
    print(" i |       a [s] |       b [s] |       m [s] | distancia [m] | cota error [s]")
    print("---+-------------+-------------+-------------+---------------+---------------")
    for i in range(1, 80):
        m = (a + b) / 2
        Ym = rk4(a, Ya, m - a)
        dm = distancia(m, Ym)
        if dm <= radio:
            b = m
        else:
            a, Ya = m, Ym
        error = (b - a) / 2
        print(f"{i:2d} | {a:11.6f} | {b:11.6f} | {m:11.6f} | {dm:13.6f} | {error:13.6e}")
        if error < tol:
            break
    Yb = rk4(a, Ya, b - a)
    return b, Yb


def simular():
    t, Y = 0.0, np.array([0.0, 0.0, vc0x, vc0y])
    ts, cohete, obj = [t], [Y[:2].copy()], [objeto(t)]
    while t < tmax:
        Y2, t2 = rk4(t, Y, dt), t + dt
        if distancia(t, Y) > radio >= distancia(t2, Y2):
            th, Yh = biseccion(t, Y, t2)
            ts.append(th); cohete.append(Yh[:2].copy()); obj.append(objeto(th))
            return ts, cohete, obj, (th, Yh)
        t, Y = t2, Y2
        ts.append(t); cohete.append(Y[:2].copy()); obj.append(objeto(t))
    return ts, cohete, obj, None


def main():
    ts, cohete, obj, hit = simular()
    if hit:
        th, Yh = hit
        print("\nIntercepcion: si")
        print(f"tiempo = {th:.6f} s")
        print(f"posicion = ({Yh[0]:.3f}, {Yh[1]:.3f}) m")
        print(f"rapidez cohete = {np.linalg.norm(Yh[2:]):.3f} m/s")
    else:
        print("\nIntercepcion: no")

    cohete = np.array(cohete)
    obj = np.array(obj)

    plt.plot(cohete[:, 0], cohete[:, 1], label="cohete")
    plt.plot(obj[:, 0], obj[:, 1], "--", label="objeto")
    plt.scatter([obj[-1, 0]], [obj[-1, 1]], color="black", s=25)
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.title("Persecucion con gravedad")
    plt.grid(True)
    plt.axis("equal")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
