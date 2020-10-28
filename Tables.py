# -*- coding: utf-8 -*-
from scipy.interpolate import interp1d
from scipy.optimize import fminbound, bisect
import numpy as np

gamma = 1.4
table = {}


def writeTable(file, data):
    np.savetxt(file, data, delimiter=",")


def readTable(file):
    return np.genfromtxt(file, delimiter=",")


def deflectionAngleFromShock(delta, M):
    return np.arctan2(2 * np.cos(delta) * ((M * np.sin(delta))**2 - 1),
                      np.sin(delta) * (M**2 * (gamma + np.cos(2 * delta)) + 2))


def shockAngleFromDeflection(theta, M):
    # Find angle using bisection
    assert theta >= 0

    if theta <= maximumDeflectionAngle(M):
        return bisect(lambda d: (deflectionAngleFromShock(d, M) - theta),
                      *shockWaveRange(M))
    raise "Deflection angle too big. No attached shock."


def maximumDeflectionAngle(M):
    return table["MDA"](M)


def shockWaveRange(M):
    return table["SWAmin"](M), table["SWAmax"](M)


def main(plt=False):
    """ Generate tables """

    n = 200
    listM = np.logspace(0, 2, n)
    listMDA = np.empty(n)  # Maximum Deflection Angle
    listSWAmin = np.empty(n)  # Shock Wave Angle at minimum
    listSWAmax = np.empty(n)  # Shock Wave Angle at maximum

    iterM = enumerate(listM)
    next(iterM)  # Skip first iteration
    listMDA[0] = 0
    listSWAmin[0] = np.pi / 2
    listSWAmax[0] = np.pi / 2
    for i, M in iterM:
        xopt, fval, ierr, numfunc = fminbound(
            lambda t: -deflectionAngleFromShock(t, M), 0, np.pi / 2,
            full_output=True)
        assert ierr == 0
        listMDA[i] = -fval
        listSWAmax[i] = xopt

        x0, r = bisect(deflectionAngleFromShock, 0, xopt, args=M,
                       full_output=True)
        assert r.converged
        listSWAmin[i] = x0

    writeTable("DSMdata.csv", (listM, listMDA, listSWAmin, listSWAmax))
    table["MDA"] = interp1d(listM, listMDA, kind='cubic')
    table["SWAmin"] = interp1d(listM, listSWAmin, kind='cubic')
    table["SWAmax"] = interp1d(listM, listSWAmax, kind='cubic')

    # Tables Prandtl-Meyer


    if plt:
        from matplotlib.ticker import ScalarFormatter

        fig1, ax1 = plt.subplots()
        ax1.plot(listM, listMDA * 180 / np.pi,
                 label="Maximum deflection")
        ax1.plot(listM, listSWAmin * 180 / np.pi,
                 label="Shock wave minimum")
        ax1.plot(listM, listSWAmax * 180 / np.pi,
                 label="Shock wave maximum")
        ax1.legend()
        ax1.set_title("Shock wave angle at maximum and at no deflection"
                      " (γ = " + str(gamma) + ")")
        ax1.set_xlabel("Mach")
        ax1.set_ylabel("Angle (°)")
        ax1.set_xscale("log")
        ax1.grid(True, which="both")
        ax1.set_xticks([1, 2, 3, 5, 10, 20, 30, 50, 100])
        yticks = list(range(0, 90, 20))
        yticks.append(round(listSWAmax[-1] * 180 / np.pi))
        yticks.append(round(listMDA[-1] * 180 / np.pi))
        ax1.set_yticks(yticks)
        ax1.xaxis.set_major_formatter(ScalarFormatter())

        plt.show()


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    main(plt)
else:
    try:
        listM, listMDA, listSWAmin, listSWAmax = readTable("DSMdata.csv")
        table["MDA"] = interp1d(listM, listMDA, kind='cubic')
        table["SWAmin"] = interp1d(listM, listSWAmax, kind='cubic')
        table["SWAmax"] = interp1d(listM, listSWAmin, kind='cubic')
    except OSError as e:
        print(e)
        main()
