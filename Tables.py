# -*- coding: utf-8 -*-
from scipy.interpolate import interp1d
from scipy.optimize import fminbound, bisect
import numpy as np

gamma = 1.4
table = {}

fileDSM = "DSMdata.csv"  # Deflection-Shock-Mach (θ-β-M)
fileIFP = "IFPdata.csv"  # Isentropic Flow Properties
fileMNuMU = "MNuMUdata.csv"  # Prandtl-Meyer


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

# Partie Olivier


def p_sur_pisentropique(M):
    return (1 + ((gamma - 1) * M**2 / 2))**(-gamma / (gamma - 1))


def T_sur_Tisentropique(M):
    return (1 + ((gamma - 1) * M**2 / 2))**(-1)


def rho_sur_rhoisentropique(M):
    return (1 + ((gamma - 1) * M**2 / 2))**(-1 / (gamma - 1))

# Partie Anthony


def prandtlMeyerNuMuFromMach(M):
    return (np.sqrt((gamma + 1) / (gamma - 1)) * np.arctan(
            np.sqrt((gamma - 1) / (gamma + 1) * (M**2 - 1))) - np.arctan(
            np.sqrt(M**2 - 1))), np.arcsin(1 / M)

# def prandtlMeyerMachFromMu (Mu):
#   return 1/np.sin(Mu)

# def prandltMeyerMachFromNu (Nu):
#    return


def main(plt=False):
    """ Generate tables """

    n = 200
    listM1 = np.logspace(0, 2, n)
    listMDA = np.empty(n)  # Maximum Deflection Angle
    listSWAmin = np.empty(n)  # Shock Wave Angle at minimum
    listSWAmax = np.empty(n)  # Shock Wave Angle at maximum
    listMu = np.empty(n)  # Mu Prandlt Meyer
    listNu = np.empty(n)  # Nu Prandlt Meyer

    iterM = enumerate(listM1)
    next(iterM)  # Skip first iteration
    listMDA[0] = 0
    listSWAmin[0] = np.pi / 2
    listSWAmax[0] = np.pi / 2
    listMu[0] = 0
    listNu[0] = np.pi / 2
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
        listMu[i], listNu[i] = prandtlMeyerNuMuFromMach(M)

    # Partie Olivier
    n1 = 200
    n2 = 100
    listM21 = np.logspace(0, 1, n1, endpoint=False)  # Mach
    listM22 = np.logspace(1, 2, n2)
    listM2 = np.concatenate((listM21, listM22))

    # pressure/isentropic pressure
    listp_sur_pi = p_sur_pisentropique(listM2)
    # temperature/isentropic temperature
    listT_sur_Ti = T_sur_Tisentropique(listM2)
    # density/isentropic density
    listrho_sur_rhoi = rho_sur_rhoisentropique(listM2)

    writeTable(fileDSM, (listM1, listMDA, listSWAmin, listSWAmax))
    table["MDA"] = interp1d(listM1, listMDA, kind='cubic')
    table["SWAmin"] = interp1d(listM1, listSWAmin, kind='cubic')
    table["SWAmax"] = interp1d(listM1, listSWAmax, kind='cubic')

    writeTable(fileIFP, (listM2, listp_sur_pi,
                         listT_sur_Ti, listrho_sur_rhoi))
    table["p_sur_pi"] = interp1d(listM2, listp_sur_pi, kind='cubic')
    table["T_sur_Ti"] = interp1d(listM2, listT_sur_Ti, kind='cubic')
    table["rho_sur_rhoi"] = interp1d(listM2, listrho_sur_rhoi, kind='cubic')

    writeTable(fileMNuMU, (listM1, listMu, listNu))
    table["Mu"] = interp1d(listM1, listMu, kind='cubic')
    table["Nu"] = interp1d(listM1, listNu, kind='cubic')

    # Tables Prandtl-Meyer

    if plt:
        from matplotlib.ticker import ScalarFormatter

        fig1, ax1 = plt.subplots()
        ax1.plot(listM1, listMDA * 180 / np.pi,
                 label="Maximum deflection")
        ax1.plot(listM1, listSWAmin * 180 / np.pi,
                 label="Shock wave minimum")
        ax1.plot(listM1, listSWAmax * 180 / np.pi,
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

        fig2, ax2 = plt.subplots()
        ax2.plot(listM2, listp_sur_pi, label="P/P\u2080")
        ax2.plot(listM2, listT_sur_Ti, label="T/P\u2080")
        ax2.plot(listM2, listrho_sur_rhoi, label="ρ/ρ\u2080")
        ax2.legend()
        ax2.set_title("Isentropic flow properties (γ = " + str(gamma) + ")")
        ax2.set_xlabel("Mach")
        ax2.set_xscale("log")
        ax2.grid(True, which="both")
        ax2.set_xticks([1, 2, 3, 5, 10, 20, 30, 50, 100])
        ax2.xaxis.set_major_formatter(ScalarFormatter())

        fig3, ax3 = plt.subplots()
        ax3.plot(listM1, listMu * 180 / np.pi, label="Mu")
        ax3.plot(listM1, listNu * 180 / np.pi, label="Nu")
        ax3.legend()
        ax3.set_title("Prandtl-Meyer function and Mach angle"
                      " (γ = " + str(gamma) + ")")
        ax3.set_xlabel("Mach")
        ax3.set_ylabel("Angle (°)")
        ax3.set_xscale("log")
        ax3.grid(True, which="both")
        ax3.set_xticks([1, 2, 3, 5, 10, 20, 30, 50, 100])
        ax3.xaxis.set_major_formatter(ScalarFormatter())

        plt.show()


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    main(plt)
else:
    try:
        listM1, listMDA, listSWAmin, listSWAmax = readTable(fileDSM)
        table["MDA"] = interp1d(listM1, listMDA, kind='cubic')
        table["SWAmin"] = interp1d(listM1, listSWAmax, kind='cubic')
        table["SWAmax"] = interp1d(listM1, listSWAmin, kind='cubic')

        listM2, listp_sur_pi, listT_sur_Ti, listrho_sur_rhoi = readTable(
            fileIFP)
        table["p_sur_pi"] = interp1d(listM2, listp_sur_pi, kind='cubic')
        table["T_sur_Ti"] = interp1d(listM2, listT_sur_Ti, kind='cubic')
        table["rho_sur_rhoi"] = interp1d(listM2, listrho_sur_rhoi,
                                         kind='cubic')

        listM1, listMu, listNu = readTable(fileMNuMU)
        table["Mu"] = interp1d(listM1, listMu, kind='cubic')
        table["Nu"] = interp1d(listM1, listNu, kind='cubic')
    except OSError as e:
        print(e)
        main()
