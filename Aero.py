# -*- coding: utf-8 -*-
import tkinter as tk
import tkinter.filedialog as tkfd
import numpy as np
import json as js
import sys

import matplotlib.pyplot as plt
from Tables import main as tablePlot

from Tables import (
    maximumDeflectionAngle, shockAngleFromDeflection, machWaveAngle, maximumNu,
    prandtlMeyerMachFromNu, prandtlMeyerNuFromMach, prandltMeyerMuFromMach,
    p_sur_pisentropique, T_sur_Tisentropique, rho_sur_rhoisentropique,
    p2surp1fromMachSh, rho2surrho1fromMachSh, t2surt1fromMachSh,
    pressureFromAlt, temperatureFromAlt, densityFromAlt, soundSpeedFromAlt,
    cpFromM0p2surp1, mach2frommach1, soundSpeedFromTemp,
    readTable, plotAtmo, interp
)


def pointToBorder(x, y, angle, width, height):  # angle in rad
    # Get line from point and angle to border of screen
    angle %= 2 * np.pi
    if angle < np.pi / 4 or angle > 7 * np.pi / 4:  # Right (x = width) 315-45°
        return width, y + np.tan(angle) * (width - x)
    elif angle < 3 * np.pi / 4:  # Bottom (y = height) 45-135°
        return x + (height - y) * np.cos(angle) / np.sin(angle), height
    elif angle < 5 * np.pi / 4:  # Left (x = 0) 135-225°
        return 0, y - np.tan(angle) * x
    else:  # Top (y = 0) 135-225°
        return x - y * np.cos(angle) / np.sin(angle), 0


def round_sig(x, sig):
    return np.round(x, sig - int(np.floor(np.log10(np.abs(x) + 1e-6))) - 1)


class Zone:
    def __init__(self, mach, angle, parent=None,
                 z=None, p=None, t=None, rho=None):
        self.mach = mach
        self.angle = angle
        if z is None:
            self.altitude = z
            self.pressure = p
            self.temperature = t
            self.density = rho
        else:
            self.updateAlt(z)
        self.parent = parent

        if parent is None:
            self.i = 0
        else:
            self.i = parent.i + 1
        self.up = None

        self.lines = []

    def getForce(self, area):
        fx = self.pressure * area * np.cos(self.angle)
        fy = -self.pressure * area * np.sin(self.angle)
        return fx, fy

    def setCp(self, Cp):
        self.Cp = Cp

    def updateAlt(self, z):
        self.altitude = z
        self.pressure = pressureFromAlt(z)
        self.temperature = temperatureFromAlt(z)
        self.density = densityFromAlt(z)

    def addIncidentAngle(self, x, y, angle, color):
        self.lines.append((x, y, angle, color))

    def displayData(self, canvas, sig):
        text = "P = " + str(round_sig(self.pressure, sig)) + \
            "\nT = " + str(round_sig(self.temperature, sig)) + \
            "\nρ = " + str(round_sig(self.density, sig)) + \
            "\nM = " + str(round_sig(self.mach, sig))
        le = 70 if self.lines[0][3] == "red" else 140
        x = self.lines[0][0]
        y = self.lines[0][1] + (le if self.lines[0][2] > self.angle else -le)
        c = self.lines[0][3]

        return [canvas.createTextSpace(x, y + 10, text, fill=c),
                canvas.createVectorSpace(x, y - 30, 40, self.angle, c)]


class Camera:
    """ Handles coordinate transform between workspace and screen """

    def __init__(self, x=0, y=0, z=1):
        self.x, self.y, self.z = x, y, z

    def move(self, x, y):
        self.x += x
        self.y += y
        print(self.x, self.y)

    def updateZoom(self, z):
        self.z = z
        print(z)

    def screenToSpace(self, x, y):
        return x, -y

    def spaceToScreen(self, x, y):
        return x, -y

    def angleToScreen(self, a):
        return -a


def getZoneC(front, theta, up, x, y):
    # Get air after shock using clasic method
    dAngle = up * (theta - front.angle)

    # Get maximum attached angle
    maxAngle = maximumDeflectionAngle(front.mach)

    # TODO: dAngle > maxAngle
    if dAngle >= 0 and dAngle < maxAngle:
        if dAngle == 0:
            sAngle = machWaveAngle(front.mach)
        else:
            sAngle = shockAngleFromDeflection(dAngle, front.mach)

        Mn = front.mach * np.sin(sAngle)
        MnC = mach2frommach1(Mn)
        MachC = MnC / np.sin(sAngle - dAngle)
        if MachC <= 1:  # TODO: Deal with this
            return None

        # No temperature change:
        # MachC = front.mach * np.cos(sAngle) / np.cos(sAngle - dAngle)

        c = soundSpeedFromTemp(front.temperature)
        Vt1 = c * front.mach * np.cos(sAngle)

        # TODO: Put in Tables.py
        Rstar = 0.0083144621  # kg(km/s)²/(kmol.K)
        M_0 = 28.9644  # kg/kmol
        gamma = 1.4
        tC = M_0 * Vt1**2 / (1e6 * gamma * Rstar * (MachC**2 - MnC**2))

        Vn = c * Mn
        VnC = c * MnC

        rhoC = front.density * Vn / VnC

        pC = front.pressure + front.density * Vn**2 - rhoC * VnC**2

        zCompr = Zone(MachC, theta, front, p=pC, t=tC, rho=rhoC)

        zCompr.addIncidentAngle(x, y, front.angle + up * sAngle, "red")
        zCompr.setCp(cpFromM0p2surp1(front.mach, pC / front.pressure))
        return zCompr
    elif dAngle < 0 and -dAngle < maxAngle:
        mu1 = prandltMeyerMuFromMach(front.mach)
        newNu = prandtlMeyerNuFromMach(front.mach) - dAngle

        if newNu > maximumNu():  # TODO: Deal with this
            return None
        MachE = prandtlMeyerMachFromNu(newNu)
        mu2 = prandltMeyerMuFromMach(MachE)

        pf_sur_pi = p_sur_pisentropique(front.mach)
        pE_sur_pi = p_sur_pisentropique(MachE)
        Tf_sur_Ti = T_sur_Tisentropique(front.mach)
        TE_sur_Ti = T_sur_Tisentropique(MachE)
        rhof_sur_rhoi = rho_sur_rhoisentropique(front.mach)
        rhoE_sur_rhoi = rho_sur_rhoisentropique(MachE)

        zExp = Zone(MachE, theta, front,
                    p=front.pressure * pE_sur_pi / pf_sur_pi,
                    t=front.temperature * TE_sur_Ti / Tf_sur_Ti,
                    rho=front.density * rhoE_sur_rhoi / rhof_sur_rhoi)
        zExp.addIncidentAngle(x, y, front.angle + up * mu1, "green")
        zExp.addIncidentAngle(x, y, theta + up * mu2, "lime")
        zExp.setCp(cpFromM0p2surp1(front.mach, pE_sur_pi / pf_sur_pi))
        return zExp


def getZoneTTP(front, theta, up, x, y):
    dAngle = up * (theta - front.angle)

    # Get maximum attached angle
    maxAngle = maximumDeflectionAngle(front.mach)

    # TODO: dAngle > maxAngle
    if dAngle >= 0 and dAngle < maxAngle:
        if dAngle == 0:
            sAngle = machWaveAngle(front.mach)
        else:
            sAngle = shockAngleFromDeflection(dAngle, front.mach)

        # TODO: Put in Tables.py
        gamma = 1.4
        sAngle = dAngle * ((gamma + 1) / 4
                           + np.sqrt((gamma + 1)**2 / 16
                                     + 1 / (front.mach**2 * dAngle**2)))

        MSA2 = front.mach**2 * sAngle**2
        MachC = np.sqrt((gamma + 1)**2 * front.mach**2 * MSA2
                        - 4 * (MSA2 - 1) * (gamma * MSA2 + 1)
                        / ((2 * gamma * MSA2 - gamma + 1)
                            * ((gamma - 1) * MSA2 + 2)))
        if MachC <= 1:  # TODO: Deal with this
            return None

        pC = front.pressure + front.pressure * 2 * gamma\
            * (front.mach**2 * np.sin(sAngle)**2 - 1) / (gamma + 1)

        rhoC = front.density * (gamma + 1) * front.mach**2 * np.sin(sAngle)**2\
            / ((gamma - 1) * front.mach**2 * np.sin(sAngle)**2 + 2)

        tC = front.temperature * (pC / front.pressure) / (rhoC / front.density)

        zCompr = Zone(MachC, theta, front, p=pC, t=tC, rho=rhoC)

        zCompr.addIncidentAngle(x, y, front.angle + up * sAngle, "red")
        zCompr.setCp(cpFromM0p2surp1(front.mach, pC / front.pressure))
        return zCompr
    elif dAngle < 0 and -dAngle < maxAngle:
        mu1 = prandltMeyerMuFromMach(front.mach)
        newNu = prandtlMeyerNuFromMach(front.mach) - dAngle

        if newNu > maximumNu():  # TODO: Deal with this
            return None
        MachE = prandtlMeyerMachFromNu(newNu)
        mu2 = prandltMeyerMuFromMach(MachE)

        pf_sur_pi = p_sur_pisentropique(front.mach)
        pE_sur_pi = p_sur_pisentropique(MachE)
        Tf_sur_Ti = T_sur_Tisentropique(front.mach)
        TE_sur_Ti = T_sur_Tisentropique(MachE)
        rhof_sur_rhoi = rho_sur_rhoisentropique(front.mach)
        rhoE_sur_rhoi = rho_sur_rhoisentropique(MachE)

        zExp = Zone(MachE, theta, front,
                    p=front.pressure * pE_sur_pi / pf_sur_pi,
                    t=front.temperature * TE_sur_Ti / Tf_sur_Ti,
                    rho=front.density * rhoE_sur_rhoi / rhof_sur_rhoi)
        zExp.addIncidentAngle(x, y, front.angle + up * mu1, "green")
        zExp.addIncidentAngle(x, y, theta + up * mu2, "lime")
        zExp.setCp(cpFromM0p2surp1(front.mach, pE_sur_pi / pf_sur_pi))
        return zExp


def getZoneDT(front, theta, up, x, y):
    dAngle = up * (theta - front.angle)

    # Get maximum attached angle
    maxAngle = maximumDeflectionAngle(front.mach)

    # TODO: dAngle > maxAngle
    if dAngle >= 0 and dAngle < maxAngle:
        if dAngle == 0:
            sAngle = machWaveAngle(front.mach)
        else:
            sAngle = shockAngleFromDeflection(dAngle, front.mach)

        Mn = front.mach * np.sin(sAngle)
        MnC = mach2frommach1(Mn)
        MachC = MnC / np.sin(sAngle - dAngle)
        if MachC <= 1:  # TODO: Deal with this
            return None

        pCsurpf = p2surp1fromMachSh(front.mach, sAngle)
        rhoCsurrhof = rho2surrho1fromMachSh(front.mach, sAngle)
        TCsurTf = t2surt1fromMachSh(front.mach, sAngle)

        zCompr = Zone(MachC, theta, front,
                      p=front.pressure * pCsurpf,
                      t=front.temperature * TCsurTf,
                      rho=front.density * rhoCsurrhof)

        zCompr.addIncidentAngle(x, y, front.angle + up * sAngle, "red")
        zCompr.setCp(cpFromM0p2surp1(front.mach, pCsurpf))
        return zCompr
    elif dAngle < 0 and -dAngle < maxAngle:
        mu1 = prandltMeyerMuFromMach(front.mach)
        newNu = prandtlMeyerNuFromMach(front.mach) - dAngle

        if newNu > maximumNu():  # TODO: Deal with this
            return None
        MachE = prandtlMeyerMachFromNu(newNu)
        mu2 = prandltMeyerMuFromMach(MachE)

        pf_sur_pi = p_sur_pisentropique(front.mach)
        pE_sur_pi = p_sur_pisentropique(MachE)
        Tf_sur_Ti = T_sur_Tisentropique(front.mach)
        TE_sur_Ti = T_sur_Tisentropique(MachE)
        rhof_sur_rhoi = rho_sur_rhoisentropique(front.mach)
        rhoE_sur_rhoi = rho_sur_rhoisentropique(MachE)

        zExp = Zone(MachE, theta, front,
                    p=front.pressure * pE_sur_pi / pf_sur_pi,
                    t=front.temperature * TE_sur_Ti / Tf_sur_Ti,
                    rho=front.density * rhoE_sur_rhoi / rhof_sur_rhoi)
        zExp.addIncidentAngle(x, y, front.angle + up * mu1, "green")
        zExp.addIncidentAngle(x, y, theta + up * mu2, "lime")
        zExp.setCp(cpFromM0p2surp1(front.mach, pE_sur_pi / pf_sur_pi))
        return zExp


def calculateZones(front, getZone, angleBack, x1, y1, x2, y2):
    theta = np.arctan2(y2 - y1, x2 - x1)  # rad
    # if theta > np.pi / 2 or theta < -np.pi / 2:
    #     x1, x2 = x2, x1
    #     theta = (theta + np.pi / 2) % np.pi - np.pi / 2
    area = np.sqrt((x2 - x1)**2 + (y2 - y1)**2) * 1  # m^2 (unit depth)

    zUp1 = getZone(front, theta, 1, x1, y1)
    zDown1 = getZone(front, theta, -1, x1, y1)
    if zUp1 is not None and zDown1 is not None:
        fxU, fyU = zUp1.getForce(area)
        fxD, fyD = zDown1.getForce(area)
        C_N = zDown1.Cp - zUp1.Cp
        fx = fxD - fxU
        fy = fyD - fyU
    else:
        C_N = fx = fy = None

    newAngle = front.angle + angleBack

    if zUp1 is None:
        zUp2 = None
    else:
        zUp2 = getZone(zUp1, newAngle, 1, x2, y2)

    if zDown1 is None:
        zDown2 = None
    else:
        zDown2 = getZone(zDown1, newAngle, -1, x2, y2)

    if zUp2 is not None:
        zUp2.addIncidentAngle(x2, y2, newAngle, "gray")
    elif zDown2 is not None:
        zDown2.addIncidentAngle(x2, y2, newAngle, "gray")

    return ([zUp1, zDown1, zUp2, zDown2], C_N,
            fx, fy, (x1 + x2) / 2, (y1 + y2) / 2)


def calculateGeomZones(front, getZone, angleBack, conv,
                       points, vertex, geometry, depth):
    newAngle = front.angle + angleBack

    def sequence(seq, prec):
        zones, C_N, fx, fy, posfx, posfy = [], 0, 0, 0, None, None

        for s in seq:
            if type(s) == list:
                z, c, x, y, pfx, pfy = sequence(s, prec)
                zones += z
                if None in [c, C_N, x, y, pfx, pfy]:
                    C_N, fx, fy, posfx, posfy = None, None, None, None, None
                else:
                    C_N += c
                    fx += x
                    fy += y
                    posfx = pfx if posfx is None else (posfx + pfx) / 2
                    posfy = pfy if posfy is None else (posfy + pfy) / 2
            else:
                p1, p2, up = vertex[s]
                p1, p2 = points[p1], points[p2]

                d = 1 if depth is None or up == 0 else depth[s]
                area = np.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2) * d
                pfx = (p2[0] + p1[0]) / 2
                pfy = (p2[1] + p1[1]) / 2

                theta = np.arctan2(p2[1] - p1[1], p2[0] - p1[0])  # rad
                if theta > np.pi / 2 or theta < -np.pi / 2:
                    p1, p2 = p2, p1
                    theta = (theta + np.pi / 2) % np.pi - np.pi / 2

                zone = getZone(prec, theta, up, *p1)
                if zone is None:
                    return zones, None, None, None, None, None

                zone.up = up
                zones.append(zone)
                C_N += -up * zone.Cp
                x, y = zone.getForce(area)
                fx += -up * x
                fy += -up * y
                posfx = pfx if posfx is None else (posfx + pfx) / 2
                posfy = pfy if posfy is None else (posfy + pfy) / 2
                prec = zone

        if type(s) != list:
            p1, p2, up = vertex[s]
            p1, p2 = points[p1], points[p2]
            theta = np.arctan2(p2[1] - p1[1], p2[0] - p1[0])  # rad
            if theta > np.pi / 2 or theta < -np.pi / 2:
                p1, p2 = p2, p1
                theta = (theta + np.pi / 2) % np.pi - np.pi / 2

            lastZone = getZone(prec, newAngle - up * conv, up, *p2)
            if lastZone is not None:
                lastZone.addIncidentAngle(*p2, newAngle - up * conv, "gray")
            zones.append(lastZone)

        return zones, C_N, fx, fy, posfx, posfy

    return sequence(geometry, front)


class Canvas(tk.Canvas):
    def __init__(self, parent, camera, *args, **kwargs):
        tk.Canvas.__init__(self, parent, *args, **kwargs,
                           width=900, height=500, bd=0, highlightthickness=0)
        self.parent = parent
        self.camera = camera
        self.text = {}

    def createLine(self, x, y, a, c):
        x1, y1 = self.camera.spaceToScreen(x, y)
        x2, y2 = pointToBorder(x1, y1, self.camera.angleToScreen(a),
                               self.winfo_width(), self.winfo_height())
        return self.create_line(x1, y1, x2, y2, fill=c)

    def createTextSpace(self, x, y, text, **kwargs):
        return self.create_text(*self.camera.spaceToScreen(x, y),
                                text=text, **kwargs)

    def createTarget(self, x, y, r):
        x, y = self.camera.spaceToScreen(x, y)
        draw = []
        draw.append(self.create_circle(x, y, r, fill="black"))
        draw.append(self.create_circle_arc(x, y, r, fill="yellow",
                                           start=0, end=90))
        draw.append(self.create_circle_arc(x, y, r, fill="yellow",
                                           start=180, end=270))
        return draw

    def createVector(self, x, y, n, a, c=None, h=0.5):
        a = self.camera.angleToScreen(a)
        n1, n2 = n * h, n * (1 - h)
        return self.create_line(
            x - n1 * np.cos(a), y - n1 * np.sin(a),
            x + n2 * np.cos(a), y + n2 * np.sin(a),
            arrow=tk.LAST, fill=c)

    def createVectorSpace(self, x, y, n, a, c=None, h=0.5):
        return self.createVector(*self.camera.spaceToScreen(x, y), n, a, c, h)

    def updateVector(self, vid, x, y, n, a):
        a = self.camera.angleToScreen(a)
        self.coords(vid, x - n * np.cos(a) / 2, y - n * np.sin(a) / 2,
                    x + n * np.cos(a) / 2, y + n * np.sin(a) / 2)

    def updateVectorSpace(self, vid, x, y, n, a):
        self.updateVector(vid, *self.camera.spaceToScreen(x, y), n, a)

    def create_circle(self, x, y, r, **kwargs):
        return self.create_oval(x - r, y - r, x + r, y + r, **kwargs)

    def create_circle_arc(self, x, y, r, **kwargs):
        if "start" in kwargs and "end" in kwargs:
            kwargs["extent"] = kwargs["end"] - kwargs["start"]
            del kwargs["end"]
        return self.create_arc(x - r, y - r, x + r, y + r, **kwargs)


class AtmoApplication(tk.Frame):
    def __init__(self, parent, mApp, altInput, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        parent.title("Sélection atmosphère")
        parent.geometry("500x200")
        self.parent = parent
        self.mApp = mApp
        self.altInput = altInput

        buttons = tk.Frame(self)
        sigmas = tk.LabelFrame(self, text="Sigma")
        frame = tk.Frame(self)

        self.buttonOk = tk.Button(
            buttons, text="Ok", command=self.onOk)
        self.buttonCancel = tk.Button(
            buttons, text="Cancel", command=self.onCancel)

        self.buttonOk.grid(row=0, column=0)
        self.buttonCancel.grid(row=0, column=1)

        self.sig_default = mApp.text_sig
        self.sig = tk.IntVar(value=self.sig_default)
        tk.Label(sigmas, text="Displayed data").pack(side=tk.LEFT, padx=5)
        spinbox = tk.Spinbox(sigmas, from_=1, to=25, width=8,
                             textvariable=self.sig, command=self.updateSig)
        spinbox.pack(side=tk.LEFT, pady=5)
        spinbox.bind('<Return>', lambda _: self.updateSig())

        self.inputs = [
            UserInput(frame, "Altitude    ", "km   ", mApp.air.altitude,
                      0, 100, self.updateZ),
            UserInput(frame, "Density     ", "kg/m\u00B3", mApp.air.density,
                      0.002, 20, self.updateD, True),
            UserInput(frame, "Temperature ", "K    ", mApp.air.temperature,
                      0, 500, self.updateT),
            UserInput(frame, "Pressure    ", "Pa   ", mApp.air.pressure,
                      0.1, 1e6, self.updateP, True),
        ]
        for u in self.inputs:
            u.pack(fill=tk.X)

        buttons.pack(side=tk.BOTTOM, fill=tk.X, pady=5)
        buttons.grid_columnconfigure(0, weight=1)
        buttons.grid_columnconfigure(1, weight=1)

        sigmas.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=5)

        frame.pack(fill=tk.X, padx=10, pady=5)

    def updateZ(self, z):
        self.altInput.updateValue(z, False)
        self.inputs[1].updateValue(densityFromAlt(z), False)
        self.inputs[2].updateValue(temperatureFromAlt(z), False)
        self.inputs[3].updateValue(pressureFromAlt(z), False)
        self.mApp.userUpdateZ(z)

    def updateD(self, d):
        self.mApp.userUpdateD(d)

    def updateT(self, t):
        self.mApp.userUpdateT(t)

    def updateP(self, p):
        self.mApp.userUpdateP(p)

    def updateSig(self):
        self.mApp.userUpdateSig(self.sig.get())

    def onOk(self):
        for u in self.inputs:
            u.updateValue(u.var.get())
        self.mApp.userUpdateSig(self.sig.get())
        self.parent.destroy()

    def onCancel(self):
        z = self.inputs[0].default
        d = self.inputs[1].default
        t = self.inputs[2].default
        p = self.inputs[3].default
        self.altInput.updateValue(z, False)
        self.mApp.userUpdateZ(z)
        self.mApp.userUpdateD(d)
        self.mApp.userUpdateT(t)
        self.mApp.userUpdateP(p)
        self.mApp.userUpdateSig(self.sig_default)
        self.parent.destroy()


class TrApplication(tk.Frame):
    def __init__(self, parent, mApp, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        parent.title("Modification transformation")
        parent.geometry("500x220")
        self.parent = parent
        self.mApp = mApp

        buttons = tk.Frame(self)
        sigmas = tk.Frame(self)
        frame = tk.Frame(self)

        self.buttonOk = tk.Button(
            buttons, text="Ok", command=self.onOk)
        self.buttonCancel = tk.Button(
            buttons, text="Cancel", command=self.onCancel)

        self.buttonOk.grid(row=0, column=0)
        self.buttonCancel.grid(row=0, column=1)

        self.pRange = tk.IntVar(value=5 * round_sig(
            max(abs(mApp.dataT["x0"]), abs(mApp.dataT["y0"]), 100), 1))
        spinbox = tk.Spinbox(sigmas, from_=1, to=1e6, width=8,
                             textvariable=self.pRange,
                             command=self.updatePRange)
        spinbox.pack(side=tk.RIGHT)
        spinbox.bind('<Return>', lambda _: self.updatePRange())
        tk.Label(sigmas, text="Position range").pack(side=tk.RIGHT)

        self.lock = tk.IntVar(value=(
            1 if mApp.dataT["xs"] == mApp.dataT["ys"] else 0))
        tk.Checkbutton(sigmas, text="Lock Scale", onvalue=1, offvalue=0,
                       variable=self.lock, command=self.updateLock).pack(
            side=tk.LEFT, padx=5, pady=5)

        self.inverted = False
        tk.Button(sigmas, text="Invert Vertex",
                  command=self.invertVertex).pack(padx=5, pady=5)

        pRange = self.pRange.get()
        self.inputs = [  # Update PRange, Xs and Ys if order modified
            UserInput(frame, "Pos X    ", "px", mApp.dataT["x0"],
                      -pRange / 2, pRange / 2, self.updateX0),
            UserInput(frame, "Pos Y    ", "px", mApp.dataT["y0"],
                      -pRange / 2, pRange / 2, self.updateY0),
            UserInput(frame, "Scale X  ", "  ", mApp.dataT["xs"],
                      0.1, 10, self.updateXs, True),
            UserInput(frame, "Scale Y  ", "  ", mApp.dataT["ys"],
                      0.1, 10, self.updateYs, True),
            UserInput(frame, "Rotation ", "° ", mApp.dataT["rot"],
                      -180, 180, self.updateRot),
        ]
        for u in self.inputs:
            u.pack(fill=tk.X)

        buttons.pack(side=tk.BOTTOM, fill=tk.X, pady=5)
        buttons.grid_columnconfigure(0, weight=1)
        buttons.grid_columnconfigure(1, weight=1)

        sigmas.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=5)

        frame.pack(fill=tk.X, padx=10, pady=5)

    def updateX0(self, x0):
        self.mApp.dataT["x0"] = x0
        self.mApp.loadGeom()

    def updateY0(self, y0):
        self.mApp.dataT["y0"] = y0
        self.mApp.loadGeom()

    def updateXs(self, xs):
        if self.lock.get() == 1:
            self.inputs[3].updateValue(xs, False)
            self.mApp.dataT["ys"] = xs
        self.mApp.dataT["xs"] = xs
        self.mApp.loadGeom()

    def updateYs(self, ys):
        if self.lock.get() == 1:
            self.inputs[2].updateValue(ys, False)
            self.mApp.dataT["xs"] = ys
        self.mApp.dataT["ys"] = ys
        self.mApp.loadGeom()

    def updateRot(self, rot):
        self.mApp.dataT["rot"] = rot
        self.mApp.loadGeom()

    def updateLock(self):
        if self.lock.get() == 1:
            self.inputs[3].updateValue(self.mApp.dataT["xs"], False)
            self.mApp.dataT["ys"] = self.mApp.dataT["xs"]

    def updatePRange(self):
        pRange = self.pRange.get()
        self.inputs[0].updateScale(-pRange / 2, pRange / 2)
        self.inputs[1].updateScale(-pRange / 2, pRange / 2)

    def invertVertex(self):
        self.mApp.invertVertex()
        self.inverted = not self.inverted

    def onOk(self):
        for u in self.inputs:
            u.updateValue(u.var.get())
        self.parent.destroy()

    def onCancel(self):
        self.lock.set(0)
        if self.inverted:
            self.mApp.invertVertex()
        for u in self.inputs:
            u.updateValue(u.default)
        self.parent.destroy()


class CGApplication(tk.Frame):
    def __init__(self, parent, mApp, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        parent.title("Plot Center of Gravity")
        parent.geometry("500x200")
        self.parent = parent
        self.mApp = mApp

        self.defaultI = mApp.air.angle * 180 / np.pi
        self.defaultM = mApp.air.mach
        self.defaultA = mApp.air.altitude
        if self.defaultA is None:
            self.defaultPTD = (mApp.air.pressure,
                               mApp.air.temperature,
                               mApp.air.density)

        self.t, self.cgx, self.cgy = [0], [0], [0]
        self.currentT, self.x, self.y = 0, 0, 0
        self.resetIMA()

        load = tk.Button(self, text="Load CG File", command=self.loadCGFile)

        self.scale = UserInput(self, "Time ", "s", self.currentT,
                               0, 0, self.updateCG, fun=self.getTime)

        coords = tk.Frame(self)

        self.varX = tk.DoubleVar(value=0)
        self.varY = tk.DoubleVar(value=0)
        sx = tk.Spinbox(coords, textvariable=self.varX,
                        command=self.updateX)
        sy = tk.Spinbox(coords, textvariable=self.varY,
                        command=self.updateY)
        sx.bind('<Return>', self.updateX)
        sy.bind('<Return>', self.updateY)

        tk.Label(coords, text="X : ").grid(row=0, column=0)
        sx.grid(row=0, column=1)
        tk.Label(coords, text="Y : ").grid(row=0, column=2)
        sy.grid(row=0, column=3)

        plot = tk.Button(self, text="Plot CG over time", command=self.plotCG)

        files = tk.Frame(self)

        tk.Button(files, text="Load Angle of Attack File",
                  command=self.loadIncFile).grid(row=0, column=0, padx=5)
        tk.Button(files, text="Load Mach File",
                  command=self.loadMachFile).grid(row=0, column=1, padx=5)
        tk.Button(files, text="Load Altitude File",
                  command=self.loadAltFile).grid(row=0, column=2, padx=5)
        tk.Button(files, text="Reset Files",
                  command=self.resetIMA).grid(row=0, column=3, padx=5)

        load.pack(pady=10)
        self.scale.pack(fill=tk.X, pady=5)
        coords.pack(pady=5)
        plot.pack(pady=5)
        files.pack(pady=10)

        self.draw = []
        self.parent.protocol("WM_DELETE_WINDOW", self.onClose)

        self.drawCG()

    def updateX(self, e=None):
        self.x = self.varX.get()
        self.drawCG()

    def updateY(self, e=None):
        self.y = self.varY.get()
        self.drawCG()

    def drawCG(self):
        [[x, y]] = self.mApp.transformPoints([[self.x, self.y]])
        self.mApp.canvas.delete(*self.draw)
        self.draw = self.mApp.canvas.createTarget(x, y, 10)

    def plotCG(self):
        fig1, ax1 = plt.subplots(figsize=(10, 7))

        ax1.plot(self.t, self.cgx, label="X")
        ax1.plot(self.t, self.cgy, label="Y")
        ax1.legend()
        ax1.set_title("Position du CG")
        ax1.set_xlabel("Temps (s)")
        ax1.set_ylabel("Position (same units as json)")
        ax1.grid(True, which="both")

        plt.show()

    def getTime(self, value, inv=False):
        if inv:
            value = min(self.t, key=lambda x: abs(x - float(value)))
            return np.where(self.t == value)[0][0]

        return self.t[int(value)]

    def updateCG(self, value):
        self.currentT = value
        t = self.getTime(value, inv=True)

        if self.inc is not None:
            self.mApp.userUpdateA(self.inc(value))
        if self.mach is not None:
            self.mApp.userUpdateM(self.mach(value))
        if self.alt is not None:
            self.mApp.userUpdateZ(self.alt(value))

        self.x, self.y = self.cgx[t], self.cgy[t]
        self.varX.set(self.x)
        self.varY.set(self.y)
        self.drawCG()

    def onClose(self):
        self.mApp.canvas.delete(*self.draw)
        self.resetIMA()
        self.parent.destroy()

    def loadCGFile(self):
        fileName = tkfd.askopenfile(mode="r", filetypes=[
            ("CSV Files", "*.csv"), ("Any File", "*.*")],
            title="CSV file containing Time (s) "
                  "and CG position (same units as json).")

        if fileName is not None:
            self.t, self.cgx, self.cgy = readTable(fileName)

            self.scale.updateScale(0, len(self.t) - 1)
            self.updateCG(0)

    def interp(self, tlist, xlist):
        mint = min(tlist)
        maxt = max(tlist)
        interpfun = interp(tlist, xlist)

        def interpBorder(t):
            if t < mint:
                return xlist[mint]
            elif t > maxt:
                return xlist[maxt]
            return interpfun(t)

        return interpBorder

    def loadIncFile(self):
        fileName = tkfd.askopenfile(mode="r", filetypes=[
            ("CSV Files", "*.csv"), ("Any File", "*.*")],
            title="CSV file containing Time (s) and Angle of Attack (°).")

        if fileName is not None:
            self.inc = self.interp(*readTable(fileName))
            self.mApp.userUpdateA(self.inc(self.currentT))

    def loadMachFile(self):
        fileName = tkfd.askopenfile(mode="r", filetypes=[
            ("CSV Files", "*.csv"), ("Any File", "*.*")],
            title="CSV file containing Time (s) and Mach.")

        if fileName is not None:
            self.mach = self.interp(*readTable(fileName))
            self.mApp.userUpdateM(self.mach(self.currentT))

    def loadAltFile(self):
        fileName = tkfd.askopenfile(mode="r", filetypes=[
            ("CSV Files", "*.csv"), ("Any File", "*.*")],
            title="CSV file containing Time (s) and Altitude (km).")

        if fileName is not None:
            self.alt = self.interp(*readTable(fileName))
            self.mApp.userUpdateZ(self.alt(self.currentT))

    def resetIMA(self):
        self.inc, self.mach, self.alt = None, None, None
        self.mApp.userUpdateA(self.defaultI)
        self.mApp.userUpdateM(self.defaultM)
        if self.defaultA is None:
            p, t, d = self.defaultPTD
            self.mApp.userUpdateP(p)
            self.mApp.userUpdateT(t)
            self.mApp.userUpdateD(d)
        else:
            self.mApp.userUpdateZ(self.defaultA)


class Menubar(tk.Menu):
    def __init__(self, parent, *args, **kwargs):
        tk.Menu.__init__(self, parent, *args, **kwargs)
        self.parent = parent

        fileMenu = tk.Menu(self, tearoff=0)
        fileMenu.add_command(label="New", command=parent.newGeom)
        fileMenu.add_command(label="Open", command=parent.openGeom)
        fileMenu.add_command(label="Save", command=parent.saveGeom)
        fileMenu.add_separator()
        exampleNumber = 2
        for k in range(1, exampleNumber + 1):
            fileMenu.add_command(label="Example %i" % k,
                                 command=lambda k=k: parent.exampleGeom(k))
        fileMenu.add_separator()
        fileMenu.add_command(label="Quit", command=parent.onExit)

        toolsMenu = tk.Menu(self, tearoff=0)
        toolsMenu.add_command(label="Set atmosphere",
                              command=parent.openAtmo)
        self.trLabel = "Set transformation"
        toolsMenu.add_command(label=self.trLabel,
                              command=parent.openTr)
        toolsMenu.add_separator()
        toolsMenu.add_command(label="Reset inputs", command=parent.resetInputs)

        var = tk.IntVar(value=1)
        methodMenu = tk.Menu(self, tearoff=0)
        methodMenu.add_radiobutton(label="Classical",
                                   variable=var, value=1, command=self.updateM)
        methodMenu.add_radiobutton(label="Small-Disturbance Theory",
                                   variable=var, value=2, command=self.updateM)
        methodMenu.add_radiobutton(label="Tangent Wedge",
                                   variable=var, value=3, command=self.updateM)

        self.autoUpdate = tk.BooleanVar()
        self.autoUpdate.set(True)
        updateMenu = tk.Menu(self, tearoff=1)
        updateMenu.add_checkbutton(label="Auto Update", onvalue=1, offvalue=0,
                                   variable=self.autoUpdate,
                                   command=parent.updateZones)
        updateMenu.add_separator()
        updateMenu.add_command(label="Manual Update",
                               command=lambda: parent.updateZones(True))

        plotMenu = tk.Menu(self, tearoff=0)
        plotMenu.add_command(label="Plot Tables",
                             command=lambda: tablePlot(plt))
        plotMenu.add_command(label="Plot Atmosphere",
                             command=lambda: plotAtmo(plt))
        self.CGLabel = "Plot CG"
        plotMenu.add_command(label=self.CGLabel,
                             command=parent.openCG)

        plotMenu.add_command(label="Plot Everything",
                             command=self.plotAll)

        self.method = var
        self.toolsMenu = toolsMenu
        self.plotMenu = plotMenu
        self.disableTr()

        self.add_cascade(label="File", menu=fileMenu)
        self.add_cascade(label="Tools", menu=toolsMenu)
        self.add_cascade(label="Method", menu=methodMenu)
        self.add_cascade(label="Update", menu=updateMenu)
        self.add_cascade(label="Plot", menu=plotMenu)

    def enableTr(self):
        self.toolsMenu.entryconfig(self.trLabel, state=tk.NORMAL)
        self.plotMenu.entryconfig(self.CGLabel, state=tk.NORMAL)

    def disableTr(self):
        self.toolsMenu.entryconfig(self.trLabel, state=tk.DISABLED)
        self.plotMenu.entryconfig(self.CGLabel, state=tk.DISABLED)

    def plotAll(self):
        fileInc = tkfd.askopenfile(mode="r", filetypes=[
            ("CSV Files", "*.csv"), ("Any File", "*.*")],
            title="CSV file containing Time (s) and Angle of Attack (°).")

        fileSpeed = tkfd.askopenfile(mode="r", filetypes=[
            ("CSV Files", "*.csv"), ("Any File", "*.*")],
            title="CSV file containing Time (s) and Speed (m/s).")

        fileAlt = tkfd.askopenfile(mode="r", filetypes=[
            ("CSV Files", "*.csv"), ("Any File", "*.*")],
            title="CSV file containing Time (s) and Altitude (km).")

        if (fileInc is not None and fileSpeed is not None
                and fileAlt is not None):
            figsize = (10, 7)

            ti, i = readTable(fileInc)
            tv, v = readTable(fileSpeed)
            ta, a = readTable(fileAlt)
            inc = interp(ti, i)
            vel = interp(tv, v)
            alt = interp(ta, a)

            tmini = max(min(ti), min(tv), min(ta))
            tmaxi = min(max(ti), max(tv), max(ta))

            t = np.linspace(tmini, tmaxi, 500)

            fig1, ax1 = plt.subplots(figsize=figsize)

            ax1.plot(t, inc(t))
            ax1.set_title("Incidence en fonction du temps")
            ax1.set_xlabel("Temps (s)")
            ax1.set_ylabel("Incidence (°)")
            ax1.grid(True, which="both")

            fig2, ax2 = plt.subplots(figsize=figsize)

            ax2.plot(t, vel(t))
            ax2.set_title("Vitesse en fonction du temps")
            ax2.set_xlabel("Temps (s)")
            ax2.set_ylabel("Vitesse (m/s)")
            ax2.grid(True, which="both")

            fig3, ax3 = plt.subplots(figsize=figsize)

            ax3.plot(t, alt(t))
            ax3.set_title("Altitude en fonction du temps")
            ax3.set_xlabel("Temps (s)")
            ax3.set_ylabel("Altitude (km)")
            ax3.grid(True, which="both")

            # Mach

            fig4, ax4 = plt.subplots(figsize=figsize)

            ax4.plot(t, vel(t) / soundSpeedFromAlt(alt(t)))
            ax4.set_title("Mach en fonction du temps")
            ax4.set_xlabel("Temps (s)")
            ax4.set_ylabel("Mach")
            ax4.grid(True, which="both")

            # Pressure

            fig5, ax5 = plt.subplots(figsize=figsize)

            ax5.plot(t, pressureFromAlt(alt(t)))
            ax5.set_title("Pression en fonction du temps")
            ax5.set_xlabel("Temps (s)")
            ax5.set_ylabel("Pression (Pa)")
            ax5.grid(True, which="both")

            # Temperature

            fig6, ax6 = plt.subplots(figsize=figsize)

            ax6.plot(t, temperatureFromAlt(alt(t)))
            ax6.set_title("Température en fonction du temps")
            ax6.set_xlabel("Temps (s)")
            ax6.set_ylabel("Température (K)")
            ax6.grid(True, which="both")

            # Density

            fig7, ax7 = plt.subplots(figsize=figsize)

            ax7.plot(t, densityFromAlt(alt(t)))
            ax7.set_title("Densité en fonction du temps")
            ax7.set_xlabel("Temps (s)")
            ax7.set_ylabel("Densité (kg/m\u00B3)")
            ax7.grid(True, which="both")

            if self.parent.shape is None and self.parent.shapes != []:
                C_NList = []
                tList = []
                fxList = []
                fyList = []
                pfxList = []
                pfyList = []

                # testtimeList = {}
                # testTList = {}
                for h in t:
                    air = Zone(vel(h) / soundSpeedFromAlt(alt(h)),
                               inc(h) * np.pi / 180, z=alt(h))
                    try:
                        z, C_N, fx, fy, pfx, pfy = calculateGeomZones(
                            air, self.parent.method, self.parent.angleBack,
                            self.parent.conv, self.parent.points,
                            self.parent.vertex, self.parent.geometry,
                            self.parent.dataD)
                    except Exception as e:
                        print(h, "s", e)
                    else:
                        if None not in [C_N, fx, fy, pfx, pfy]:
                            tList.append(h)
                            C_NList.append(C_N)
                            fxList.append(fx)
                            fyList.append(fy)
                            pfxList.append(pfx)
                            pfyList.append(pfy)

                        # for k in z:
                        #     if k is not None:
                        #         i = k.i
                        #         if i in testtimeList:
                        #             if testtimeList[i][-1] == h:
                        #                 testTList[i][-1] = max(
                        #                     testTList[i][-1], k.pressure)
                        #             else:
                        #                 testtimeList[i].append(h)
                        #                 testTList[i].append(k.pressure)
                        #         else:
                        #             testtimeList[i] = [h]
                        #             testTList[i] = [k.pressure]
                        print(h, "s", C_N)

                # figt, axt = plt.subplots(figsize=figsize)
                # label = ["0", "Pointe", "Tour de sauvetage", "Capsule",
                #          "Étage 4", "Inter Étage 3-4", "Étage 3",
                #          "Inter Étage 2-3", "Étage 1 et 2", ""]

                # for k in testtimeList.keys():
                #     axt.plot(testtimeList[k], testTList[k], label=label[k])
                # axt.legend()
                # # axt.set_title("Température maximum sur "
                # #               "la fusée en fonction du temps")
                # axt.set_title("Pression maximum sur "
                #               "la fusée en fonction du temps")
                # axt.set_xlabel("Temps (s)")
                # # axt.set_ylabel("Température (K)")
                # axt.set_ylabel("Pression (Pa)")
                # axt.grid(True, which="both")

                fig8, ax8 = plt.subplots(figsize=figsize)

                ax8.plot(tList, C_NList)
                ax8.set_title(r"$C_N$ en fonction du temps")
                ax8.set_xlabel("Temps (s)")
                ax8.set_ylabel(r"$C_N$")
                ax8.grid(True, which="both")

                tm = self.parent.dataT["tm"]
                xs = self.parent.dataT["xs"]
                ys = self.parent.dataT["ys"]

                fxList = np.array(fxList) * tm**2 / xs
                fyList = np.array(fyList) * tm**2 / ys

                fig9, ax9 = plt.subplots(figsize=figsize)

                ax9.plot(tList, fxList, label="Selon X")
                ax9.plot(tList, fyList, label="Selon Y")
                ax9.set_title("Forces appliquées en fonction du temps")
                ax9.legend()
                ax9.set_xlabel("Temps (s)")
                ax9.set_ylabel("Force (N)")
                ax9.grid(True, which="both")

                fileCG = tkfd.askopenfile(mode="r", filetypes=[
                    ("CSV Files", "*.csv"), ("Any File", "*.*")],
                    title="CSV file containing Time (s) and CG position (in).")
                if fileCG is not None:
                    t, x, y = readTable(fileCG)

                    fx = interp(tList, fxList)
                    fy = interp(tList, fyList)
                    pfx = interp(tList, pfxList)
                    pfy = interp(tList, pfyList)
                    cgx = interp(t, x * 0.0254)  # in -> m
                    cgy = interp(t, y * 0.0254)  # in -> m

                    tmini = max(min(t), min(tList))
                    tmaxi = min(max(t), max(tList))
                    t = np.linspace(tmini, tmaxi, 500)

                    m = (pfx(t) - cgx(t)) * fx(t) + (pfy(t) - cgy(t)) * fy(t)

                    fig0, ax0 = plt.subplots(figsize=figsize)

                    ax0.plot(t, m)
                    ax0.set_title("Moments appliquées en fonction du temps")
                    ax0.set_xlabel("Temps (s)")
                    ax0.set_ylabel("Moment (N.m)")
                    ax0.grid(True, which="both")

            plt.show()

    def updateM(self):
        # TODO: Move to MainApplication
        m = self.method.get()
        if m == 1:
            self.parent.method = getZoneC
        elif m == 2:
            self.parent.method = getZoneTTP
        elif m == 3:
            self.parent.method = getZoneDT
        self.parent.updateZones()


class UserInput(tk.Frame):
    def __init__(self, parent, label, unit, default,
                 fro, to, update, log=False, fun=None, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.default = default
        self.update = update
        # TODO: Simplify by using log function as fun
        # def logfun(x, inv=False):
        #     if inv:
        #         return np.log10(x)
        #     return 10**float(x)
        self.log = log  # If the scale is logarithmic
        self.fun = fun  # If a lookup table is used (inv=True for inverse fun)

        self.label = tk.Label(self, text=label)
        self.label.config(font="TkFixedFont")  # Monospace label
        self.unit = tk.Label(self, text=unit)
        self.unit.config(font="TkFixedFont")

        self.var = tk.DoubleVar(value=default)
        self.spinbox = tk.Spinbox(self, from_=fro, to=to, width=8,
                                  textvariable=self.var, command=lambda:
                                  self.updateValue(self.var.get()))
        self.spinbox.bind('<Return>', lambda _:
                          self.updateValue(self.var.get()))

        self.scale = tk.Scale(self, from_=np.log10(fro) if log else fro,
                              to=np.log10(to) if log else to,
                              orient=tk.HORIZONTAL, showvalue=False,
                              resolution=(((np.log10(to) if log else to)
                                           - (np.log10(fro) if log else fro))
                                          if fun is None else 1)
                              / 1e6, command=lambda x: self.updateValue(
            (10**float(x) if log else float(x))
            if self.fun is None else self.fun(x)))
        self.scale.set(np.log10(default) if log else default)

        self.label.pack(side=tk.LEFT)
        self.unit.pack(side=tk.RIGHT)
        self.spinbox.pack(side=tk.RIGHT)
        self.scale.pack(fill=tk.X, expand=True)

    def updateScale(self, fro, to):
        self.scale.configure(from_=fro, to=to, resolution=(
            (
                (np.log10(to) if self.log else to)
                - (np.log10(fro) if self.log else fro)
            ) / 1e6 if self.fun is None else 1))

    def updateValue(self, value, update=True):
        value = round_sig(value, 6)
        self.var.set(value)
        self.scale.set((np.log10(value) if self.log else value)
                       if self.fun is None else self.fun(value, inv=True))
        if update:
            self.update(value)


class MainApplication(tk.Frame):
    def __init__(self, parent, air, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        parent.title("Aérodynamique hypersonique")
        self.parent = parent
        self.air = air

        self.menubar = Menubar(self)
        parent.config(menu=self.menubar)

        self.camera = Camera()
        self.canvas = Canvas(self, self.camera)

        self.angleBack = 0
        self.conv = 0
        self.atmoAppAlt = UserInput(self, "Altitude ", "km", air.altitude,
                                    0, 100, self.userUpdateZ)
        self.userInputs = [  # Order reversed
            UserInput(self, "C Back   ", "° ", self.conv,
                      -45, 45, self.userUpdateC),
            UserInput(self, "Δ Back   ", "° ", self.angleBack,
                      -90, 90, self.userUpdateB),
            self.atmoAppAlt,
            UserInput(self, "Angle    ", "° ", air.angle,
                      -90, 90, self.userUpdateA),
            UserInput(self, "Mach     ", "  ", air.mach,
                      1, 100, self.userUpdateM, True)
        ]

        self.canvas.pack(fill=tk.BOTH, expand=True)
        for u in self.userInputs:
            u.pack(side=tk.BOTTOM, fill=tk.X)

        self.zoneAngle = self.canvas.createVector(30, 40, 40, self.air.angle)
        self.userText = [self.canvas.create_text(
            10, 10, text="M\u2081 = " + str(np.round(air.mach, 4))
            + " (" + str(np.round(air.mach
                                  * soundSpeedFromAlt(air.altitude), 1))
            + " m/s at " + str(np.round(air.altitude, 4)) + " km)",
            font="Consolas", anchor=tk.W)]

        self.mouseX, self.mouseY = 0, 0
        self.mouseUX, self.mouseUY = 0, 0

        self.method = getZoneC

        self.points = []
        self.vertex = []
        self.geometry = []

        self.shape = None  # TODO: Merge into same system
        self.shapes = []
        self.shapeText = None
        self.lines = []
        self.text_sig = 6
        self.texts = []
        self.canvas.bind("<Button-1>", self.mouseDown)
        self.canvas.bind("<B1-Motion>", self.mouseDrag)
        self.canvas.bind("<ButtonRelease-1>", self.mouseUp)
        self.canvas.bind("<Configure>", self.resizeCanvas)

        self.canvas.bind("<MouseWheel>", self.mouseWheel)
        self.mouseW = 0
        self.canvas.bind("<B2-Motion>", self.mouseDragMiddle)
        self.mouse2X, self.mouse2Y = 0, 0

        self.dataT = None
        self.dataD = None
        if len(sys.argv) == 2 and sys.argv[1] is not None:
            with open(sys.argv[1], "r") as f:
                data = js.load(f)
                if "transform" in data:
                    self.dataT = data["transform"]
                else:
                    self.dataT = {"rot": 0, "xs": 1, "ys": 1,
                                  "x0": 0, "y0": 0, "tm": 1}
                self.dataP = data["points"]
                self.dataV = data["vertex"]
                self.dataG = data["geometry"]
                if "depth" in data:
                    self.dataD = data["depth"]
                self.loadGeom()

    def mouseWheel(self, event):
        self.mouseW -= event.delta / 120
        self.camera.updateZoom(np.exp(self.mouseW / 100))

    def resizeCanvas(self, _):
        self.updateZones()

    def mouseDragMiddle(self, event):
        self.camera.move(event.x - self.mouse2X, event.y - self.mouse2Y)
        self.mouse2X, self.mouse2Y = event.x, event.y

    def userUpdateM(self, mach):
        self.air.mach = mach
        self.updateMachText()

        if self.updateZones():
            self.updateShapeText()

    def userUpdateA(self, angle):
        self.air.angle = np.pi * angle / 180
        self.canvas.updateVector(self.zoneAngle, 30, 40, 40, self.air.angle)
        if self.updateZones():
            self.updateShapeText()

    def userUpdateC(self, angle):
        self.conv = np.pi * angle / 180
        self.updateZones()

    def userUpdateB(self, angle):
        self.angleBack = np.pi * angle / 180
        self.updateZones()

    def userUpdateT(self, t):
        self.air.temperature = t
        self.updateMachText()
        self.updateZones()

    def userUpdateP(self, p):
        self.air.pressure = p
        # self.updateMachText()
        self.updateZones()

    def userUpdateD(self, d):
        self.air.density = d
        # self.updateMachText()
        self.updateZones()

    def userUpdateZ(self, z):
        self.air.updateAlt(z)
        self.updateMachText()
        self.updateZones()

    def userUpdateSig(self, sig):
        self.text_sig = sig
        # self.updateMachText()
        self.updateZones()

    def invertVertex(self):
        for i, v in enumerate(self.vertex):
            self.vertex[i] = (v[0], v[1], -v[2])
        self.updateZones()

    def updateMachText(self):
        self.canvas.itemconfig(self.userText[0],
                               text="M\u2081 = " + str(np.round(air.mach, 4))
                               + " (" + str(np.round(air.mach
                                                     * soundSpeedFromTemp(
                                                         air.temperature), 1))
                               + " m/s at " + str(np.round(air.altitude, 4))
                               + " km)")

    def updateShapeText(self):
        theta = np.arctan2(self.mouseY - self.mouseUY,
                           self.mouseUX - self.mouseX)  # rad
        maxAngle = maximumDeflectionAngle(self.air.mach)
        color = "red"
        if self.air.angle - theta <= maxAngle \
                and self.air.angle >= theta - maxAngle:
            color = "black"

        self.canvas.itemconfig(self.shapeText, text="θ = " + str(np.round(
            (theta - self.air.angle) * 180 / np.pi, 2)) + "°", fill=color)

    def updateZones(self, forced=False):
        if not self.menubar.autoUpdate.get() and not forced:
            return False
        if self.air.mach < 1:
            return False

        if self.shape is None and self.shapes != []:
            zones, C_N, fx, fy, pfx, pfy = calculateGeomZones(
                self.air, self.method, self.angleBack, self.conv,
                self.points, self.vertex, self.geometry, self.dataD)
        elif self.shape is not None and self.shapes == []:
            zones, C_N, fx, fy, pfx, pfy = calculateZones(
                self.air, self.method, self.angleBack,
                *self.camera.screenToSpace(self.mouseX, self.mouseY),
                *self.camera.screenToSpace(self.mouseUX, self.mouseUY))
        else:
            return False

        # TODO: Temp fix (Remove None):
        zones = [z for z in zones if z is not None]

        self.canvas.delete(*self.lines, *self.texts)
        self.lines = []
        self.texts = []

        if C_N is not None and C_N < 1e6:
            tm = 0.01 if self.dataT is None else self.dataT["tm"]
            xs = 1 if self.dataT is None else self.dataT["xs"]
            ys = 1 if self.dataT is None else self.dataT["ys"]

            # px*u*Pa ----> * (m/u)**2 / (px/u) ----> m*m*Pa = N
            fx *= tm**2 / xs  # TODO: Not correct when xs != yx
            fy *= tm**2 / ys

            print("C_N =", round_sig(C_N, self.text_sig), "fx =",
                  round(fx, self.text_sig), "N fy =",
                  round(fy, self.text_sig), "N")

            if self.dataT is not None:
                m = np.arctan(np.sqrt(fx**2 + fy**2) / 30) * 30
                a = np.arctan2(fx, fy)

                self.texts.append(self.canvas.createVectorSpace(
                    pfx, pfy, m, a, c="#57E", h=0))

        for z in zones:
            self.texts += z.displayData(self.canvas, self.text_sig)
            for x, y, a, c in z.lines:
                self.lines.append(self.canvas.createLine(x, y, a, c))

        return True

    def resetCanvas(self):
        self.canvas.delete(self.shape, self.shapeText)
        self.canvas.delete(*self.shapes, *self.lines, *self.texts)

        self.shapes = []
        self.lines = []
        self.texts = []
        self.shapeText = self.shape = None

        self.points = []
        self.vertex = []
        self.geometry = []

    def newGeom(self):
        self.dataT = self.dataP = self.dataV = self.dataG = self.dataD = None
        self.resetCanvas()
        self.menubar.disableTr()

    def transformPoints(self, points):
        a = -np.pi * self.dataT["rot"] / 180
        xs = self.dataT["xs"]
        ys = self.dataT["ys"]
        x0 = self.dataT["x0"]
        y0 = self.dataT["y0"]
        R = np.array([[np.cos(a) * xs, -np.sin(a) * xs],
                      [np.sin(a) * ys, np.cos(a) * ys]])

        return np.round(np.dot(points, R)
                        + [[x0, y0]] * len(points), 9).tolist()

    def loadGeom(self):
        self.resetCanvas()
        self.menubar.enableTr()

        self.points = self.transformPoints(self.dataP)
        self.vertex = self.dataV
        self.geometry = self.dataG

        for p1, p2, _ in self.vertex:
            self.shapes.append(self.canvas.create_line(
                *self.camera.spaceToScreen(*self.points[p1]),
                *self.camera.spaceToScreen(*self.points[p2])))

        self.updateZones()

    def openGeom(self):
        file = tkfd.askopenfile(mode="r", filetypes=[
            ("JSON Files", "*.json"), ("Any File", "*.*")])
        if file is not None:
            data = js.load(file)
            file.close()
            if "transform" in data:
                self.dataT = data["transform"]
            else:
                self.dataT = {"rot": 0, "xs": 1, "ys": 1,
                              "x0": 0, "y0": 0, "tm": 1}
            self.dataP = data["points"]
            self.dataV = data["vertex"]
            self.dataG = data["geometry"]
            if "depth" in data:
                self.dataD = data["depth"]
            self.loadGeom()

    def saveGeom(self):
        file = tkfd.asksaveasfile(mode="w", filetypes=[
            ("JSON Files", "*.json"), ("Any File", "*.*")],
            defaultextension=".json")
        if file is not None:
            if self.shape is None and self.shapes != []:
                data = {
                    "transform": self.dataT,
                    "points": self.dataP,
                    "vertex": self.dataV,
                    "geometry": self.dataG,
                }
                if self.dataD is not None:
                    data["depth"] = self.dataD
                js.dump(data, file)
            elif self.shape is not None and self.shapes == []:
                x1, y1, x2, y2 = self.canvas.coords(self.shape)
                js.dump({
                    "points": [self.camera.screenToSpace(x1, y1),
                               self.camera.screenToSpace(x2, y2)],
                    "vertex": [(0, 1, 1), (0, 1, -1)],
                    "geometry": [[0], [1]],
                }, file)
            else:
                raise "Could not save data, invalid internal format."
            file.close()

    def exampleGeom(self, k):
        if k == 1:
            self.dataT = {"rot": 0, "xs": 1, "ys": 1,
                          "x0": 400, "y0": -250, "tm": 0.01}
            self.dataP = [(-200, 0), (0, -10), (0, 10), (200, 0)]
            self.dataV = [(0, 1, -1), (1, 3, -1), (0, 2, 1), (2, 3, 1)]
            self.dataG = [[0, 1], [2, 3]]
            self.dataD = None
        elif k == 2:
            self.dataT = {"rot": 0, "xs": 1, "ys": 1,
                          "x0": 200, "y0": -260, "tm": 0.01}
            self.dataP = [[0, 0], [300, -40], [150, 10], [450, -30]]
            self.dataV = [[0, 1, 1], [0, 1, -1], [2, 3, 1], [2, 3, -1]]
            self.dataG = [[0, [[2], [3]]], [1]]
            self.dataD = None
        else:
            print(k)
            return

        self.loadGeom()

    def openAtmo(self):
        t = tk.Toplevel(self)
        t.grab_set()
        AtmoApplication(t, self, self.atmoAppAlt).pack(
            side=tk.TOP, fill=tk.BOTH, expand=True)

    def openTr(self):
        if self.dataT is not None:
            t = tk.Toplevel(self)
            t.grab_set()
            TrApplication(t, self).pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def openCG(self):
        if self.dataT is not None:
            t = tk.Toplevel(self)
            t.grab_set()
            CGApplication(t, self).pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def resetInputs(self):
        for u in self.userInputs:
            u.updateValue(u.default)

    def mouseDown(self, event):
        self.mouseX, self.mouseY = event.x, event.y

        self.newGeom()

        self.shape = self.canvas.create_line(
            self.mouseX, self.mouseY, event.x, event.y)
        self.shapeText = self.canvas.create_text(
            self.canvas.winfo_width() / 2, 10, text="θ = 0°")

        self.updateZones(True)

        self.canvas.tag_raise(self.shape)

    def mouseDrag(self, event):
        self.mouseUX, self.mouseUY = event.x, event.y

        self.canvas.coords(
            self.shape, self.mouseX, self.mouseY, self.mouseUX, self.mouseUY)

        self.updateZones(True)
        self.updateMachText()
        self.updateShapeText()

    def mouseUp(self, event):
        self.mouseUX, self.mouseUY = event.x, event.y
        if self.shape is not None and self.shapes == []:
            x, y = self.camera.screenToSpace(self.mouseX, self.mouseY)
            dx, dy = self.camera.screenToSpace(self.mouseUX, self.mouseUY)
            dx -= x
            dy -= y
            self.dataT = {"rot": 0, "xs": 1, "ys": 1,
                          "x0": x, "y0": y, "tm": 0.01}
            self.dataP = [[0, 0], [dx, dy]]
            self.dataV = [[0, 1, 1], [0, 1, -1]]
            self.dataG = [[0], [1]]
            self.dataD = None
            self.loadGeom()

    def onExit(self):
        self.parent.destroy()


if __name__ == "__main__":
    air = Zone(2, 0, z=10)

    root = tk.Tk()
    MainApplication(root, air).pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    root.mainloop()
