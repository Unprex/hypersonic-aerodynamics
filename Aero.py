# -*- coding: utf-8 -*-
import tkinter as tk
import tkinter.filedialog as tkfd
import numpy as np
import json as js
import sys

from Tables import (
    maximumDeflectionAngle, shockAngleFromDeflection, machWaveAngle, maximumNu,
    prandtlMeyerMachFromNu, prandtlMeyerNuFromMach, prandltMeyerMuFromMach,
    p_sur_pisentropique, T_sur_Tisentropique, rho_sur_rhoisentropique,
    p2surp1fromMachSh, rho2surrho1fromMachSh, t2surt1fromMachSh,
    pressureFromAlt, temperatureFromAlt, densityFromAlt, soundSpeedFromAlt,
    cpFromM0p2surp1, mach2frommach1, soundSpeedFromTemp
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
    return np.round(x, sig - np.int(np.floor(np.log10(np.abs(x) + 1e-6))) - 1)


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

        self.lines = []

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
        x = self.lines[0][0]
        y = self.lines[0][1] + (70 if self.lines[0][2] > self.angle else -70)
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

    zUp1 = getZone(front, theta, 1, x1, y1)
    zDown1 = getZone(front, theta, -1, x1, y1)
    if zUp1 is not None and zDown1 is not None:
        C_N = zDown1.Cp - zUp1.Cp
    else:
        C_N = None

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

    return [zUp1, zDown1, zUp2, zDown2], C_N


def calculateGeomZones(front, getZone, angleBack, points, vertex, geometry):
    newAngle = front.angle + angleBack

    def sequence(seq, prec):
        zones = []
        for s in seq:
            if type(s) == list:
                zones += sequence(s, prec)
            else:
                p1, p2, up = vertex[s]
                p1, p2 = points[p1], points[p2]
                theta = np.arctan2(p2[1] - p1[1], p2[0] - p1[0])  # rad
                zone = getZone(prec, theta, up, *p1)
                if zone is None:
                    return zones
                zones.append(zone)
                prec = zone
        if type(s) != list:
            _, p2, up = vertex[s]
            lastZone = getZone(prec, newAngle, up, *points[p2])
            if lastZone is not None:
                lastZone.addIncidentAngle(*points[p2], newAngle, "gray")
            zones.append(lastZone)
        return zones

    return sequence(geometry, front), None


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

    def createVector(self, x, y, n, a, c=None):
        a = self.camera.angleToScreen(a)
        return self.create_line(
            x - n * np.cos(a) / 2, y - n * np.sin(a) / 2,
            x + n * np.cos(a) / 2, y + n * np.sin(a) / 2,
            arrow=tk.LAST, fill=c)

    def createVectorSpace(self, x, y, n, a, c=None):
        return self.createVector(*self.camera.spaceToScreen(x, y), n, a, c)

    def updateVector(self, vid, x, y, n, a):
        a = self.camera.angleToScreen(a)
        self.coords(vid, x - n * np.cos(a) / 2, y - n * np.sin(a) / 2,
                    x + n * np.cos(a) / 2, y + n * np.sin(a) / 2)

    def updateVectorSpace(self, vid, x, y, n, a):
        self.updateVector(vid, *self.camera.spaceToScreen(x, y), n, a)


class AtmoApplication(tk.Frame):
    def __init__(self, parent, mApp, uInput, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        parent.title("Sélection atmosphère")
        parent.geometry("500x200")
        self.parent = parent
        self.mApp = mApp
        self.uInput = uInput

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
        tk.Spinbox(sigmas, from_=1, to=25, width=8, textvariable=self.sig,
                   command=self.updateSig).pack(side=tk.LEFT, pady=5)

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
        self.uInput[1].updateValue(z, False)
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
        self.uInput[1].updateValue(z, False)
        self.mApp.userUpdateZ(z)
        self.mApp.userUpdateD(d)
        self.mApp.userUpdateT(t)
        self.mApp.userUpdateP(p)
        self.mApp.userUpdateSig(self.sig_default)
        self.parent.destroy()


class Menubar(tk.Menu):
    def __init__(self, parent, *args, **kwargs):
        tk.Menu.__init__(self, parent, *args, **kwargs)
        self.parent = parent

        fileMenu = tk.Menu(self, tearoff=0)
        fileMenu.add_command(label="New", command=parent.resetCanvas)
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
        toolsMenu.add_separator()
        toolsMenu.add_command(label="Reset inputs", command=parent.resetInputs)

        var = tk.IntVar(value=1)
        methodMenu = tk.Menu(self, tearoff=0)
        methodMenu.add_radiobutton(label="Clasical",
                                   variable=var, value=1, command=self.updateM)
        methodMenu.add_radiobutton(label="Small-Disturbance Theory",
                                   variable=var, value=2, command=self.updateM)
        methodMenu.add_radiobutton(label="Tangent Wedge",
                                   variable=var, value=3, command=self.updateM)
        self.method = var

        self.add_cascade(label="File", menu=fileMenu)
        self.add_cascade(label="Tools", menu=toolsMenu)
        self.add_cascade(label="Method", menu=methodMenu)

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
                 fro, to, update, log=False, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.default = default
        self.update = update
        self.log = log  # If the scale is logarithmic

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
                              resolution=((np.log10(to) if log else to)
                                          - (np.log10(fro) if log else fro))
                              / 1e6, command=lambda x: self.updateValue(
            10**float(x) if log else float(x)))
        self.scale.set(np.log10(default) if log else default)

        self.label.pack(side=tk.LEFT)
        self.unit.pack(side=tk.RIGHT)
        self.spinbox.pack(side=tk.RIGHT)
        self.scale.pack(fill=tk.X, expand=True)

    def updateValue(self, value, update=True):
        value = round_sig(value, 6)
        self.var.set(value)
        self.scale.set(np.log10(value) if self.log else value)
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
        self.userInputs = [  # Order reversed
            UserInput(self, "Δ Back   ", "° ", self.angleBack,
                      -90, 90, self.userUpdateB),
            UserInput(self, "Altitude ", "km", air.altitude,
                      0, 100, self.userUpdateZ),
            UserInput(self, "Angle    ", "° ", air.angle,
                      -90, 90, self.userUpdateA),
            UserInput(self, "Mach     ", "  ", air.mach,
                      1, 100, self.userUpdateM, True)
        ]  # If updated, also update AtmoApplication

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

        self.canvas.bind("<MouseWheel>", self.mouseWheel)
        self.mouseW = 0
        self.canvas.bind("<B2-Motion>", self.mouseDragMiddle)
        self.mouse2X, self.mouse2Y = 0, 0

        if len(sys.argv) == 2 and sys.argv[1] is not None:
            with open(sys.argv[1], "r") as f:
                data = js.load(f)
                self.loadGeom(data["points"], data["vertex"], data["geometry"])

    def mouseWheel(self, event):
        self.mouseW -= event.delta / 120
        self.camera.updateZoom(np.exp(self.mouseW / 100))

    def mouseDragMiddle(self, event):
        self.camera.move(event.x - self.mouse2X, event.y - self.mouse2Y)
        self.mouse2X, self.mouse2Y = event.x, event.y

    def userUpdateM(self, mach):
        self.air.mach = mach
        self.updateMachText()

        if not self.updateZones():
            self.updateShapeText()

    def userUpdateA(self, angle):
        self.air.angle = np.pi * angle / 180
        self.canvas.updateVector(self.zoneAngle, 30, 40, 40, self.air.angle)
        if self.updateZones():
            self.updateShapeText()

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

    def updateZones(self):
        if self.shape is None and self.shapes != []:
            zones, C_N = calculateGeomZones(
                self.air, self.method, self.angleBack,
                self.points, self.vertex, self.geometry)
        elif self.shape is not None and self.shapes == []:
            zones, C_N = calculateZones(
                self.air, self.method, self.angleBack,
                *self.camera.screenToSpace(self.mouseX, self.mouseY),
                *self.camera.screenToSpace(self.mouseUX, self.mouseUY))
        else:
            return False

        if C_N is not None and C_N < 1e6:
            print("C_N:", round_sig(C_N, self.text_sig))

        # TODO: Temp fix (Remove None):
        zones = [z for z in zones if z is not None]

        self.canvas.delete(*self.lines, *self.texts)
        self.lines = []
        self.texts = []

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

    def loadGeom(self, points, vertex, geometry):
        self.resetCanvas()

        self.points = points
        self.vertex = vertex
        self.geometry = geometry

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
            self.loadGeom(data["points"], data["vertex"], data["geometry"])

    def saveGeom(self):
        file = tkfd.asksaveasfile(mode="w", filetypes=[
            ("JSON Files", "*.json"), ("Any File", "*.*")],
            defaultextension=".json")
        if file is not None:
            if self.points != [] and self.vertex != [] and self.geometry != []:
                js.dump({
                    "points": self.points,
                    "vertex": self.vertex,
                    "geometry": self.geometry,
                }, file)
            elif self.shape is not None:
                x1, y1, x2, y2 = self.canvas.coords(self.shape)
                js.dump({
                    "points": [self.camera.screenToSpace(x1, y1),
                               self.camera.screenToSpace(x2, y2)],
                    "vertex": [(0, 1, 1), (0, 1, -1)],
                    "geometry": [[0], [1]],
                }, file)
            file.close()

    def exampleGeom(self, k):
        if k == 1:
            points = [(200, -250), (400, -260), (400, -240), (600, -250)]
            vertex = [(0, 1, -1), (1, 3, -1), (0, 2, 1), (2, 3, 1)]
            geometry = [[0, 1], [2, 3]]
        elif k == 2:
            points = [[200, -260], [500, -300], [350, -250], [650, -290]]
            vertex = [[0, 1, 1], [0, 1, -1], [2, 3, 1], [2, 3, -1]]
            geometry = [[0, [[2], [3]]], [1]]
        else:
            print(k)
            return

        self.loadGeom(points, vertex, geometry)

    def openAtmo(self):
        t = tk.Toplevel(self)
        t.grab_set()
        AtmoApplication(t, self, self.userInputs).pack(
            side=tk.TOP, fill=tk.BOTH, expand=True)

    def resetInputs(self):
        for u in self.userInputs:
            u.updateValue(u.default)

    def mouseDown(self, event):
        self.mouseX, self.mouseY = event.x, event.y

        self.resetCanvas()

        self.shape = self.canvas.create_line(
            self.mouseX, self.mouseY, event.x, event.y)
        self.shapeText = self.canvas.create_text(
            self.canvas.winfo_width() / 2, 10, text="θ = 0°")

        self.updateZones()

        self.canvas.tag_raise(self.shape)

    def mouseDrag(self, event):
        self.mouseUX, self.mouseUY = event.x, event.y

        self.canvas.coords(
            self.shape, self.mouseX, self.mouseY, self.mouseUX, self.mouseUY)

        self.updateZones()
        self.updateMachText()
        self.updateShapeText()

    def mouseUp(self, event):
        self.mouseUX, self.mouseUY = event.x, event.y
        # TODO: updateCalculations()

    def onExit(self):
        self.parent.destroy()


if __name__ == "__main__":
    air = Zone(2, 0, z=10)

    root = tk.Tk()
    MainApplication(root, air).pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    root.mainloop()
