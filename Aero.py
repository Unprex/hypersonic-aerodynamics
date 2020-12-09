# -*- coding: utf-8 -*-
import tkinter as tk
import numpy as np

from Tables import (
    maximumDeflectionAngle, shockAngleFromDeflection, machWaveAngle, maximumNu,
    prandtlMeyerMachFromNu, prandtlMeyerNuFromMach, prandltMeyerMuFromMach

)


def pointToBorder(x, y, angle, width, height):  # angle in rad
    # Get line from point and angle to border of screen

    angle %= 2 * np.pi
    if angle < np.pi / 4 or angle > 7 * np.pi / 4:  # Right (x = width)
        return width, np.tan(angle) * (x - width) + y
    elif angle < 3 * np.pi / 4:  # Top (y = 0)
        return x + y * np.cos(angle) / np.sin(angle), 0
    elif angle < 5 * np.pi / 4:  # Left (x = 0)
        return 0, np.tan(angle) * x + y
    else:  # Bottom (y = height)
        return x - (height - y) * np.cos(angle) / np.sin(angle), height


class Zone:
    def __init__(self, mach, angle, parent=None):
        self.mach = mach
        self.angle = angle
        self.parent = parent

        self.lines = []

    def addIncidentAngle(self, x, y, angle, color):
        self.lines.append((x, y, angle, color))

    def displayData(self, canvas):
        text = "M = %f" % np.round(self.mach, 4)
        x = self.lines[0][0]
        y = self.lines[0][1] - np.sign(self.lines[0][2]) * 25
        c = self.lines[0][3]

        return [canvas.createTextSpace(x, y + 10, text, fill=c),
                canvas.createVectorSpace(x, y - 10, 40, self.angle, c)]


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


def calculateZones(front, x1, y1, x2, y2):
    theta = np.arctan2(y2 - y1, x2 - x1)  # rad

    # Get maximum attached angle
    maxAngle = maximumDeflectionAngle(front.mach)

    # def getZoneCompression(theta):
    #     sAngle = shockAngleFromDeflection(theta, M)

    # def getZoneExpansion():
    #     pass
    dAngle = theta - front.angle

    if dAngle > 0 and theta < front.angle + maxAngle:
        # If the deflection is "up" regarding the stream with attached shock

        # Compression
        sAngle = shockAngleFromDeflection(dAngle, front.mach)
        MachC = front.mach * np.cos(sAngle) / np.cos(sAngle - dAngle)

        zCompr = Zone(MachC, theta, front)
        zCompr.addIncidentAngle(x1, y1, -sAngle - front.angle, "red")

        # Expansion
        MachLine1 = np.arcsin(1 / front.mach)

        newNu = prandtlMeyerNuFromMach(front.mach) + dAngle
        assert newNu <= maximumNu()
        MachE = prandtlMeyerMachFromNu(newNu)

        MachLine2 = np.arcsin(1 / MachE)
        zExp = Zone(MachE, theta, front)
        zExp.addIncidentAngle(x1, y1, MachLine1 - front.angle, "green")
        zExp.addIncidentAngle(x1, y1, MachLine2 - theta, "lime")

        return [zCompr, zExp]

    elif dAngle < 0 and dAngle > -maxAngle:
        # If the deflection is "down" regarding the stream with attached shock
        # print("OkD", -theta)
        # Compression
        sAngle = shockAngleFromDeflection(-dAngle, front.mach)
        MachC = front.mach * np.cos(sAngle) / np.cos(sAngle - dAngle)
        zCompr = Zone(MachC, theta, front)
        zCompr.addIncidentAngle(x1, y1, sAngle - front.angle, "red")

        # Expansion
        MachLine1 = prandltMeyerMuFromMach(front.mach)

        newNu = prandtlMeyerNuFromMach(front.mach) - dAngle
        assert newNu <= maximumNu()
        MachE = prandtlMeyerMachFromNu(newNu)

        MachLine2 = prandltMeyerMuFromMach(MachE)
        zExp = Zone(MachE, theta, front)
        zExp.addIncidentAngle(x1, y1, -MachLine1 - front.angle, "green")
        zExp.addIncidentAngle(x1, y1, -MachLine2 - theta, "lime")

        return [zCompr, zExp]

    elif dAngle == 0:
        sAngle = machWaveAngle(front.mach)

        zUp = Zone(front.mach, theta, front)
        zUp.addIncidentAngle(x1, y1, -sAngle - front.angle, "red")
        zDo = Zone(front.mach, theta, front)
        zDo.addIncidentAngle(x1, y1, sAngle - front.angle, "red")
        return [zUp, zDo]

    return []


class Canvas(tk.Canvas):
    def __init__(self, parent, camera, *args, **kwargs):
        tk.Canvas.__init__(self, parent, *args, **kwargs,
                           bd=0, highlightthickness=0)
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


class UserInput(tk.Frame):
    def __init__(self, parent, label, default,
                 from_, to, update, log=False, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.update = update
        self.log = log  # If the scale is logarithmic

        self.label = tk.Label(self, text=label)

        self.var = tk.DoubleVar(value=default)
        self.spinbox = tk.Spinbox(self, from_=from_, to=to, width=8,
                                  textvariable=self.var, command=lambda:
                                  self.updateValue(self.var.get()))
        self.spinbox.bind('<Return>', lambda _:
                          self.updateValue(self.var.get()))

        self.scale = tk.Scale(self, from_=np.log10(from_) if log else from_,
                              to=np.log10(to) if log else to,
                              orient=tk.HORIZONTAL, showvalue=False,
                              resolution=(np.log10(to) if log else to
                                          - np.log10(from_) if log else from_)
                              / 1e6, command=lambda x: self.updateValue(
            10**float(x) if log else float(x)))
        self.scale.set(np.log10(default) if log else default)

        self.label.pack(side=tk.LEFT)
        self.spinbox.pack(side=tk.RIGHT)
        self.scale.pack(fill=tk.X, expand=True)

    def updateValue(self, value):
        value = np.round(value, 4)
        self.var.set(value)
        self.scale.set(np.log10(value) if self.log else value)
        self.update(value)


class MainApplication(tk.Frame):
    def __init__(self, parent, air, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.air = air

        self.camera = Camera()
        self.canvas = Canvas(self, self.camera)

        self.userMach = UserInput(
            self, "Mach  ", 2, 1, 100, self.userUpdateM, True)
        self.userAngle = UserInput(
            self, "Angle  ", 0, -90, 90, self.userUpdateA)

        self.canvas.pack(fill=tk.BOTH, expand=True)
        self.userAngle.pack(side=tk.BOTTOM, fill=tk.X)
        self.userMach.pack(side=tk.BOTTOM, fill=tk.X)

        self.zoneAngle = self.canvas.createVector(30, 40, 40, self.air.angle)
        self.userText = [self.canvas.create_text(
            10, 10, text="M\u2081 = " + str(np.round(air.mach, 4)),
            font="Consolas", anchor=tk.W)]

        self.mouseX, self.mouseY = 0, 0
        self.mouseUX, self.mouseUY = 0, 0
        self.shape = None
        self.shapeText = None
        self.lines = []
        self.texts = []
        self.canvas.bind("<Button-1>", self.mouseDown)
        self.canvas.bind("<B1-Motion>", self.mouseDrag)
        self.canvas.bind("<ButtonRelease-1>", self.mouseUp)

        self.canvas.bind("<MouseWheel>", self.mouseWheel)
        self.mouseW = 0
        self.canvas.bind("<B2-Motion>", self.mouseDragMiddle)
        self.mouse2X, self.mouse2Y = 0, 0

    def mouseWheel(self, event):
        self.mouseW -= event.delta / 120
        self.camera.updateZoom(np.exp(self.mouseW / 100))

    def mouseDragMiddle(self, event):
        self.camera.move(event.x - self.mouse2X, event.y - self.mouse2Y)
        self.mouse2X, self.mouse2Y = event.x, event.y

    def userUpdateM(self, mach):
        self.air.mach = mach

        self.canvas.itemconfig(self.userText[0],
                               text="M\u2081 = " + str(np.round(air.mach, 4)))
        if self.shape is None:
            return

        self.updateZones()

        theta = np.arctan2(self.mouseY - self.mouseUY,
                           self.mouseUX - self.mouseX)  # rad
        maxAngle = maximumDeflectionAngle(self.air.mach)
        color = "red"
        if self.air.angle - theta <= maxAngle \
                and self.air.angle >= theta - maxAngle:
            color = "black"

        self.canvas.itemconfig(self.shapeText, text="θ = " + str(np.round(
            (theta - self.air.angle) * 180 / np.pi, 2)) + "°", fill=color)

    def userUpdateA(self, angle):
        self.air.angle = np.pi * angle / 180

        self.canvas.itemconfig(self.userText[0],
                               text="M\u2081 = " + str(np.round(air.mach, 4)))

        self.canvas.updateVector(self.zoneAngle, 30, 40, 40, self.air.angle)
        if self.shape is None:
            return

        self.updateZones()

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
        zones = calculateZones(
            self.air,
            *self.camera.screenToSpace(self.mouseX, self.mouseY),
            *self.camera.screenToSpace(self.mouseUX, self.mouseUY))

        self.canvas.delete(*self.lines, *self.texts)
        self.lines = []
        self.texts = []

        for z in zones:
            self.texts += z.displayData(self.canvas)
            for x, y, a, c in z.lines:
                self.lines.append(self.canvas.createLine(x, y, a, c))

    def mouseDown(self, event):
        self.mouseX, self.mouseY = event.x, event.y

        self.canvas.delete(self.shape, self.shapeText)
        self.shapeText = None

        self.updateZones()

        self.shape = self.canvas.create_line(
            self.mouseX, self.mouseY, event.x, event.y)
        self.shapeText = self.canvas.create_text(
            self.canvas.winfo_width() / 2, 10, text="θ = 0°")

        self.canvas.tag_raise(self.shape)

    def mouseDrag(self, event):
        self.mouseUX, self.mouseUY = event.x, event.y

        self.canvas.coords(
            self.shape, self.mouseX, self.mouseY, self.mouseUX, self.mouseUY)

        self.updateZones()

        M = self.air.mach

        theta = np.arctan2(self.mouseY - event.y, event.x - self.mouseX)  # rad
        maxAngle = maximumDeflectionAngle(M)
        color = "red"
        if self.air.angle - theta <= maxAngle \
                and self.air.angle >= theta - maxAngle:
            color = "black"

        self.canvas.itemconfig(self.userText[0],
                               text="M\u2081 = " + str(np.round(air.mach, 4)))

        self.canvas.itemconfig(self.shapeText, text="θ = " + str(np.round(
            (theta - self.air.angle) * 180 / np.pi, 2)) + "°", fill=color)

    def mouseUp(self, event):
        self.mouseUX, self.mouseUY = event.x, event.y
        # TODO: updateCalculations()


if __name__ == "__main__":
    air = Zone(2, 0)  # Air at mach 2

    root = tk.Tk()
    MainApplication(root, air).pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    root.mainloop()
