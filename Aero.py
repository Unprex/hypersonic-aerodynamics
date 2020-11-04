# -*- coding: utf-8 -*-
import tkinter as tk
import numpy as np

from Tables import maximumDeflectionAngle, shockAngleFromDeflection


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


class Camera:
    """ Handles coordinate transform between workspace and screen """

    def __init__(self, screen, x=0, y=0, z=1):
        self.screen = screen
        self.x, self.y, self.z = x, y, z

    def createLine(self, x, y, a, c):
        x1, y1 = self.spaceToScreen(x, y)
        x2, y2 = pointToBorder(x1, y1, self.angleToScreen(a),
                               self.screen.winfo_width(),
                               self.screen.winfo_height())
        return self.screen.create_line(x1, y1, x2, y2, fill=c)

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

    if theta > front.angle and theta < front.angle + maxAngle:
        # If the deflection is "up" regarding the stream with attached shock

        # Compression
        sAngle = shockAngleFromDeflection(theta - front.angle, front.mach)
        MachC = front.mach * np.cos(sAngle) \
            / np.cos(front.angle + sAngle - theta)
        zCompr = Zone(MachC, theta, front)
        zCompr.addIncidentAngle(x1, y1, -sAngle - front.angle, "red")

        # Expansion
        MachLine1 = np.arcsin(1 / front.mach)
        MachE = front.mach  # TODO : Prandlt-Meyer
        MachLine2 = np.arcsin(1 / MachE)
        zExp = Zone(MachE, theta, front)
        zExp.addIncidentAngle(x1, y1, MachLine1 - front.angle, "green")
        zExp.addIncidentAngle(x1, y1, MachLine2 - theta, "lime")

        return [zCompr, zExp]

    elif theta < front.angle and theta > front.angle - maxAngle:
        # If the deflection is "down" regarding the stream with attached shock
        # print("OkD", -theta)
        # Compression
        sAngle = shockAngleFromDeflection(-theta + front.angle, front.mach)
        MachC = front.mach * np.cos(sAngle) \
            / np.cos(front.angle + sAngle + theta)
        zCompr = Zone(MachC, -theta, front)
        zCompr.addIncidentAngle(x1, y1, sAngle - front.angle, "red")

        # Expansion
        MachLine1 = np.arcsin(1 / front.mach)
        MachE = front.mach  # TODO : Prandlt-Meyer
        MachLine2 = np.arcsin(1 / MachE)
        zExp = Zone(MachE, -theta, front)
        zExp.addIncidentAngle(x1, y1, -MachLine1 - front.angle, "green")
        zExp.addIncidentAngle(x1, y1, -MachLine2 - theta, "lime")

        return [zCompr, zExp]

    return []


class Canvas(tk.Canvas):
    def __init__(self, parent, *args, **kwargs):
        tk.Canvas.__init__(self, parent, *args, **kwargs,
                           bd=0, highlightthickness=0)
        self.parent = parent


class MainApplication(tk.Frame):
    def __init__(self, parent, air, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.air = air
        self.canvas = Canvas(self)

        self.canvas.pack(fill=tk.BOTH, expand=True)

        self.canvas.create_line(10, 10, 70, 10, arrow=tk.LAST)
        self.canvas.create_text(40, 15,
                                text="M\u2081 = " + str(round(air.mach, 2)),
                                font="Consolas", anchor=tk.N)

        self.camera = Camera(self.canvas)

        self.mouseX, self.mouseY = 0, 0
        self.shape, self.text = None, []
        self.lines = [None] * 6
        self.linesDelete = []
        self.canvas.bind("<Button-1>", self.mouseDown)
        self.canvas.bind("<B1-Motion>", self.mouseDrag)
        self.canvas.bind("<ButtonRelease-1>", self.mouseUp)

    def createText(self, x, y, text, anchor=tk.W, color="black"):
        self.text.append(self.canvas.create_text(x, y, text=text,
                                                 anchor=anchor, fill=color))
        return len(self.text) - 1

    def updateText(self, tId, text, color="black"):
        self.canvas.itemconfig(self.text[tId], text=text, fill=color)

    def mouseDown(self, event):
        self.mouseX, self.mouseY = event.x, event.y

        zones = calculateZones(
            self.air,
            *self.camera.screenToSpace(self.mouseX, self.mouseY),
            *self.camera.screenToSpace(event.x, event.y))

        self.canvas.delete(self.shape, *self.text, *self.lines,
                           *self.linesDelete)
        self.text = []

        for z in zones:
            for x, y, a, c in z.lines:
                self.linesDelete.append(self.camera.createLine(x, y, a, c))

        self.shape = self.canvas.create_line(
            self.mouseX, self.mouseY, event.x, event.y)
        self.createText(80, 10, "θ = 0°")

        self.canvas.tag_raise(self.shape)

    def mouseDrag(self, event):
        self.canvas.coords(self.shape,
                           self.mouseX, self.mouseY, event.x, event.y)

        zones = calculateZones(
            self.air,
            *self.camera.screenToSpace(self.mouseX, self.mouseY),
            *self.camera.screenToSpace(event.x, event.y))

        self.canvas.delete(*self.linesDelete)

        for z in zones:
            for x, y, a, c in z.lines:
                self.linesDelete.append(self.camera.createLine(x, y, a, c))

        M = self.air.mach

        theta = np.arctan2(self.mouseY - event.y, event.x - self.mouseX)  # rad
        maxAngle = maximumDeflectionAngle(M)
        color = "red"
        if self.air.angle - theta <= maxAngle \
                and self.air.angle >= theta - maxAngle:
            color = "black"

        self.updateText(0, "θ = " + str(round(
            (theta - self.air.angle) * 180 / np.pi, 2)) + "°", color)

    def mouseUp(self, event):
        # Create aero surface and start calculations
        pass  # TODO: updateCalculations()


if __name__ == "__main__":
    air = Zone(2, 0)  # Air at mach 2

    root = tk.Tk()
    MainApplication(root, air).pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    root.mainloop()
