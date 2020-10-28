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


def calculateZones(front, x1, y1, x2, y2):
    theta = np.arctan2(y1 - y2, x1 - x2)  # rad

    maxAngle = maximumDeflectionAngle(front.mach)

    if theta > front.angle:
        print("Ok", theta)


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

        self.mouseX, self.mouseY = 0, 0
        self.shape, self.text = None, None
        self.lines = [None] * 6
        self.canvas.bind("<Button-1>", self.mouseDown)
        self.canvas.bind("<B1-Motion>", self.mouseDrag)
        self.canvas.bind("<ButtonRelease-1>", self.mouseUp)

    def mouseDown(self, event):
        self.mouseX, self.mouseY = event.x, event.y

        calculateZones(self.air, self.mouseX, self.mouseY, event.x, event.y)
        M = self.air.mach

        self.canvas.delete(self.shape, self.text, *self.lines)
        self.shape = self.canvas.create_line(
            self.mouseX, self.mouseY, event.x, event.y)
        self.text = self.canvas.create_text(80, 10, text="θ = 0°", anchor=tk.W)

        theta = 0
        self.lines[0] = self.canvas.create_line(
            self.mouseX, self.mouseY, pointToBorder(
                self.mouseX, self.mouseY, np.arcsin(1/M),
                self.winfo_width(), self.winfo_height()), fill="green")

        self.lines[1] = self.canvas.create_line(
            self.mouseX, self.mouseY, pointToBorder(
                self.mouseX, self.mouseY, np.arcsin(1/M) + theta,
                self.winfo_width(), self.winfo_height()), fill="lime")

        self.lines[2] = self.canvas.create_line(
            self.mouseX, self.mouseY, pointToBorder(
                self.mouseX, self.mouseY, -shockAngleFromDeflection(-theta, M),
                self.winfo_width(), self.winfo_height()), fill="red")

        self.lines[3] = self.canvas.create_line(
            event.x, event.y, pointToBorder(
                event.x, event.y, shockAngleFromDeflection(-theta, M) + theta,
                self.winfo_width(), self.winfo_height()), fill="red")

        self.lines[4] = self.canvas.create_line(
            event.x, event.y, pointToBorder(
                event.x, event.y, -np.arcsin(1/M) + theta,
                self.winfo_width(), self.winfo_height()), fill="green")

        self.lines[5] = self.canvas.create_line(
            event.x, event.y, pointToBorder(
                event.x, event.y, -np.arcsin(1/M),
                self.winfo_width(), self.winfo_height()), fill="lime")

        self.canvas.tag_raise(self.shape)

    def mouseDrag(self, event):
        self.canvas.coords(self.shape,
                           self.mouseX, self.mouseY, event.x, event.y)

        calculateZones(self.air, self.mouseX, self.mouseY, event.x, event.y)
        M = self.air.mach

        theta = np.arctan2(self.mouseY - event.y, event.x - self.mouseX)  # rad
        maxAngle = maximumDeflectionAngle(M)
        color = "red"
        if -theta <= maxAngle and -theta >= 0:
            color = "black"

            self.canvas.coords(self.lines[0], self.mouseX, self.mouseY,
                               *pointToBorder(
                                   self.mouseX, self.mouseY,
                                   np.arcsin(1/M),
                                   self.winfo_width(), self.winfo_height()))

            self.canvas.coords(self.lines[1], self.mouseX, self.mouseY,
                               *pointToBorder(
                                   self.mouseX, self.mouseY,
                                   np.arcsin(1/M) + theta,
                                   self.winfo_width(), self.winfo_height()))

            self.canvas.coords(self.lines[2], self.mouseX, self.mouseY,
                               *pointToBorder(
                                   self.mouseX, self.mouseY,
                                   -shockAngleFromDeflection(-theta, M),
                                   self.winfo_width(), self.winfo_height()))

            self.canvas.coords(self.lines[3], event.x, event.y, *pointToBorder(
                event.x, event.y, shockAngleFromDeflection(-theta, M) + theta,
                self.winfo_width(), self.winfo_height()))

            self.canvas.coords(self.lines[4], event.x, event.y, *pointToBorder(
                event.x, event.y, -np.arcsin(1/M) + theta,
                self.winfo_width(), self.winfo_height()))

            self.canvas.coords(self.lines[5], event.x, event.y, *pointToBorder(
                event.x, event.y, -np.arcsin(1/M),
                self.winfo_width(), self.winfo_height()))

            for line in self.lines:
                self.canvas.itemconfigure(line, state="normal")
        else:
            for line in self.lines:
                self.canvas.itemconfigure(line, state="hidden")

        self.canvas.itemconfig(
            self.text, text="θ = " + str(round(theta * 180 / np.pi, 2)) + "°",
            fill=color)

    def mouseUp(self, event):
        # Create aero surface and start calculations
        pass  # TODO: updateCalculations()


if __name__ == "__main__":
    air = Zone(2, 0)  # Air at mach 2

    root = tk.Tk()
    MainApplication(root, air).pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    root.mainloop()
