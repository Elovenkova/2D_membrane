from math import*
import tkinter
import random
import numpy as np
from numpy.linalg import inv

WIDTH = 1540.
HEIGHT = 1480.
BG_COLOR = 'white'
DELAY = 1.
SPEED_X = 0
SPEED_Y = 0
SPEED_Z = 0
X_0 = 270
Y_0 = 40
Z_0 = 0.
SELF_COLOR = 'black'
CX = (1/3)**0.5
CY = (1/3)**0.5
CZ = (1/3)**0.5
ANGLE = 0
LENGTH = 0.1
PIXEL = 10000
DIST = 0.0005
H = 0.0005
DIAMETER = 0.0001
E = 125000000000
DENSITY = 1430
PRESSURE = 5000000000
FORCE = PRESSURE*DIAMETER*DIST
MASS = DENSITY*(0.25*pi*(DIAMETER**2)*DIST)
S_SPEED = (E/DENSITY)**0.5
POINTS_PER_SIZE = int(LENGTH/DIST)
RA = LA = UA = DA = 1

class Point():
    def __init__(self, x, y, z, color, dx, dy, dz, cx, cy, cz, ang, ra, la, da, ua):
        self.x = x
        self.y = y
        self.z = z
        self.color = color
        self.dx = dx
        self.dy = dy	
        self.dz = dz
        self.cx = cx
        self.cy = cy
        self.cz = cz
        self.ang = ang
        self.la = la
        self.ra = ra
        self.da = da
        self.ua = ua
    def matrix(self, d, x0, y0, z0): 
        T = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [-x0, -y0, -z0, 1]])
        Rx = np.array([[1, 0, 0, 0], [0, (self.cz/d), (self.cy/d), 0], [0, -(self.cy/d), (self.cz/d), 0], [0, 0, 0, 1]])
        Ry = np.array([[d, 0, self.cx, 0], [0, 1, 0, 0], [-self.cx, 0, d, 0], [0, 0, 0, 1]])
        R0 = np.array([[cos(self.ang), sin(self.ang), 0, 0], [-sin(self.ang), cos(self.ang), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]) 
        T_inv = inv(T)
        Rx_inv = inv(Rx)
        Ry_inv = inv(Ry)
        M = np.dot(T, Rx)
        M = np.dot(M, Ry)
        M = np.dot(M, R0)
        M = np.dot(M, Ry_inv)
        M = np.dot(M, Rx_inv)
        M = np.dot(M, T_inv)
        return M

    def draw(self):
        m = self.matrix((self.cy**2 + self.cz**2)**0.5, 200, 100, 100)
        canvas.create_oval(self.x*m[0][0] + self.y*m[1][0] + self.z*m[2][0] + m[3][0] - 1, self.x*m[0][1] + self.y*m[1][1] + self.z*m[2][1] + m[3][1] - 1, self.x*m[0][0] + self.y*m[1][0] + self.z*m[2][0] + m[3][0] + 1, self.x*m[0][1] + self.y*m[1][1] + self.z*m[2][1] + m[3][1] + 1, fill=self.color, outline=self.color)

    def hide(self):
        m = self.matrix((self.cy**2 + self.cz**2)**0.5, 200, 100, 100)
        canvas.create_oval(self.x*m[0][0] + self.y*m[1][0] + self.z*m[2][0] + m[3][0] - 1, self.x*m[0][1] + self.y*m[1][1] + self.z*m[2][1] + m[3][1] - 1, self.x*m[0][0] + self.y*m[1][0] + self.z*m[2][0] + m[3][0] + 1, self.x*m[0][1] + self.y*m[1][1] + self.z*m[2][1] + m[3][1] + 1, fill=BG_COLOR, outline=BG_COLOR)

    def move(self):
        self.hide()
        self.x += self.dx*DELAY
        self.y += self.dy*DELAY
        self.z += self.dz*DELAY
        self.draw()

def P(i, j, d):
    if (((i-POINTS_PER_SIZE//2)**2+(j-POINTS_PER_SIZE//2)**2)**0.5 <= POINTS_PER_SIZE//5):
        return (FORCE/MASS)*d*exp(-((i-POINTS_PER_SIZE//2)**2+(j-POINTS_PER_SIZE//2)**2)/30)
    else:
        return 0

def main():

    global DELAY, VIS_DELAY

    def motion(event):
        if event.keysym == "R1ight":	     
            for i in range(0, POINTS_PER_SIZE):
                for j in range(0, POINTS_PER_SIZE):
                    points[i][j].ang += 0.03
                    canvas.create_rectangle(0, 0, WIDTH, HEIGHT, fill = BG_COLOR)
        elif event.keysym == "Left":	     
            for i in range(0, POINTS_PER_SIZE):
                for j in range(0, POINTS_PER_SIZE):
                    points[i][j].ang -= 0.03
                    canvas.create_rectangle(0, 0, WIDTH, HEIGHT, fill = BG_COLOR)
 
    for i in range(1, POINTS_PER_SIZE):
        for j in range(1, POINTS_PER_SIZE):
             m1 = points[i][j].matrix((points[i][j].cy**2 + points[i][j].cz**2)**0.5, 200, 100, 100)
             m2 = points[i-1][j].matrix((points[i-1][j].cy**2 + points[i-1][j].cz**2)**0.5, 200, 100, 100)
             m3 = points[i][j-1].matrix((points[i][j-1].cy**2 + points[i][j-1].cz**2)**0.5, 200, 100, 100)
             canvas.create_line(points[i][j].x*m1[0][0] + points[i][j].y*m1[1][0] + points[i][j].z*m1[2][0] + m1[3][0], points[i][j].x*m1[0][1] + points[i][j].y*m1[1][1] + points[i][j].z*m1[2][1] + m1[3][1], points[i-1][j].x*m2[0][0] + points[i-1][j].y*m2[1][0] + points[i-1][j].z*m2[2][0] + m2[3][0], points[i-1][j].x*m2[0][1] + points[i-1][j].y*m2[1][1] + points[i-1][j].z*m2[2][1] + m2[3][1], fill = BG_COLOR, width = 3)
             canvas.create_line(points[i][j].x*m1[0][0] + points[i][j].y*m1[1][0] + points[i][j].z*m1[2][0] + m1[3][0], points[i][j].x*m1[0][1] + points[i][j].y*m1[1][1] + points[i][j].z*m1[2][1] + m1[3][1], points[i][j-1].x*m3[0][0] + points[i][j-1].y*m3[1][0] + points[i][j-1].z*m3[2][0] + m3[3][0], points[i][j-1].x*m3[0][1] + points[i][j-1].y*m3[1][1] + points[i][j-1].z*m3[2][1] + m3[3][1], fill = BG_COLOR, width = 3)

    for j in range(1, POINTS_PER_SIZE):
        m1 = points[0][j].matrix((points[0][j].cy**2 + points[0][j].cz**2)**0.5, 200, 100, 100)
        m3 = points[0][j-1].matrix((points[0][j-1].cy**2 + points[0][j-1].cz**2)**0.5, 200, 100, 100)
        canvas.create_line(points[0][j].x*m1[0][0] + points[0][j].y*m1[1][0] + points[0][j].z*m1[2][0] + m1[3][0], points[0][j].x*m1[0][1] + points[0][j].y*m1[1][1] + points[0][j].z*m1[2][1] + m1[3][1], points[0][j-1].x*m3[0][0] + points[0][j-1].y*m3[1][0] + points[0][j-1].z*m3[2][0] + m3[3][0], points[0][j-1].x*m3[0][1] + points[0][j-1].y*m3[1][1] + points[0][j-1].z*m3[2][1] + m3[3][1], fill = BG_COLOR, width = 3)

    for i in range(1, POINTS_PER_SIZE):
        m1 = points[i][0].matrix((points[i][0].cy**2 + points[i][0].cz**2)**0.5, 200, 100, 100)
        m2 = points[i-1][0].matrix((points[i-1][0].cy**2 + points[i-1][0].cz**2)**0.5, 200, 100, 100)
        canvas.create_line(points[i][0].x*m1[0][0] + points[i][0].y*m1[1][0] + points[i][0].z*m1[2][0] + m1[3][0], points[i][0].x*m1[0][1] + points[i][0].y*m1[1][1] + points[i][0].z*m1[2][1] + m1[3][1], points[i-1][0].x*m2[0][0] + points[i-1][0].y*m2[1][0] + points[i-1][0].z*m2[2][0] + m2[3][0], points[i-1][0].x*m2[0][1] + points[i-1][0].y*m2[1][1] + points[i-1][0].z*m2[2][1] + m2[3][1], fill = BG_COLOR, width = 3)
    
    for i in range(0, POINTS_PER_SIZE):
        for j in range(0, POINTS_PER_SIZE):
             points[i][j].move()

    for i in range(1, POINTS_PER_SIZE):
        if ((((points[0][i-1].x - points[0][i].x))**2 + ((points[0][i-1].y - points[0][i].y))**2 + (points[0][i-1].z - points[0][i].z)**2)**0.5 != 0):
            D[0][i] = (((points[0][i-1].x - points[0][i].x))**2 + ((points[0][i-1].y - points[0][i].y))**2 + (points[0][i-1].z - points[0][i].z)**2)**0.5/S_SPEED
    for i in range(1, POINTS_PER_SIZE):
        if ((((points[i-1][0].x - points[i][0].x))**2 + ((points[i-1][0].y - points[i][0].y))**2 + (points[i-1][0].z - points[i][0].z)**2)**0.5 != 0):
            D[j][0] = (((points[i-1][0].x - points[i][0].x))**2 + ((points[i-1][0].y - points[i][0].y))**2 + (points[i-1][0].z - points[i][0].z)**2)**0.5/S_SPEED
    for i in range(1, POINTS_PER_SIZE):
        for j in range(1, POINTS_PER_SIZE):
            d_l = (((points[i][j-1].x - points[i][j].x))**2 + ((points[i][j-1].y - points[i][j].y))**2 + (points[i][j-1].z - points[i][j].z)**2)**0.5 
            d_u = (((points[i-1][j].x - points[i][j].x))**2 + ((points[i-1][j].y - points[i][j].y))**2 + (points[i-1][j].z - points[i][j].z)**2)**0.5
            if (d_l != 0 and d_u != 0):
                D[i][j] = min(d_l, d_u)/S_SPEED

    for i in range(0, POINTS_PER_SIZE):
        D_S[i] = min(D[i])
    DELAY = min(D_S)

    for i in range(0, POINTS_PER_SIZE):
        for j in range(0, POINTS_PER_SIZE):
            points_new[i][j].dz += P(i, j, DELAY)

    for i in range(0, POINTS_PER_SIZE):
        for j in range(0, POINTS_PER_SIZE):
            if (j > 0 and j < POINTS_PER_SIZE - 1 and i > 0 and i < POINTS_PER_SIZE - 1):
                d_l = (((points[i][j-1].x - points[i][j].x))**2 + ((points[i][j-1].y - points[i][j].y))**2 + (points[i][j-1].z - points[i][j].z)**2)**0.5 
                d_r = (((points[i][j+1].x - points[i][j].x))**2 + ((points[i][j+1].y - points[i][j].y))**2 + (points[i][j+1].z - points[i][j].z)**2)**0.5 
                d_d = (((points[i+1][j].x - points[i][j].x))**2 + ((points[i+1][j].y - points[i][j].y))**2 + (points[i+1][j].z - points[i][j].z)**2)**0.5 
                d_u = (((points[i-1][j].x - points[i][j].x))**2 + ((points[i-1][j].y - points[i][j].y))**2 + (points[i-1][j].z - points[i][j].z)**2)**0.5
                e_l = 1 - H/d_l
                e_r = 1 - H/d_r
                e_d = 1 - H/d_d
                e_u = 1 - H/d_u
                l_z = (points[i][j-1].z - points[i][j].z)/d_l
                l_y = (points[i][j-1].y - points[i][j].y)/(d_l*(1-l_z**2)**0.5)
                l_x = (points[i][j-1].x - points[i][j].x)/(d_l*(1-l_z**2)**0.5) 
                r_z = (points[i][j+1].z - points[i][j].z)/d_r
                r_y = (points[i][j+1].y - points[i][j].y)/(d_r*(1-r_z**2)**0.5)
                r_x = (points[i][j+1].x - points[i][j].x)/(d_r*(1-r_z**2)**0.5)
                d_z = (points[i+1][j].z - points[i][j].z)/d_d
                d_y = (points[i+1][j].y - points[i][j].y)/(d_d*(1-d_z**2)**0.5)
                d_x = (points[i+1][j].x - points[i][j].x)/(d_d*(1-d_z**2)**0.5)
                u_z = (points[i-1][j].z - points[i][j].z)/d_u
                u_y = (points[i-1][j].y - points[i][j].y)/(d_u*(1-u_z**2)**0.5)
                u_x = (points[i-1][j].x - points[i][j].x)/(d_u*(1-u_z**2)**0.5)
                if (e_l >= 0.4):
                    points_new[i][j].la = 0
                    e_l = 0
                if (e_r >= 0.4):
                    points_new[i][j].ra = 0
                    e_r = 0
                if (e_d >= 0.4):
                    points_new[i][j].da = 0
                    e_d = 0
                if (e_u >= 0.4):
                    points_new[i][j].ua = 0
                    e_u = 0
                points_new[i][j].dx += (e_l*l_x + e_r*r_x + e_d*d_x + e_u*u_x)*E*DELAY*DIAMETER**2*pi/(4*MASS) 
                points_new[i][j].dy += (e_l*l_y + e_r*r_y + e_d*d_y + e_u*u_y)*E*DELAY*DIAMETER**2*pi/(4*MASS)
                points_new[i][j].dz += (e_l*l_z + e_r*r_z + e_d*d_z + e_u*u_z)*E*DELAY*DIAMETER**2*pi/(4*MASS) 
            elif (i == 0 and j == 0): 
                d_r = (((points[i][j+1].x - points[i][j].x))**2 + ((points[i][j+1].y - points[i][j].y))**2 + (points[i][j+1].z - points[i][j].z)**2)**0.5 
                d_d = (((points[i+1][j].x - points[i][j].x))**2 + ((points[i+1][j].y - points[i][j].y))**2 + (points[i+1][j].z - points[i][j].z)**2)**0.5 
                e_r = 1 - H/d_r
                e_d = 1 - H/d_d
                r_z = (points[i][j+1].z - points[i][j].z)/d_r
                r_y = (points[i][j+1].y - points[i][j].y)/(d_r*(1-r_z**2)**0.5)
                r_x = (points[i][j+1].x - points[i][j].x)/(d_r*(1-r_z**2)**0.5)
                d_z = (points[i+1][j].z - points[i][j].z)/d_d
                d_y = (points[i+1][j].y - points[i][j].y)/(d_d*(1-d_z**2)**0.5)
                d_x = (points[i+1][j].x - points[i][j].x)/(d_d*(1-d_z**2)**0.5)
                if (e_r >= 0.4):
                    points_new[i][j].ra = 0
                    e_r = 0
                if (e_d >= 0.4):
                    points_new[i][j].da = 0
                    e_d = 0
                points_new[i][j].dx += (e_r*r_x + e_d*d_x)*E*DELAY*DIAMETER**2*pi/(4*MASS) 
                points_new[i][j].dy += (e_r*r_y + e_d*d_y)*E*DELAY*DIAMETER**2*pi/(4*MASS)
                points_new[i][j].dz += (e_r*r_z + e_d*d_z)*E*DELAY*DIAMETER**2*pi/(4*MASS)
            elif (j == POINTS_PER_SIZE - 1 and i == 0): 
                d_l = (((points[i][j-1].x - points[i][j].x))**2 + ((points[i][j-1].y - points[i][j].y))**2 + (points[i][j-1].z - points[i][j].z)**2)**0.5 
                d_d = (((points[i+1][j].x - points[i][j].x))**2 + ((points[i+1][j].y - points[i][j].y))**2 + (points[i+1][j].z - points[i][j].z)**2)**0.5 
                e_l = 1 - H/d_l
                e_d = 1 - H/d_d
                l_z = (points[i][j-1].z - points[i][j].z)/d_l
                l_y = (points[i][j-1].y - points[i][j].y)/(d_l*(1-l_z**2)**0.5)
                l_x = (points[i][j-1].x - points[i][j].x)/(d_l*(1-l_z**2)**0.5)
                d_z = (points[i+1][j].z - points[i][j].z)/d_d
                d_y = (points[i+1][j].y - points[i][j].y)/(d_d*(1-d_z**2)**0.5)
                d_x = (points[i+1][j].x - points[i][j].x)/(d_d*(1-d_z**2)**0.5)
                if (e_l >= 0.4):
                    points_new[i][j].la = 0
                    e_l = 0
                if (e_d >= 0.4):
                    points_new[i][j].da = 0
                    e_d = 0
                points_new[i][j].dx += (e_l*l_x + e_d*d_x)*E*DELAY*DIAMETER**2*pi/(4*MASS)
                points_new[i][j].dy += (e_l*l_y + e_d*d_y)*E*DELAY*DIAMETER**2*pi/(4*MASS)
                points_new[i][j].dz += (e_l*l_z + e_d*d_z)*E*DELAY*DIAMETER**2*pi/(4*MASS)
            elif (j == 0 and i == POINTS_PER_SIZE - 1): 
                d_r = (((points[i][j+1].x - points[i][j].x))**2 + ((points[i][j+1].y - points[i][j].y))**2 + (points[i][j+1].z - points[i][j].z)**2)**0.5 
                d_u = (((points[i-1][j].x - points[i][j].x))**2 + ((points[i-1][j].y - points[i][j].y))**2 + (points[i-1][j].z - points[i][j].z)**2)**0.5 
                e_r = 1 - H/d_r
                e_u = 1 - H/d_u
                r_z = (points[i][j+1].z - points[i][j].z)/d_r
                r_y = (points[i][j+1].y - points[i][j].y)/(d_r*(1-r_z**2)**0.5)
                r_x = (points[i][j+1].x - points[i][j].x)/(d_r*(1-r_z**2)**0.5)
                u_z = (points[i-1][j].z - points[i][j].z)/d_u
                u_y = (points[i-1][j].y - points[i][j].y)/(d_u*(1-u_z**2)**0.5)
                u_x = (points[i-1][j].x - points[i][j].x)/(d_u*(1-u_z**2)**0.5)
                if (e_r >= 0.4):
                    points_new[i][j].ra = 0
                    e_r = 0
                if (e_u >= 0.4):
                    points_new[i][j].ua = 0
                    e_u = 0
                points_new[i][j].dx += (e_r*r_x + e_u*u_x)*E*DELAY*DIAMETER**2*pi/(4*MASS) 
                points_new[i][j].dy += (e_r*r_y + e_u*u_y)*E*DELAY*DIAMETER**2*pi/(4*MASS)
                points_new[i][j].dz += (e_r*r_z + e_u*u_z)*E*DELAY*DIAMETER**2*pi/(4*MASS)
            elif (j == POINTS_PER_SIZE - 1 and i == POINTS_PER_SIZE - 1): 
                d_l = (((points[i][j-1].x - points[i][j].x))**2 + ((points[i][j-1].y - points[i][j].y))**2 + (points[i][j-1].z - points[i][j].z)**2)**0.5 
                d_u = (((points[i-1][j].x - points[i][j].x))**2 + ((points[i-1][j].y - points[i][j].y))**2 + (points[i-1][j].z - points[i][j].z)**2)**0.5 
                e_l = 1 - H/d_l
                e_u = 1 - H/d_u
                l_z = (points[i][j-1].z - points[i][j].z)/d_l
                l_y = (points[i][j-1].y - points[i][j].y)/(d_l*(1-l_z**2)**0.5)
                l_x = (points[i][j-1].x - points[i][j].x)/(d_l*(1-l_z**2)**0.5)
                u_z = (points[i-1][j].z - points[i][j].z)/d_u
                u_y = (points[i-1][j].y - points[i][j].y)/(d_u*(1-u_z**2)**0.5)
                u_x = (points[i-1][j].x - points[i][j].x)/(d_u*(1-u_z**2)**0.5)
                if (e_l >= 0.4):
                    points_new[i][j].la = 0
                    e_l = 0
                if (e_u >= 0.4):
                    points_new[i][j].ua = 0
                    e_u = 0
                points_new[i][j].dx += (e_l*l_x + e_u*u_x)*E*DELAY*DIAMETER**2*pi/(4*MASS) 
                points_new[i][j].dy += (e_l*l_y + e_u*u_y)*E*DELAY*DIAMETER**2*pi/(4*MASS)
                points_new[i][j].dz += (e_l*l_z + e_u*u_z)*E*DELAY*DIAMETER**2*pi/(4*MASS)
            elif (i == 0 and j < POINTS_PER_SIZE - 1 and j > 0):
                d_l = (((points[i][j-1].x - points[i][j].x))**2 + ((points[i][j-1].y - points[i][j].y))**2 + (points[i][j-1].z - points[i][j].z)**2)**0.5 
                d_r = (((points[i][j+1].x - points[i][j].x))**2 + ((points[i][j+1].y - points[i][j].y))**2 + (points[i][j+1].z - points[i][j].z)**2)**0.5 
                d_d = (((points[i+1][j].x - points[i][j].x))**2 + ((points[i+1][j].y - points[i][j].y))**2 + (points[i+1][j].z - points[i][j].z)**2)**0.5 
                e_l = 1 - H/d_l
                e_r = 1 - H/d_r
                e_d = 1 - H/d_d
                l_z = (points[i][j-1].z - points[i][j].z)/d_l
                l_y = (points[i][j-1].y - points[i][j].y)/(d_l*(1-l_z**2)**0.5)
                l_x = (points[i][j-1].x - points[i][j].x)/(d_l*(1-l_z**2)**0.5) 
                r_z = (points[i][j+1].z - points[i][j].z)/d_r
                r_y = (points[i][j+1].y - points[i][j].y)/(d_r*(1-r_z**2)**0.5)
                r_x = (points[i][j+1].x - points[i][j].x)/(d_r*(1-r_z**2)**0.5)
                d_z = (points[i+1][j].z - points[i][j].z)/d_d
                d_y = (points[i+1][j].y - points[i][j].y)/(d_d*(1-d_z**2)**0.5)
                d_x = (points[i+1][j].x - points[i][j].x)/(d_d*(1-d_z**2)**0.5)
                if (e_l >= 0.4):
                    points_new[i][j].la = 0
                    e_l = 0
                if (e_r >= 0.4):
                    points_new[i][j].ra = 0
                    e_r = 0
                if (e_d >= 0.4):
                    points_new[i][j].da = 0
                    e_d = 0
                points_new[i][j].dx += (e_l*l_x + e_r*r_x + e_d*d_x)*E*DELAY*DIAMETER**2*pi/(4*MASS)  
                points_new[i][j].dy += (e_l*l_y + e_r*r_y + e_d*d_y)*E*DELAY*DIAMETER**2*pi/(4*MASS)
                points_new[i][j].dz += (e_l*l_z + e_r*r_z + e_d*d_z)*E*DELAY*DIAMETER**2*pi/(4*MASS)
            elif (i == POINTS_PER_SIZE - 1 and j < POINTS_PER_SIZE - 1 and j > 0): 
                d_l = (((points[i][j-1].x - points[i][j].x))**2 + ((points[i][j-1].y - points[i][j].y))**2 + (points[i][j-1].z - points[i][j].z)**2)**0.5 
                d_r = (((points[i][j+1].x - points[i][j].x))**2 + ((points[i][j+1].y - points[i][j].y))**2 + (points[i][j+1].z - points[i][j].z)**2)**0.5 
                d_u = (((points[i-1][j].x - points[i][j].x))**2 + ((points[i-1][j].y - points[i][j].y))**2 + (points[i-1][j].z - points[i][j].z)**2)**0.5
                e_l = 1 - H/d_l
                e_r = 1 - H/d_r
                e_u = 1 - H/d_u
                l_z = (points[i][j-1].z - points[i][j].z)/d_l
                l_y = (points[i][j-1].y - points[i][j].y)/(d_l*(1-l_z**2)**0.5)
                l_x = (points[i][j-1].x - points[i][j].x)/(d_l*(1-l_z**2)**0.5) 
                r_z = (points[i][j+1].z - points[i][j].z)/d_r
                r_y = (points[i][j+1].y - points[i][j].y)/(d_r*(1-r_z**2)**0.5)
                r_x = (points[i][j+1].x - points[i][j].x)/(d_r*(1-r_z**2)**0.5)
                u_z = (points[i-1][j].z - points[i][j].z)/d_u
                u_y = (points[i-1][j].y - points[i][j].y)/(d_u*(1-u_z**2)**0.5)
                u_x = (points[i-1][j].x - points[i][j].x)/(d_u*(1-u_z**2)**0.5)
                if (e_l >= 0.4):
                    points_new[i][j].la = 0
                    e_l = 0
                if (e_r >= 0.4):
                    points_new[i][j].ra = 0
                    e_r = 0
                if (e_u >= 0.4):
                    points_new[i][j].ua = 0
                    e_u = 0
                points_new[i][j].dx += (e_l*l_x + e_r*r_x + e_u*u_x)*E*DELAY*DIAMETER**2*pi/(4*MASS) 
                points_new[i][j].dy += (e_l*l_y + e_r*r_y + e_u*u_y)*E*DELAY*DIAMETER**2*pi/(4*MASS)
                points_new[i][j].dz += (e_l*l_z + e_r*r_z + e_u*u_z)*E*DELAY*DIAMETER**2*pi/(4*MASS)
            elif (j == 0 and i > 0 and i < POINTS_PER_SIZE - 1):
                d_r = (((points[i][j+1].x - points[i][j].x))**2 + ((points[i][j+1].y - points[i][j].y))**2 + (points[i][j+1].z - points[i][j].z)**2)**0.5 
                d_d = (((points[i+1][j].x - points[i][j].x))**2 + ((points[i+1][j].y - points[i][j].y))**2 + (points[i+1][j].z - points[i][j].z)**2)**0.5 
                d_u = (((points[i-1][j].x - points[i][j].x))**2 + ((points[i-1][j].y - points[i][j].y))**2 + (points[i-1][j].z - points[i][j].z)**2)**0.5
                e_r = 1 - H/d_r
                e_d = 1 - H/d_d
                e_u = 1 - H/d_u
                r_z = (points[i][j+1].z - points[i][j].z)/d_r
                r_y = (points[i][j+1].y - points[i][j].y)/(d_r*(1-r_z**2)**0.5)
                r_x = (points[i][j+1].x - points[i][j].x)/(d_r*(1-r_z**2)**0.5)
                d_z = (points[i+1][j].z - points[i][j].z)/d_d
                d_y = (points[i+1][j].y - points[i][j].y)/(d_d*(1-d_z**2)**0.5)
                d_x = (points[i+1][j].x - points[i][j].x)/(d_d*(1-d_z**2)**0.5)
                u_z = (points[i-1][j].z - points[i][j].z)/d_u
                u_y = (points[i-1][j].y - points[i][j].y)/(d_u*(1-u_z**2)**0.5)
                u_x = (points[i-1][j].x - points[i][j].x)/(d_u*(1-u_z**2)**0.5)
                if (e_r >= 0.4):
                    points_new[i][j].ra = 0
                    e_r = 0
                if (e_d >= 0.4):
                    points_new[i][j].da = 0
                    e_d = 0
                if (e_u >= 0.4):
                    points_new[i][j].ua = 0
                    e_u = 0
                points_new[i][j].dx += (e_r*r_x + e_d*d_x + e_u*u_x)*E*DELAY*DIAMETER**2*pi/(4*MASS)
                points_new[i][j].dy += (e_r*r_y + e_d*d_y + e_u*u_y)*E*DELAY*DIAMETER**2*pi/(4*MASS)
                points_new[i][j].dz += (e_r*r_z + e_d*d_z + e_u*u_z)*E*DELAY*DIAMETER**2*pi/(4*MASS)
            elif (j == POINTS_PER_SIZE - 1 and i > 0 and i < POINTS_PER_SIZE - 1):
                d_l = (((points[i][j-1].x - points[i][j].x))**2 + ((points[i][j-1].y - points[i][j].y))**2 + (points[i][j-1].z - points[i][j].z)**2)**0.5 
                d_d = (((points[i+1][j].x - points[i][j].x))**2 + ((points[i+1][j].y - points[i][j].y))**2 + (points[i+1][j].z - points[i][j].z)**2)**0.5 
                d_u = (((points[i-1][j].x - points[i][j].x))**2 + (points[i-1][j].y - points[i][j].y)**2 + (points[i-1][j].z - points[i][j].z)**2)**0.5
                e_l = 1 - H/d_l
                e_d = 1 - H/d_d
                e_u = 1 - H/d_u
                l_z = (points[i][j-1].z - points[i][j].z)/d_l
                l_y = (points[i][j-1].y - points[i][j].y)/(d_l*(1-l_z**2)**0.5)
                l_x = (points[i][j-1].x - points[i][j].x)/(d_l*(1-l_z**2)**0.5) 
                d_z = (points[i+1][j].z - points[i][j].z)/d_d
                d_y = (points[i+1][j].y - points[i][j].y)/(d_d*(1-d_z**2)**0.5)
                d_x = (points[i+1][j].x - points[i][j].x)/(d_d*(1-d_z**2)**0.5)
                u_z = (points[i-1][j].z - points[i][j].z)/d_u
                u_y = (points[i-1][j].y - points[i][j].y)/(d_u*(1-u_z**2)**0.5)
                u_x = (points[i-1][j].x - points[i][j].x)/(d_u*(1-u_z**2)**0.5)
                if (e_l >= 0.4):
                    points_new[i][j].la = 0
                    e_l = 0
                if (e_d >= 0.4):
                    points_new[i][j].da = 0
                    e_d = 0
                if (e_u >= 0.4):
                    points_new[i][j].ua = 0
                    e_u = 0
                points_new[i][j].dx += (e_l*l_x + e_d*d_x + e_u*u_x)*E*DELAY*DIAMETER**2*pi/(4*MASS) 
                points_new[i][j].dy += (e_l*l_y + e_d*d_y + e_u*u_y)*E*DELAY*DIAMETER**2*pi/(4*MASS)
                points_new[i][j].dz += (e_l*l_z + e_d*d_z + e_u*u_z)*E*DELAY*DIAMETER**2*pi/(4*MASS) 
    for i in range(1, POINTS_PER_SIZE):
        for j in range(1, POINTS_PER_SIZE):
            m1 = points[i][j].matrix((points[i][j].cy**2 + points[i][j].cz**2)**0.5, 200, 100, 100)
            m2 = points[i-1][j].matrix((points[i-1][j].cy**2 + points[i-1][j].cz**2)**0.5, 200, 100, 100)
            m3 = points[i][j-1].matrix((points[i][j-1].cy**2 + points[i][j-1].cz**2)**0.5, 200, 100, 100)
            if (points[i][j].ua == 1):
                canvas.create_line(points[i][j].x*m1[0][0] + points[i][j].y*m1[1][0] + points[i][j].z*m1[2][0] + m1[3][0], points[i][j].x*m1[0][1] + points[i][j].y*m1[1][1] + points[i][j].z*m1[2][1] + m1[3][1], points[i-1][j].x*m2[0][0] + points[i-1][j].y*m2[1][0] + points[i-1][j].z*m2[2][0] + m2[3][0], points[i-1][j].x*m2[0][1] + points[i-1][j].y*m2[1][1] + points[i-1][j].z*m2[2][1] + m2[3][1], fill = SELF_COLOR, width = 3)
            if (points[i][j].la == 1):
                canvas.create_line(points[i][j].x*m1[0][0] + points[i][j].y*m1[1][0] + points[i][j].z*m1[2][0] + m1[3][0], points[i][j].x*m1[0][1] + points[i][j].y*m1[1][1] + points[i][j].z*m1[2][1] + m1[3][1], points[i][j-1].x*m3[0][0] + points[i][j-1].y*m3[1][0] + points[i][j-1].z*m3[2][0] + m3[3][0], points[i][j-1].x*m3[0][1] + points[i][j-1].y*m3[1][1] + points[i][j-1].z*m3[2][1] + m3[3][1], fill = SELF_COLOR, width = 3)

    for j in range(1, POINTS_PER_SIZE):
        m1 = points[0][j].matrix((points[0][j].cy**2 + points[0][j].cz**2)**0.5, 200, 100, 100)
        m2 = points[0][j-1].matrix((points[0][j-1].cy**2 + points[0][j-1].cz**2)**0.5, 200, 100, 100)
        if (points[0][j].la == 1):
            canvas.create_line(points[0][j].x*m1[0][0] + points[0][j].y*m1[1][0] + points[0][j].z*m1[2][0] + m1[3][0], points[0][j].x*m1[0][1] + points[0][j].y*m1[1][1] + points[0][j].z*m1[2][1] + m1[3][1], points[0][j-1].x*m2[0][0] + points[0][j-1].y*m2[1][0] + points[0][j-1].z*m2[2][0] + m2[3][0], points[0][j-1].x*m2[0][1] + points[0][j-1].y*m2[1][1] + points[0][j-1].z*m2[2][1] + m2[3][1], fill = SELF_COLOR, width = 3)

    for i in range(1, POINTS_PER_SIZE):
        m1 = points[i][0].matrix((points[i][0].cy**2 + points[i][0].cz**2)**0.5, 200, 100, 100)
        m2 = points[i-1][0].matrix((points[i-1][0].cy**2 + points[i-1][0].cz**2)**0.5, 200, 100, 100)
        if (points[i][0].ua == 1):
            canvas.create_line(points[i][0].x*m1[0][0] + points[i][0].y*m1[1][0] + points[i][0].z*m1[2][0] + m1[3][0], points[i][0].x*m1[0][1] + points[i][0].y*m1[1][1] + points[i][0].z*m1[2][1] + m1[3][1], points[i-1][0].x*m2[0][0] + points[i-1][0].y*m2[1][0] + points[i-1][0].z*m2[2][0] + m2[3][0], points[i-1][0].x*m2[0][1] + points[i-1][0].y*m2[1][1] + points[i-1][0].z*m2[2][1] + m2[3][1], fill = SELF_COLOR, width = 3)

    for i in range(0, POINTS_PER_SIZE):
        for j in range(0, POINTS_PER_SIZE):
            points[i][j].dx = points_new[i][j].dx
            points[i][j].dy = points_new[i][j].dy
            points[i][j].dz = points_new[i][j].dz
            points[i][j].la = points_new[i][j].la
            points[i][j].ra = points_new[i][j].ra
            points[i][j].da = points_new[i][j].da
            points[i][j].ua = points_new[i][j].ua      
    VIS_DELAY = DELAY*10**9
    canvas.bind('<KeyPress>', motion)
    root.after(int(VIS_DELAY), main)


root = tkinter.Tk()
canvas = tkinter.Canvas(root, width=WIDTH, height=HEIGHT, bg=BG_COLOR)
canvas.pack()

points = [[0 for i in range(0, POINTS_PER_SIZE)] for j in range(0, POINTS_PER_SIZE)]
points_new = [[0 for i in range(0, POINTS_PER_SIZE)] for j in range(0, POINTS_PER_SIZE)]
D = [[0 for i in range(0, POINTS_PER_SIZE)] for j in range(0, POINTS_PER_SIZE)]
D_S = []

points[0][0] = Point(X_0, Y_0, Z_0, SELF_COLOR, SPEED_X, SPEED_Y, SPEED_Z, CX, CY, CZ, ANGLE, RA, 0, DA, 0)
points_new[0][0] = Point(X_0, Y_0, Z_0, SELF_COLOR, SPEED_X, SPEED_Y, SPEED_Z, CX, CY, CZ, ANGLE, RA, 0, DA, 0)
D[0][0] = 1

points[POINTS_PER_SIZE-1][POINTS_PER_SIZE-1] = Point(X_0 + DIST*(POINTS_PER_SIZE-1), Y_0 + DIST*POINTS_PER_SIZE, Z_0, SELF_COLOR, SPEED_X, SPEED_Y, SPEED_Z, CX, CY, CZ, ANGLE, 0, LA, 0, UA)
points_new[POINTS_PER_SIZE-1][POINTS_PER_SIZE-1] = Point(X_0 + DIST*(POINTS_PER_SIZE-1), Y_0 + DIST*POINTS_PER_SIZE, Z_0, SELF_COLOR, SPEED_X, SPEED_Y, SPEED_Z, CX, CY, CZ, ANGLE, 0, LA, 0, UA)
D[POINTS_PER_SIZE-1][POINTS_PER_SIZE-1] = 1

points[0][POINTS_PER_SIZE-1] = Point(X_0, Y_0 + DIST*(POINTS_PER_SIZE-1), Z_0, SELF_COLOR, SPEED_X, SPEED_Y, SPEED_Z, CX, CY, CZ, ANGLE, 0, LA, DA, 0)
points_new[0][POINTS_PER_SIZE-1] = Point(X_0, Y_0 + DIST*(POINTS_PER_SIZE-1), Z_0, SELF_COLOR, SPEED_X, SPEED_Y, SPEED_Z, CX, CY, CZ, ANGLE, 0, LA, DA, 0)
D[0][POINTS_PER_SIZE-1] = 1

points[POINTS_PER_SIZE-1][0] = Point(X_0 + DIST*(POINTS_PER_SIZE-1), Y_0, Z_0, SELF_COLOR, SPEED_X, SPEED_Y, SPEED_Z, CX, CY, CZ, ANGLE, RA, 0, 0, UA)
points_new[POINTS_PER_SIZE-1][0] = Point(X_0 + DIST*(POINTS_PER_SIZE-1), Y_0, Z_0, SELF_COLOR, SPEED_X, SPEED_Y, SPEED_Z, CX, CY, CZ, ANGLE, RA, 0, 0, UA)
D[POINTS_PER_SIZE-1][0] = 1

for i in range(1, POINTS_PER_SIZE-1):
    points[0][i] = Point(X_0, Y_0 + DIST*i, Z_0, SELF_COLOR, SPEED_X, SPEED_Y, SPEED_Z, CX, CY, CZ, ANGLE, RA, LA, DA, 0)
    points_new[0][i] = Point(X_0, Y_0 + DIST*i, Z_0, SELF_COLOR, SPEED_X, SPEED_Y, SPEED_Z, CX, CY, CZ, ANGLE, RA, LA, DA, 0)
    D[0][i] = 1

for i in range(1, POINTS_PER_SIZE-1):
    points[POINTS_PER_SIZE-1][i] = Point(X_0 + DIST*(POINTS_PER_SIZE-1), Y_0 + DIST*i, Z_0, SELF_COLOR, SPEED_X, SPEED_Y, SPEED_Z, CX, CY, CZ, ANGLE, RA, LA, 0, UA)
    points_new[POINTS_PER_SIZE-1][i] = Point(X_0 + DIST*(POINTS_PER_SIZE-1), Y_0 + DIST*i, Z_0, SELF_COLOR, SPEED_X, SPEED_Y, SPEED_Z, CX, CY, CZ, ANGLE, RA, LA, 0, UA)
    D[POINTS_PER_SIZE-1][i] = 1

for i in range(1, POINTS_PER_SIZE-1):
    points[i][0] = Point(X_0 + DIST*i, Y_0, Z_0, SELF_COLOR, SPEED_X, SPEED_Y, SPEED_Z, CX, CY, CZ, ANGLE, RA, 0, DA, UA)
    points_new[i][0] = Point(X_0 + DIST*i, Y_0, Z_0, SELF_COLOR, SPEED_X, SPEED_Y, SPEED_Z, CX, CY, CZ, ANGLE, RA, 0, DA, UA)
    D[i][0] = 1

for i in range(1, POINTS_PER_SIZE-1):
    points[i][POINTS_PER_SIZE-1] = Point(X_0 + DIST*i, Y_0 + DIST*(POINTS_PER_SIZE-1), Z_0, SELF_COLOR, SPEED_X, SPEED_Y, SPEED_Z, CX, CY, CZ, ANGLE, 0, LA, DA, UA)
    points_new[i][POINTS_PER_SIZE-1] = Point(X_0 + DIST*i, Y_0 + DIST*(POINTS_PER_SIZE-1), Z_0, SELF_COLOR, SPEED_X, SPEED_Y, SPEED_Z, CX, CY, CZ, ANGLE, 0, LA, DA, UA)
    D[i][POINTS_PER_SIZE-1] = 1

for i in range(0, POINTS_PER_SIZE):
    D_S.append(1)

for i in range(1, POINTS_PER_SIZE-1):
    for j in range(1, POINTS_PER_SIZE-1):
        points[i][j] = Point(X_0 + DIST*i, Y_0 + DIST*j, Z_0, SELF_COLOR, SPEED_X, SPEED_Y, SPEED_Z, CX, CY, CZ, ANGLE, RA, LA, DA, UA)
        points_new[i][j] = Point(X_0 + DIST*i, Y_0 + DIST*j, Z_0, SELF_COLOR, SPEED_X, SPEED_Y, SPEED_Z, CX, CY, CZ, ANGLE, RA, LA, DA, UA)
        D[i][j] = 1

canvas.focus_set()
main()
root.mainloop()