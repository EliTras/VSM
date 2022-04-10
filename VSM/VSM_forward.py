#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VSM - Version 1.0

   Volcanic and Seismic source Modelling

   Author:  Elisa Trasatti, elisa.trasatti@ingv.it
   Istituto Nazionale di Geofisica e Vulcanologia - Rome (Italy)

   Last update:  March 2022

   License:  E. Trasatti, covered by GNU-GPL License

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

      * Redistributions of source code must retain the above copyright 
        notice, this list of conditions and the following disclaimer.
      * Redistributions in binary form must reproduce the above copyright 
 	    notice, this list of conditions and the following disclaimer in 
 	    the documentation and/or other materials provided with the distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPY RIGHT OWNER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.

*******************************************************************************

   BIBLIOGRAPHY

   Davis, P.M. (1986). Surface deformation due to inflation of an arbitrarily 
       oriented triaxial ellipsoidal cavity in an elastic half-space, with
       reference to Kilauea volcano, Hawaii, J. Geophys. Res. 91, 7429-7438.
   Fialko, Y. Khazan, Y., and Simons, M. (2001). Deformation due to a  
       pressurized horizontal circular crack in an elastic half-space, with 
       applications to volcano geodesy. Geophys. J. Int. 146(1), 181-190.
   McTigue, D.F. (1987). Elastic stress and deformation near a finite spherical 
       magma body: resolution of the point-source paradox. J. Geophys. Res. 92,
        12931-12940.
   Mogi, K. (1958). Relation between the eruptions of various volcanoes and the
       deformations of the ground surfaces around them. Bull. Earth Res. Inst. 
       36, 99-134.
   Okada, Y. (1985). Surface deformation due to shear and tensile faults in a 
       half-space. B. Seismol. Soc. Am. 75, 1135-1154. 
   Yang, X., Davis, P.M., and Dieterich, J.H. (1988). Deformation from inflation
       of a dippint finite prolate spheroid in an elastic half-space as a model
       for volcanic stressing. J. Geophys. Res. 93, 4249-4257.

*******************************************************************************

    Alphabetical order

    # DAVIS   by Davis (1986)
    # MCTIGUE by McTigue (1986)
    # MOGI    by Mogi (1958)
    # OKADA   by Okada (1985)
    # PENNY   by Fialko et al. (2001)
    # YANG    by Yang et al. (1988)

"""

import numpy as np
import copy

# =====================
# Utilities
# =====================
def cart2pol(x1,x2):
    theta = np.arctan2(x2,x1)
    r = np.hypot(x2,x1)
    return theta, r

def pol2cart(theta,r):
    x1 = r * np.cos(theta)
    x2 = r * np.sin(theta)
    return x1,x2


# =====================
# DAVIS
# =====================

def davis(x,y,xcen,ycen,depth,Pxx,Pyy,Pzz,Pxy,Pyz,Pzx,nu):
    
    fact = 1./2./np.pi

    x1 = x - xcen
    x2 = y - ycen
    d = depth
    	
    R = np.sqrt(x1*x1 + x2*x2 + d*d)
    R3 = R**3
    R5 = R**5
    alfa = (3.*R + d)/R3/((R + d)**3)
    beta  = (2.*R + d)/R3/((R + d)**2)
    eta = 1./R/((R + d)**2)
    psi = 1./R/(R + d)
	      
# diagonal components
    x11a = x1*x1*alfa
    x11b = x1*x1*beta
    x22a = x2*x2*alfa
    x22b = x2*x2*beta
    anu = 1. - 2*nu
    
    U11_x = x1/2.*(-1./R3 + 3.*x1*x1/R5 + anu*(3.*eta - x11a))
    U11_y = x2/2.*(-1./R3 + 3.*x1*x1/R5 + anu*(   eta - x11a))
    U11_z = -0.5* (d/R3 - 3.*x1*x1*d/R5 - anu*(   psi - x11b))
     
    U22_x = x1/2.*(-1./R3 + 3.*x2*x2/R5 + anu*(   eta - x22a))
    U22_y = x2/2.*(-1./R3 + 3.*x2*x2/R5 + anu*(3.*eta - x22a))
    U22_z = -0.5* (d/R3 - 3.*x2*x2*d/R5 - anu*(   psi - x22b))
     
    U33_x = x1/2.*( 3.*d*d/R5 - 2.*nu/R3)
    U33_y = x2/2.*( 3.*d*d/R5 - 2.*nu/R3)
    U33_z = -d/2.*(-3.*d*d/R5 + 2.*nu/R3)
	      	      
# off-diagonal components      
    W12_x = x2*(3.*x1*x1/R5 + anu*(eta - x11a))
    W12_y = x1*(3.*x2*x2/R5 + anu*(eta - x22a))
    W12_z = -x1*x2*(-3*d/R5 + anu*beta)
	
    W13_x = -x1*x1*(-3.*d/R5)
    W13_y = -x1*x2*(-3.*d/R5)
    W13_z =  x1*d*(  3.*d/R5)
	
    W23_x = W13_y
    W23_y = -x2*x2*(-3.*d/R5)
    W23_z =  x2*d*(  3.*d/R5)
    
# surface displacements computed as a linear combinantion
    Ue = fact*(Pxx*U11_x + Pyy*U22_x + Pzz*U33_x + 
               Pxy*W12_x + Pyz*W23_x + Pzx*W13_x)
    Un = fact*(Pxx*U11_y + Pyy*U22_y + Pzz*U33_y + 
               Pxy*W12_y + Pyz*W23_y + Pzx*W13_y)
    Uz = fact*(Pxx*U11_z + Pyy*U22_z + Pzz*U33_z + 
               Pxy*W12_z + Pyz*W23_z + Pzx*W13_z)  
    
    return Ue, Un, Uz


# =====================
# McTIGUE
# =====================
def mctigue(x,y,xcen,ycen,depth,radius,dP_mu,nu):

    # Center coordinate grid on point source
    x = x - xcen
    y = y - ycen

    # Convert to surface cylindrical coordinates
    th, rho = cart2pol(x,y)
    r = rho/depth
    
    # McTigue displacement calculation
    R3 = (r**2 + 1.)**1.5
    R5 = (r**2 + 1.)**2.5
    a_d = radius/depth
    a_d3 = a_d**3
    fact1 = (1.-nu) * (1.+nu)/(2.*(7. - 5.*nu))
    fact2 = 15.*(2.-nu) * (1.-nu)/(4.*(7. - 5.*nu))
    C = dP_mu * depth * a_d3
    
    uz = C * (1.-nu)/R3 - C * a_d3*(fact1/R3) - C*a_d3*fact2/R5
    ur = uz*r

    ux, uy = pol2cart(th, ur)

    return np.array([ux,uy,uz])



# =====================
# MOGI
# =====================
def mogi(x,y,xcen,ycen,depth,dVol,nu):

    # Center coordinate grid on point source
    x = x - xcen
    y = y - ycen

    # Convert to surface cylindrical coordinates
    th, rho = cart2pol(x,y)
    R = np.hypot(depth,rho)

    # Mogi displacement calculation
    C = dVol * ((1.-nu) / np.pi)
    ur = C * rho / R**3
    uz = C * depth / R**3

    ux, uy = pol2cart(th, ur)

    return np.array([ux,uy,uz])


# =====================
# OKADA
# =====================

eps = 1e-14 #numerical constant

def okada(x, y, xtlc, ytlc, dtlc, length, width,
            strike, dip, param1,param2, opening, opt,nu):
    
    L = length
    W = width

    if(abs(strike) < 0.0001):
        if(strike == 0.):
            strike = 0.0001
        else:
            strike = np.sign(strike)*0.0001

    strike = np.deg2rad(strike) #transformations accounted for below
    dip = np.deg2rad(dip)
    
    #print('OPT',opt)
    
    if(opt == 'S'):
        U1 = param1
        U2 = param2
    else:
        slip = param1
        rake = np.deg2rad(param2)
        U1 = np.cos(rake) * slip
        U2 = np.sin(rake) * slip
    U3 = opening

    xbrc = xtlc + W*np.cos(dip)*np.cos(strike) + L*np.sin(strike)
    ybrc = ytlc - W*np.cos(dip)*np.sin(strike) + L*np.cos(strike)
    
    e = x - xbrc
    n = y - ybrc

    d = dtlc + np.sin(dip) * W #fault bottom edge

    ec = e + np.cos(strike) * np.cos(dip) * W
    nc = n - np.sin(strike) * np.cos(dip) * W
    x = np.cos(strike) * nc + np.sin(strike) * ec + L
    y = np.sin(strike) * nc - np.cos(strike) * ec + np.cos(dip) * W
    
    p = y * np.cos(dip) + d * np.sin(dip)
    q = y * np.sin(dip) - d * np.cos(dip)

    ux = - U1 / (2 * np.pi) * chinnery(ux_ss, x, p, L, W, q, dip, nu) - \
           U2 / (2 * np.pi) * chinnery(ux_ds, x, p, L, W, q, dip, nu) + \
           U3 / (2 * np.pi) * chinnery(ux_tf, x, p, L, W, q, dip, nu)

    uy = - U1 / (2 * np.pi) * chinnery(uy_ss, x, p, L, W, q, dip, nu) - \
           U2 / (2 * np.pi) * chinnery(uy_ds, x, p, L, W, q, dip, nu) + \
           U3 / (2 * np.pi) * chinnery(uy_tf, x, p, L, W, q, dip, nu)

    uz = - U1 / (2 * np.pi) * chinnery(uz_ss, x, p, L, W, q, dip, nu) - \
           U2 / (2 * np.pi) * chinnery(uz_ds, x, p, L, W, q, dip, nu) + \
           U3 / (2 * np.pi) * chinnery(uz_tf, x, p, L, W, q, dip, nu)

    ue = np.sin(strike) * ux - np.cos(strike) * uy
    un = np.cos(strike) * ux + np.sin(strike) * uy

    return ue,un,uz


def chinnery(f, x, p, L, W, q, dip, nu):
    u =  (f(x, p, q, dip, nu) -
          f(x, p - W, q, dip, nu) -
          f(x - L, p, q, dip, nu) +
          f(x - L, p - W, q, dip, nu))
    return u


def ux_ss(xi, eta, q, dip, nu):

    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    u = xi * q / (R * (R + eta)) + \
        I1(xi, eta, q, dip, nu, R) * np.sin(dip)
    k = (q != 0)
    u[k] = u[k] + np.arctan( (xi[k] * eta[k]) / (q[k] * R[k]) )
    return u


def uy_ss(xi, eta, q, dip, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    u = (eta * np.cos(dip) + q * np.sin(dip)) * q / (R * (R + eta)) + \
        q * np.cos(dip) / (R + eta) + \
        I2(eta, q, dip, nu, R) * np.sin(dip)
    return u


def uz_ss(xi, eta, q, dip, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    db = eta * np.sin(dip) - q * np.cos(dip)
    u = (eta * np.sin(dip) - q * np.cos(dip)) * q / (R * (R + eta)) + \
        q * np.sin(dip) / (R + eta) + \
        I4(db, eta, q, dip, nu, R) * np.sin(dip)
    return u


def ux_ds(xi, eta, q, dip, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    u = q / R - \
        I3(eta, q, dip, nu, R) * np.sin(dip) * np.cos(dip)
    return u


def uy_ds(xi, eta, q, dip, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    u = ( (eta * np.cos(dip) + q * np.sin(dip)) * q / (R * (R + xi)) -
           I1(xi, eta, q, dip, nu, R) * np.sin(dip) * np.cos(dip) )
    k = (q != 0)
    u[k] = u[k] + np.cos(dip) * np.arctan( (xi[k] * eta[k]) / (q[k] * R[k]))
    return u


def uz_ds(xi, eta, q, dip, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    db = eta * np.sin(dip) - q * np.cos(dip)
    u = ( db * q / (R * (R + xi)) -
          I5(xi, eta, q, dip, nu, R, db) * np.sin(dip) * np.cos(dip) )
    k = (q != 0)
    #u[k] = u[k] + np.sin(dip) * np.arctan2(xi[k] * eta[k] , q[k] * R[k])
    u[k] = u[k] + np.sin(dip) * np.arctan( (xi[k] * eta[k]) / (q[k] * R[k]))
    return u


def ux_tf(xi, eta, q, dip, nu):
    R = np.sqrt(xi**2 + eta**2 + q**2)
    u = q**2 / (R * (R + eta)) - \
        (I3(eta, q, dip, nu, R) * np.sin(dip)**2)
    return u


def uy_tf(xi, eta, q, dip, nu):
    R = np.sqrt(xi**2 + eta**2 + q**2)
    u = - (eta * np.sin(dip) - q * np.cos(dip)) * q / (R * (R + xi)) - \
        (np.sin(dip) * xi * q / (R * (R + eta))) - \
        (I1(xi, eta, q, dip, nu, R) * np.sin(dip) ** 2)
    k = (q != 0)
    #u[k] = u[k] + np.sin(dip) * np.arctan2(xi[k] * eta[k] , q[k] * R[k])
    u[k] = u[k] + np.sin(dip) * np.arctan( (xi[k] * eta[k]) / (q[k] * R[k]) )
    return u


def uz_tf(xi, eta, q, dip, nu):
    R = np.sqrt(xi**2 + eta**2 + q**2)
    db = eta * np.sin(dip) - q * np.cos(dip)
    u = (eta * np.cos(dip) + q * np.sin(dip)) * q / (R * (R + xi)) + \
         np.cos(dip) * xi * q / (R * (R + eta)) - \
         I5(xi, eta, q, dip, nu, R, db) * np.sin(dip)**2
    k = (q != 0)
    u[k] = u[k] - np.cos(dip) * np.arctan( (xi[k] * eta[k]) / (q[k] * R[k]) )
    return u


def I1(xi, eta, q, dip, nu, R):
    db = eta * np.sin(dip) - q * np.cos(dip)
    if np.cos(dip) > eps:
        I = (1 - 2 * nu) * (- xi / (np.cos(dip) * (R + db))) - \
            np.sin(dip) / np.cos(dip) * I5(xi, eta, q, dip, nu, R, db)
    else:
        I = -(1 - 2 * nu)/2 * xi * q / (R + db)**2
    return I


def I2(eta, q, dip, nu, R):
    I = (1 - 2 * nu) * (-np.log(R + eta)) - \
        I3(eta, q, dip, nu, R)
    return I


def I3(eta, q, dip, nu, R):
    yb = eta * np.cos(dip) + q * np.sin(dip)
    db = eta * np.sin(dip) - q * np.cos(dip)
    if np.cos(dip) > eps:
        I = (1 - 2 * nu) * (yb / (np.cos(dip) * (R + db)) - np.log(R + eta)) + \
            np.sin(dip) / np.cos(dip) * I4(db, eta, q, dip, nu, R)
    else:
        I = (1 - 2 * nu) / 2 * (eta / (R + db) + yb * q / (R + db) ** 2 - np.log(R + eta))
    return I


def I4(db, eta, q, dip, nu, R):
    if np.cos(dip) > eps:
        I = (1 - 2 * nu) * 1.0 / np.cos(dip) * \
            (np.log(R + db) - np.sin(dip) * np.log(R + eta))
    else:
        I = - (1 - 2 * nu) * q / (R + db)
    return I


def I5(xi, eta, q, dip, nu, R, db):
    X = np.sqrt(xi**2 + q**2)
    if np.cos(dip) > eps:
#        print(np.sin(dip),  min((eta * (X + q*np.cos(dip)) + X*(R + X) * np.sin(dip)) /
#                        (xi*(R + X) * np.cos(dip))),max((eta * (X + q*np.cos(dip)) + X*(R + X) * np.sin(dip)) /
#                        (xi*(R + X) * np.cos(dip))) )
        I = (1 - 2 * nu) * 2 / np.cos(dip) * \
             np.arctan( (eta * (X + q*np.cos(dip)) + X*(R + X) * np.sin(dip)) /
                        (xi*(R + X) * np.cos(dip)) )
        I[xi == 0] = 0
    else:
        I = -(1 - 2 * nu) * xi * np.sin(dip) / (R + db)
    return I


# =====================
# FIALKO
# =====================

def penny(x,y,xcen,ycen,depth, radius, dP_mu, nu):

# number of sub-intervals on [0,1] on which integration is done using a 16-point 
# Gauss quadrature (i.e., total of nis*16 points)
    nis = 2  
# solution accuracy for Fredholm integral equations (stop
# iterating if relative change is less than eps) 
    eps = 1e-6 #or 1e-8
    
    h = depth/radius

    xn = x - xcen
    yn = y - ycen
    
    xn = np.where(xn != 0.0, xn, 0.001)
    yn = np.where(yn != 0.0, yn, 0.001)
    
    fi, psi, t, Wt = FREDHOLM(h, nis, eps)

    r = np.hypot(xn,yn)
    r = (r/depth)*h
    
    U3 = np.zeros(np.size(xn))
    Ur = np.zeros(np.size(xn))

    for t_i in range(0, np.size(t)):
        fi_sub  = fi[t_i]
        psi_sub = psi[t_i]
        Wt_sub  = Wt[t_i]
        t_sub   = t[t_i]
        
        Ur_app , U3_app = INTGR(r,fi_sub,psi_sub,h,Wt_sub,t_sub)
        
        U3 += U3_app
        Ur += Ur_app
   
    theta=np.arctan(yn/xn)
    theta = np.where(xn >= 0., theta, theta + np.pi)
      
    coeffialko = 2.*(1.-nu)*dP_mu*radius
    Ue = Ur*np.cos(theta)*coeffialko   
    Un = Ur*np.sin(theta)*coeffialko
    Uz = -U3*coeffialko

    return Ue, Un, Uz


def KG(s,p):
    z = s**2
    y = p + z
    KKGG = (3*p - z)/y**3

    return KKGG

def KERN(w,p):
    u = (p + w**2)**3
    KKERN = (p**2 + 2*w**4 -p*w**2)/u
    return KKERN

def FPKERNEL(h,t,r,n):
    p = 4*h**2

    nr = np.size(r)
    y = (t + r)**2
    z = (t - r)**2
    
    a = np.zeros(nr)
    b = np.zeros(nr)
    trbl = np.ones(nr)    

    if n == 1:
        K = p*h*(KG(t-r,p)-KG(t+r,p))
    elif n == 2:
        Dlt = 1e-6
        a = t + r
        b = t - r
        g = 2*p*h*(p**2 + 6*p*(t**2 + r**2) + 5*(a*b)**2)
        s = ((p + z)*(p + y))**2
        s = g/s
        trbl *= -4*h/(p + t**2)	
        if t<Dlt:
           trbl = -4*h/(p + r**2)
        else:
           for j in range(0, nr):
              if r[j] > Dlt:
                  lpz_py = np.log((p+z[j])/(p+y[j]))
                  trbl[j] = h/t/r[j]*lpz_py
              else:
                 trbl[j] = -4*h/(p + t**2)
        KERN_bp = KERN(b,p)
        KERN_ap = KERN(a,p)
        K = trbl + s + h*(KERN_bp + KERN_ap)
    elif n == 3:
        a = ((p + y)*(p + z))**2
        c = t + r
        d = t - r
        b = p*t*((3*p**2 - (c*d)**2 + 2*p*(t**2 + r**2))/a)
        a = p/2*(c*KG(c,p)+d*KG(d,p))
        K = b - a
    elif n == 4:
        a = ((p + y)*(p + z))**2
        c = t + r
        d = -t  + r
        b = p*r*(3*p**2 - (c*d)**2 + 2*p*(t**2 + r**2))/a
        a = p/2*(c*KG(c,p) + d*KG(d,p))
        K = b - a
    return K

def FREDHOLM(h, m, er):

    lamda = 2/np.pi

    NumLegendreTerms = 16
    Root, Weight = np.polynomial.legendre.leggauss(NumLegendreTerms)

    t  = np.zeros(m*NumLegendreTerms)
    Wt = np.zeros(np.size(t))
    
    for k in range(1, m+1):
      	for i in range(0, NumLegendreTerms):
              d1 = 1./m
              t1 = d1*(k - 1.)
              r1 = d1*k
              j = NumLegendreTerms*(k - 1) + i
              t[j] = Root[i]*(r1 - t1)/2 + (r1 + t1)/2
              Wt[j] = 0.5*Weight[i]/m
    
    fi1  = -lamda*t
    psi1 = np.zeros(np.size(t))
    fi   = np.zeros(np.size(t))
    psi  = np.zeros(np.size(t))
    res = 1e9
    
    while res > er:
        for i in range(0, m*NumLegendreTerms):
            ti = t[i]
            fpkern1 = FPKERNEL(h,ti,t,1)
            fpkern2 = FPKERNEL(h,ti,t,2)
            fpkern3 = FPKERNEL(h,ti,t,3)
            fpkern4 = FPKERNEL(h,ti,t,4)
            fi[i]  = np.sum(Wt*(fi1*fpkern1 + psi1*fpkern3)) - ti
            psi[i] = np.sum(Wt*(psi1*fpkern2 + fi1*fpkern4))
        
        fi   = fi*lamda
        psi  = psi*lamda
        fim  = np.amax(np.abs(fi1-fi))
        im   = np.argmax(np.abs(fi1-fi))
        fim  = fim/np.abs(fi[im])
        psim = np.amax(np.abs(psi1 - psi))
        im   = np.argmax(np.abs(psi1-psi))
        psim = psim/np.abs(psi[im])
        res  = np.max([fim,psim])
        fi1  = copy.deepcopy(fi)
        psi1 = copy.deepcopy(psi)
        if res < er:
            break

    return  fi, psi, t, Wt

def INTGR(r, fi, psi, h, Wt, t):
#   Uz(r),Ur(r) - vertical and radial displacements
#   fi,psi: basis functions
#   t: interval of integration

    E = h**2 + r**2 - t**2
    D = np.sqrt(E**2 + 4*(h*t)**2)
    D3 = D**3

    rad2 = np.sqrt(2.)
    sdme = np.sqrt(D - E)
    sdpe = np.sqrt(D + E)
    
    Uz = Wt * (fi*(rad2*h*t/(D*sdpe) + h/rad2/D3*(h*sdme*(2*E + D) - 
              t*sdpe*(2*E - D))) + psi * (rad2*h*t/(D*sdpe)/t  - 
              1/rad2/D3*(h*sdpe*(2*E - D) + t*sdme*(2*E + D)))
              )
    
    Ur = Wt * (psi *((t/r - sdme/r/rad2 - h*(-(h*sdme - 
              t*sdpe)/D/r/rad2))/t - (1/r - (h*sdpe + t*sdme)/D/r/rad2) + 
              h*r*sdpe*(2*E - D)/D3/rad2) - h*fi*r*sdme*(2*E + D)/D3/rad2
              )

    return Ur, Uz

# =====================
# YANG 
# =====================
def yang(x,y,xcen,ycen,depth,s_axis_max,ratio,dP_mu,strike,dip,nu):
    
    A = s_axis_max
    PI = np.pi

    if(dip < 0.0001):
        dip = 0.001
    if(dip > 89.995):
        dip = 89.995
    if(strike == 0.0):
        strike = 0.001
    if(ratio > 1): ratio = 1./ratio #Implicitly works for the opposite
    
    ANG = np.deg2rad(dip) 
    strike = np.deg2rad(strike)
    cs = np.cos(strike)
    ss = np.sin(strike)
    CA=np.cos(ANG)
    

    ZZ0 = depth
    Z = 0. #free surface at Z = 0

    # Find the coordinates in the spheroid's reference system 
    xx1 = (y - ycen)*cs - (x - xcen)*ss
    xx2 = (y - ycen)*ss + (x - xcen)*cs
    
    xx1 = np.where(xx1 != 0.0, xx1, 0.001)
    xx2 = np.where(xx2 != 0.0, xx2, 0.001)
    
    B = A*ratio
    C = np.sqrt(A*A - B*B)

    AM = 1.e9   # actually not involved in the final results
    SI = nu
    AL = 2.*AM*SI/(1.-2.*SI)

    # pressure required as -ve
    P = -dP_mu*AM


    PC, PD = pf(P,AL,AM,SI,A,B,C,PI)
    
    D1=-2.*B**2*PD
    D2=3.*B**2/C**2*PD+2.*(1.-2.*SI)*PC
    V=A*B**2/C**3/(16.*AM*(1.-SI))


    # first call with +C
    FF1,FF2,FF3,W1,W2,W3 = ellip(D1,D2,PD,xx1,xx2,Z,C,ANG,ZZ0,SI)

    UFF1=FF1
    UFF2=FF2
    UFF3=FF3
    U1=W1
    U2=W2
    U3=W3

    # second call with -C
    FF1,FF2,FF3,W1,W2,W3 = ellip(D1,D2,PD,xx1,xx2,Z,-C,ANG,ZZ0,SI)
    
    UFF1 = (UFF1-FF1)/CA**2
    UFF2 = (UFF2-FF2)/CA**2
    UFF3 = (UFF3-FF3)/CA
    U1 = V*(U1-W1+UFF1)
    U2 = V*(U2-W2+UFF2)
    U3 = V*(U3-W3+UFF3)

    Ue =  -U1*ss + U2*cs
    Un =   U1*cs + U2*ss
    Uz =  -U3

    return Ue, Un, Uz

def pf(P,AL,AM,SI,A,B,C,PI):
# Function to compute source parameters due to the overpressure and the
# source boundary
    
    AC = (A - C) / (A + C)
    L1 = np.log(AC)
    AI=-2.0*PI*A*B**2*(2.0/A/C**2 +L1/C**3) 
    AAI=-2.0*PI*A*B**2*(2.0/3.0/C**2/A**3+2.0/A/C**4+L1/C**5)
    
    U=8.0*PI*(1.-SI)
    Q=3.0/U
    R=(1.-2.*SI)/U
    
    A11=2.*R*(AI-4.*PI)
    A12=-2.*R*(AI+4.*PI)
    A21=Q*A**2*AAI+R*AI-1
    A22=-(Q*A**2*AAI+(2.*R-Q)*AI)
    
    W=(A11*A22-A12*A21)*(3.*AL+2.*AM)
    E11=(3.*A22-A12)*P/W
    E22=(A11-3.*A21)*P/W
    
    PD=2.*AM*(E11-E22)
    PC=AL*E11+2.*(AL+AM)*E22
    
    return PC,PD


def ellip(D1,D2,PD,X,Y,Z,C,ANG,ZZ0,SI):
# Function to compute the primitive of the displacement field 
# at the benchmark. It is called twice (for +C and -C) 
      
    eps = 1e-15
    SA = np.sin(ANG)
    CA = np.cos(ANG)
    if np.abs(SA)< eps:
       SA = np.sign(SA)*eps
   
    C2 = CA*C
    C3 = SA*C
    
    Y3 = Z - ZZ0
    T2 = SA*Y - CA*Y3
    T3 = CA*Y + SA*Y3
    TT3 = T3 - C
    YY3 = Z + ZZ0
    Q2 =  SA*Y + CA*YY3
    Q3 = -CA*Y + SA*YY3
    QQ3 = Q3 + C

    X1 = X
    X2 = Y - C2
    X3 = Y3 - C3
    XX3 = YY3 + C3
    C0 = ZZ0/SA
    
    R1=np.sqrt(X1*X1 + X2*X2 + X3*X3)
    R2=np.sqrt(X1*X1 + X2*X2 + XX3*XX3)
    
    Rr = R1 + TT3
    Rq = R2 + QQ3
    Ry = R2 + XX3
    lRr = np.log(Rr)
    lRq = np.log(Rq)
    lRy = np.log(Ry)

    BT1=CA*Q2+(1.+SA)*Rq
    BT2=CA*X1 + eps

    atnbeta = np.arctan(BT1 / BT2)
    
    
    A1 = C/R1 + lRr
    A2 = R1 - T3*lRr
    A3 = C*TT3/R1 + R1
    
    AA1 = C/R2- lRq
    AA2 = R2-Q3*lRq
    AA3 = C*QQ3/R2-R2
    
    B = C*(C+C0)/R2-AA2 - C0*lRq
    
    AC  =  D1/R1/Rr + D2*(lRr + (T3+C)/Rr)
    AAC = -D1/R2/Rq - D2*(lRq + (Q3-C)/Rq)
    BC =  (D1/R1+2.*D2*A2)+(3.-4.*SI)*(D1/R2+2.*D2*AA2)

      
    F1 = -2.*SA*Z*(C*(C+C0)/R2**3+(R2+C+C0)/R2/Rq + 
                   4.*(1.-SI)*(R2+C)/R2/Rq
                   )
    F2 = -2.*SA*Z*(C*(C+C0)*QQ3/R2**3+C0/R2 + 
                   AA1 + 4.*(1.-SI)*AA1
                   )
    
    FC1=2.*Z*(CA*Q2*(D1*2.*Rq/R2**3/Rq**2-D2*(R2+2.*C)/R2/Rq**2) +
              SA*(D1/R2**3 - 
              2.*D2*(R2+C)/R2/Rq)
              )
    FC2 = 2.*Z*(D1*XX3/R2**3 - 
                2.*D2*(SA*AA1+CA*Q2*(R2+C)/R2/Rq)
                )
    
    FD1 = (CA**2*C*X1/Ry + 3.*(SA*X1*lRy - X1*lRq + 2.*Q2*atnbeta) + 
          CA**2*2.*X1*lRq - 4.*YY3*atnbeta*CA
          )
          
    FD2 = (CA**2*C*X2/Ry + 3.*(SA*Q2*lRq - Q2*lRy + 2.*SA*X1*atnbeta + 
           CA*(R2-XX3)) - 2.*CA**3*AA2 + 2.*CA*(YY3*lRy - Q3*lRq)
          )
      
    FD3 = (Q2*lRq - Q2*SA*lRy + 2.*X1*atnbeta) + (2.*SA*AA2+Q3*lRy - C)*CA
    
    UC1=(AC+(3.-4.*SI)*AAC+FC1)*X1
    UC2=SA*(AC*T2+(3.-4.*SI)*AAC*Q2+FC1*Q2)+CA*(BC-FC2)+2.*SA*CA*Z*AAC
    UC3=-CA*(AC*T2+(3.-4.*SI)*Q2*AAC-FC1*Q2)+SA*(BC+FC2)+2.*CA**2*Z*AAC
    
    UD1=X1*A1+X1*((3.-4.*SI)*AA1+F1)
    UD2=SA*( T2*A1+Q2*((3.-4.*SI)*AA1+F1))+4.*(1.-SI)*CA*(A2+AA2)+CA*(A3-(3.-4.*SI)*AA3-F2)
    UD3=CA*(-T2*A1+Q2*((3.-4.*SI)*AA1+F1))+4.*(1.-SI)*SA*(A2+AA2)+SA*(A3+(3.-4.*SI)*AA3+F2-2.*(3.-4.*SI)*B)
      
    U1 = UC1 + 2.*PD*UD1
    U2 = UC2 + 2.*PD*UD2
    U3 = UC3 + 2.*PD*UD3
     
    FF1 = -2.*PD*4.*(1.-SI)*(1.-2.*SI)*FD1
    FF2 = -2.*PD*4.*(1.-SI)*(1.-2.*SI)*FD2
    FF3 =  2.*PD*4.*(1.-SI)*(1.-2.*SI)*FD3

    return FF1, FF2, FF3, U1, U2, U3