#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Friday May29 2020

@email ---  aamiri@pks.mpg.de
"""

import numpy as np
import scipy.linalg as lin
from math import floor
import json
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.collections import LineCollection
from matplotlib import colors as mcolors


# ===========  Eleongation tensor ============
def rotate_triangle_into_xy_plane(r21, r31):
    epsValue=1e-20
    n = np.cross(r21, r31)
    n_norm = lin.norm(n) + epsValue
    n = n/n_norm
    triArea = 0.5 * n_norm
    theta_x, theta_y = np.arctan2(n[1], n[2]), np.arctan2(n[0], np.sqrt(1.0-n[0]*n[0]))
    cosX, cosY = np.cos(theta_x), np.cos(-theta_y)
    sinX, sinY = np.cos(theta_x), np.sin(-theta_y)

    ryrx_11, ryrx_12, ryrx_13 = cosY, sinX*sinY, sinY*cosX
    ryrx_21, ryrx_22, ryrx_23 = 0, cosX, -sinX
    ryrx_31, ryrx_32, ryrx_33 = -sinY, cosY*sinX, cosY*cosX

    new_r21x = ryrx_11*r21[0] + ryrx_12*r21[1] + ryrx_13*r21[2]
    new_r21y = ryrx_21*r21[0] + ryrx_22*r21[1] + ryrx_23*r21[2]
    new_r21z = ryrx_31*r21[0] + ryrx_32*r21[1] + ryrx_33*r21[2]

    new_r31x = ryrx_11*r31[0] + ryrx_12*r31[1] + ryrx_13*r31[2]
    new_r31y = ryrx_21*r31[0] + ryrx_22*r31[1] + ryrx_23*r31[2]
    new_r31z = ryrx_31*r31[0] + ryrx_32*r31[1] + ryrx_33*r31[2]

    return np.array([new_r21x,new_r21y,new_r21z]), np.array([new_r31x,new_r31y,new_r31z]), triArea


def get_triangle_elongation(P1, P2, P3):
    eps_val = 1e-20
    a0 = np.sqrt(3.0)/4.0;
    c0 = 2.0 * np.sqrt(a0) * 3.0**0.25
    ha, hb, sc, sd = 0, 0, 0, 0 #scaled rotation (ha, hb), and traceless symmetric (sc, sd)
    C21, C31 = np.array([c0,0]),np.array([0.5*c0,0.5*np.sqrt(3)*c0]) # C is the reference triange and R is the deformed
    detC = C21[0] * C31[1] - C21[1] * C31[0]

    R21, R31 = P2-P1, P3-P1
    A_vec = 0.5 * np.cross(R21,R31)
    triArea = lin.norm(A_vec)
    #R21, R31, triArea = rotate_triangle_into_xy_plane(r21, r31)

    if triArea<eps_val:
        return 0, 0, 0, triArea

    Sxx = (C31[1] * R21[0] - C21[1] * R31[0])/detC #deformation tensor of a triangle
    Sxy = (-C31[0] * R21[0] + C21[0] * R31[0])/detC
    Syx = (C31[1] * R21[1] - C21[1] * R31[1])/detC
    Syy = (C21[0] * R31[1] - C31[0] * R21[1])/detC

    ha = 0.5 * (Sxx + Syy)
    hb = 0.5 * (Syx - Sxy)
    sc = 0.5 * (Sxx - Syy)
    sd = 0.5 * (Sxy + Syx)
    hNorm = np.sqrt(ha*ha + hb * hb)
    S_tilda_norm = np.sqrt(sc*sc + sd*sd) + eps_val

    sTildaRnegTheta_xx = (ha * sc - hb * sd)/hNorm
    sTildaRnegTheta_xy = (hb * sc + ha * sd)/hNorm
    sTildaRnegTheta_yy = (-ha * sc + hb * sd)/hNorm

    QpreFactor = (1.0/S_tilda_norm) * np.arcsinh( np.sqrt(a0/triArea) * S_tilda_norm)
    Qtens_xx = QpreFactor * sTildaRnegTheta_xx
    Qtens_xy = QpreFactor * sTildaRnegTheta_xy
    Qtens_yy = QpreFactor * sTildaRnegTheta_yy

    return Qtens_xx, Qtens_xy, Qtens_yy, triArea



def get_cell_elongation(vertices):
    cell_cent = np.mean(vertices, axis=0)
    cell_Qxx, cell_Qxy, cell_Qyy = 0, 0, 0
    normQtens = 0
    cell_area = 0

    nv = len(vertices)
    for i in range(nv):
        next_id = i+1
        if i==nv-1:
            next_id = 0
        Qtens_xx, Qtens_xy, Qtens_yy, triArea = get_triangle_elongation(cell_cent, vertices[i], vertices[next_id])
        cell_area = cell_area + triArea
        normQtens = normQtens + triArea*np.sqrt(Qtens_xx*Qtens_xx + Qtens_xy*Qtens_xy)

        cell_Qxx = cell_Qxx + triArea*Qtens_xx
        cell_Qxy = cell_Qxy + triArea*Qtens_xy
        cell_Qyy = cell_Qyy + triArea*Qtens_yy

    cell_Qxx = cell_Qxx/cell_area
    cell_Qxy = cell_Qxy/cell_area
    cell_Qyy = cell_Qyy/cell_area
    Q_tens = np.array([[cell_Qxx, cell_Qxy], [cell_Qxy, cell_Qyy]])
    eigvals, eigvecs = lin.eig(Q_tens)
    max_index = np.argmax(eigvals)
    return eigvals[max_index], eigvecs[:,max_index], cell_cent


#=============================================
def is_point_near_bdr(px, py, system_size_x, system_size_y):
    if px<0.1*system_size_x or px>0.9*system_size_x:
        return True

    if py<0.1*system_size_y or py>0.9*system_size_y:
        return True

    return False


def get_cell_vertices(c_edges, tissue_edges, tissue_vertices, system_size_x, system_size_y):
    print(c_edges)
    num_sides = len(c_edges)
    #add_point_if_not_there(list_points, p)
    #for p in
    #print(num_sides)
    c_verts = np.zeros((num_sides,2))
    c_verts[0,:] = np.array(c_edges[0])
    v_this = np.array(c_edges[0])
    for eid in np.arange(1,num_sides):
        v_this = v_this - np.array(c_edges[eid])
        c_verts[eid,:]=v_this
    return c_verts

def plot_cell_edges(tail_edges, dl_edges, image_name):
    min_x, min_y = 1000, 1000
    max_x, max_y = -1000, -1000
    for tail in tail_edges:
        if tail[0]<min_x:
            min_x = tail[0]
        elif tail[0]>max_x:
            max_x = tail[0]

        if tail[1]<min_y:
            min_y = tail[1]
        elif tail[1]>max_y:
            max_y = tail[1]

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1, 1, 1)
    for tail, dl in zip(tail_edges, dl_edges):
        ax.arrow(tail[0], tail[1], dl[0], dl[1], width=0.02, length_includes_head=True)
    plt.xlim([min_x-1, max_x+1])
    plt.ylim([min_y-1, max_y+1])
    plt.savefig(image_name, dpi=300)
    plt.close()

def visualize_frame_jsonfile(datafile, image_name, color_type):
    if color_type=='d2W_dv2' or color_type=='dW_dh':
        fig = plt.figure(figsize=(17,9.5))
    else:
        fig = plt.figure(figsize=(14,14))
    ax = fig.add_subplot(1, 1, 1)
    patches = []
    with open(datafile, 'r') as fp:
        data = json.load(fp)

    ver_list = data["vertices"]
    tissue_edges = data["edges"]

    #plotting line tension
    line_coords = []
    edge_colors = []
    for te in tissue_edges:
        crosses, t_id, h_id, t_qx, t_qy, h_qx, h_qy, tension = te[0], te[1], te[2], te[3], te[4], te[5], te[6], te[7]
        cross2 = abs(t_qx)+abs(t_qy)+abs(h_qx)+abs(h_qy)
        if crosses or cross2:
            continue
        else:
            line_coords.append([(ver_list[t_id][0], ver_list[t_id][1]), (ver_list[h_id][0], ver_list[h_id][1])])
            edge_colors.append(tension)

    lc = LineCollection(line_coords, array=np.array(edge_colors), cmap=plt.cm.rainbow, alpha=0.5)#, linewidths=2)
    ax.add_collection(lc)
    lcc = plt.colorbar(lc)
    lcc.ax.set_title('$\gamma$')


    num_polygons = len(data["cells"])#data["n_cells"]

    nematic_coords = []
    nematic_colors = []
    for cid in range(num_polygons):
        c_edges = data["cells"][cid]
        #c_verts =  get_cell_vertices(c_edges, data["edges"], data["vertices"], data["box_size"][0], data["box_size"][1])
        Q_val, Q_vec, cell_cent = get_cell_elongation(c_edges)
        Q_val = np.real(Q_val)
        Q_vec = np.real(Q_vec)
        Q_tail, Q_head = cell_cent-3.0*Q_val*Q_vec, cell_cent+3.0*Q_val*Q_vec
        nematic_coords.append([(Q_tail[0], Q_tail[1]), (Q_head[0], Q_head[1])])
        #nematic_colors.append(np.real(Q_val))
    #print(nematic_coords)
    lqc = LineCollection(nematic_coords, cmap=plt.cm.rainbow, alpha=0.8, linewidths=2)
    lqc.set_color("black")
    ax.add_collection(lqc)

    poly_class = []
    for cid in range(num_polygons):
        c_edges = data["cells"][cid]
        num_sides = len(c_edges)
        poly_class.append(num_sides)
        #c_verts =  get_cell_vertices(c_edges, data["edges"], data["vertices"], data["box_size"][0], data["box_size"][1])
        polygon = Polygon(np.array(c_edges))#, True)
        patches.append(polygon)
    p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.05)
    if color_type=='area_color':
        colors = data["cell_area"]
        p.set_array(np.array(colors))
    else:
        colors = np.array(poly_class)
        p.set_array(np.array(colors))
        #colors = np.zeros(len(data["cell_area"]))
    p.set_ec('k')
    ax.add_collection(p)


    #now add the vertex bending rigidity estimates if  asked
    if color_type=='d2W_dv2':
        colors = data["d2W_dv2"]
        x = [v[0] for v in ver_list]
        y = [v[1] for v in ver_list]
        v_area = 40*np.ones(len(colors))
        sc = ax.scatter(x, y, s=v_area, c=colors)#, alpha=0.0)
        clb = plt.colorbar(sc)
        clb.ax.set_title('$d^2W/dv^2$')
    elif color_type=='dW_dh':
            colors = data["dW_dh"]
            x = [v[0] for v in ver_list]
            y = [v[1] for v in ver_list]
            v_area = 40*np.ones(len(colors))
            sc = ax.scatter(x, y, s=v_area, c=colors, norm=mcolors.LogNorm())#, alpha=0.0)
            clb = plt.colorbar(sc)
            clb.ax.set_title('$dW/dv$')

    # Turn off tick labels
    ax.set_yticklabels([])
    ax.set_xticklabels([])

    #cbaxes = fig.add_axes([0.8, 0.1, 0.03, 0.8])
    #clb = plt.colorbar(p, shrink=0.5, cax = cbaxes)
    #clb = plt.colorbar(p,shrink=0.5)
    #clb.ax.set_title('cell area')
    plt.xlim([-0.11,0.05+data["box_size"][0]])
    plt.ylim([-0.11,0.05+data["box_size"][1]])
    plt.savefig(image_name, dpi=300)
    plt.close()
