#!/usr/bin/env python3
import sys
import os
from Tissue_Objects import *
from CommonOperations import *

import argparse

parser=argparse.ArgumentParser(
    description='''This code simulates 2D vertex model with periodic boundary condition.
    In this model the vertices move dynamically down the gradient of a work function W.
    mu dX/dt = -dW/dX. The simulation runs for given "n_steps" steps.''',
    epilog="""run with the correct arguments. """)
parser.add_argument('sim_type', type=str, default='randomized_tension', help='choose from: "probing_vertices", "randomized_tension", "quadruple", "3fold", "4fold", and "randomized_packing"')
parser.add_argument('Nxy', type=int, default=8, help='the given integer value will be number of cells in both x and y direction. Please choose an even number!')
parser.add_argument('gamma', type=float, default=0.01, help='average line tension')
parser.add_argument('n_steps', type=int, default=5000, help='number of simulation time steps')
parser.add_argument('dt', type=float, default=0.1, help='size of time steps')
parser.add_argument('inserting_gamma', type=float, default=0.1, help='line tension for the inserting cells')
args=parser.parse_args()

if not os.path.exists("data"):
    os.makedirs("data")
if not os.path.exists("images"):
    os.makedirs("images")

#if len(sys.argv)==1:
#    print("Please provide an input file.")
#    exit(1)

def visulaise_frames(n_steps, skip, color_type):
    for i in np.arange(0,int(n_steps+1), skip):
        data_name = "data/T"+str(i)+".json"
        image_name = "images/t"+str(i).zfill(5)+".png"
        visualize_frame_jsonfile(data_name, image_name, color_type)

def run_simulation(Nx, Ny, Ta, Kc, A0c, Gc, sim_type, n_steps, dt, inserting_gamma):
    """
    First a regular hexagonal tissue with given arguments is creatred.
    Next: it runs for 1000 steps to move tissue near force balance,
    before I start the actual simulation mode, which is given as "sim_type".
    """
    print("inserting_gamma:"+str(inserting_gamma))
    if not os.path.exists(str(inserting_gamma)+"/data"):
        os.makedirs(str(inserting_gamma)+"/data")
    if not os.path.exists(str(inserting_gamma)+"/images"):
        os.makedirs(str(inserting_gamma)+"/images")

    T = Tissue(Ta, Kc, A0c, Gc)

    if sim_type=='edge_collapse':
        T.build_tissue_from_file("data/T_final.json")

        delta_gamma = 0.2 * T.Ta
        el_ratio = []
        e_lengths = []
        e_tensions = []
        for i in range(len(T.ListEdge)):
            T.is_with_T1_trans = False
            e = T.ListEdge[i]
            e.get_dl()
            e_i = e.length
            e_lengths.append(e_i)
            e_tensions.append(e.lineTension)
            e.lineTension = e.lineTension + delta_gamma
            T.minimize_dynamically(3000, 0.1, 0, False, False)
            e.get_dl()
            e_f = e.length
            print(i, e_i, e_f, e_f/e_i)
            el_ratio.append(e_f/e_i)
            T.build_tissue_from_file("data/T_final.json")
        np.savetxt("edges_collapse.txt", el_ratio, fmt="%.3f")
        np.savetxt("edges_lengths.txt", e_lengths, fmt="%.3f")
        np.savetxt("edges_tensions.txt", e_tensions, fmt="%.3f")
        return T, []

    T.create_hexagonal_tissue(Nx, Ny, Kc, A0c, Ta, Gc)
    T.update_derivatives_analytically()

    T.is_with_T1_trans = True
    T.minimize_dynamically(1000, 0.1, 0, False, False)

    delta = 0.01*np.sqrt(T.A0c) # opening length of edge
    write_freq = 1000
    W_steps = []

    if sim_type=='4fold':
        v_cent = T.find_vertex_at_center()
        shrinking_edge = v_cent.connectedEdges[0]
        T.shrink_edge_to_vertex(shrinking_edge)
        T.minimize_dynamically(2000, 0.1, write_freq, False, False)
        for v in T.ListVertex:
            if len(v.connectedEdges)==4:
                T.open_vertex_to_cell(v)
        W_steps = T.minimize_dynamically(n_steps, dt, write_freq, True, True)

    elif sim_type=='3fold':
        v_cent = T.find_vertex_at_center()
        T.open_vertex_to_cell(v_cent)
        W_steps = T.minimize_dynamically(n_steps, dt, write_freq, True, True)

    elif sim_type=='both_sides':
        v_cent = T.find_vertex_at_center()
        e = v_cent.connectedEdges[0]
        T.open_vertex_to_cell(e.tailVertex)
        T.open_vertex_to_cell(e.headVertex)
        W_steps = T.minimize_dynamically(n_steps, dt, write_freq, True, True)

    elif sim_type=='randomized_packing':
        for i in range(10):
            n_e = len(T.ListEdge)
            eid = random.randrange(n_e)
            e = T.ListEdge[eid]

            if e.crossBdry or e.c1.crossBdry or e.c2.crossBdry:
                continue
            if len(e.c1.ListEdges)>4 and len(e.c2.ListEdges)>4:
                T.shrink_edge_to_vertex(e)
            #T.shrink_random_edge()
        for v in T.ListVertex:
            dr = Point(random.uniform(-1.0,1.0), random.uniform(-1.0,1.0), 0)
            v.coord = v.coord +dr.X(0.05)
        W_steps = T.minimize_dynamically(n_steps, dt, write_freq, True, True)

    elif sim_type=='quadruple':
        v_cent = T.find_vertex_at_center()
        shrinking_edge = v_cent.connectedEdges[0]
        T.flip_edge(shrinking_edge)
        W_steps = T.minimize_dynamically(n_steps, dt, write_freq, True, True)
    elif sim_type=='randomized_tension':
        for i in range(int(0.2*len(T.ListEdge))):
            n_e = len(T.ListEdge)
            eid = random.randrange(n_e)
            e = T.ListEdge[eid]
            tx, ty = e.tailVertex.coord.x, e.tailVertex.coord.y
            if is_point_near_bdr(tx, ty, T.SystemSize.x, T.SystemSize.y):
                continue

            if e.crossBdry or e.c1.crossBdry or e.c2.crossBdry:
                continue
            if len(e.c1.ListEdges)>4 and len(e.c2.ListEdges)>4:
                T.flip_edge(e)
        T.minimize_dynamically(10000, dt, write_freq, False, False)
        gamma_mean = T.Ta
        gamma_std = T.Ta/4.0
        Tvalues = np.random.normal(gamma_mean, gamma_std, len(T.ListEdge))
        for i in range(len(Tvalues)):
            if T.ListEdge[i].c1.crossBdry or T.ListEdge[i].c2.crossBdry:
                continue
            T.ListEdge[i].lineTension = Tvalues[i]
        W_steps = T.minimize_dynamically(n_steps, dt, write_freq, True, True)

    elif sim_type=='probing_vertices':
        for i in range(int(0.15*len(T.ListEdge))):
            n_e = len(T.ListEdge)
            eid = random.randrange(n_e)
            e = T.ListEdge[eid]
            tx, ty = e.tailVertex.coord.x, e.tailVertex.coord.y
            if is_point_near_bdr(tx, ty, T.SystemSize.x, T.SystemSize.y):
                continue

            if e.crossBdry or e.c1.crossBdry or e.c2.crossBdry:
                continue
            if len(e.c1.ListEdges)>4 and len(e.c2.ListEdges)>4:
                T.flip_edge(e)
        T.minimize_dynamically(1000, dt, write_freq, False, False)
        gamma_mean = T.Ta
        gamma_std = T.Ta/5.0
        Tvalues = np.random.normal(gamma_mean, gamma_std, len(T.ListEdge))
        for i in range(len(Tvalues)):
            if T.ListEdge[i].c1.crossBdry or T.ListEdge[i].c2.crossBdry:
                continue
            T.ListEdge[i].lineTension = Tvalues[i]
        W_steps = T.minimize_dynamically(10000, dt, write_freq, True, True)

        old_coords = []
        for v in T.ListVertex:
            old_coords.append(v.coord)

        print("number of vertices: "+str(len(T.ListVertex)))
        for v in T.ListVertex:
            T.probe_vertex(v, 500, old_coords)
            prt_str = "vertex id: "+str(v.storing_order)+", k_d="+"{:.3f}".format(v.probing_k)
            prt_str = prt_str + ", polyclass of neighb cells: "
            for cor in v.ListCorners:
                prt_str = prt_str + str(len(cor.c.ListEdges))+" "
            print(prt_str)
        T.write_tissue_to_file('data/T_final.json', True)
        datafile = 'data/T_final.json'
        color_type = 'dW_dh'
        image_name = 'gamma'+str(Ta)+color_type+'.pdf'
        visualize_frame_jsonfile(datafile, image_name, color_type)
        color_type = 'probing'
        image_name = 'gamma'+str(Ta)+color_type+'.pdf'
        visualize_frame_jsonfile(datafile, image_name, color_type)

        color_type = 'd2W_dv2'
        image_name = 'gamma'+str(Ta)+color_type+'.pdf'
        visualize_frame_jsonfile(datafile, image_name, color_type)
        visulaise_frames(n_steps, write_freq, color_type)
    elif sim_type=='check_vertex_opening':
        #now we check the vertex opening
        print("number of vertices: "+str(len(T.ListVertex)))
        area_opening = []
        T.build_tissue_from_file("data/T_final.json")
        for vid in range(len(T.ListVertex)):
            v = T.ListVertex[vid]
            prt_str = "vertex id: "+str(v.storing_order)
            prt_str = prt_str + ", polyclass of neighb cells: "
            for cor in v.ListCorners:
                prt_str = prt_str + str(len(cor.c.ListEdges))+" "

            T.check_vertex_opening(v, 1000, inserting_gamma)
            prt_str = prt_str +", a_f="+"{:.3f}".format(v.final_opening_area)
            print(prt_str)
            area_opening.append(v.final_opening_area)
            T.build_tissue_from_file("data/T_final.json")

        for iv in range(len(T.ListVertex)):
            T.ListVertex[iv].final_opening_area = area_opening[iv]

        T.write_tissue_to_file(str(inserting_gamma)+'/data/rebuilt_T_final.json', True)
        datafile = str(inserting_gamma)+'/data/rebuilt_T_final.json'
        color_type = 'opening_area'
        image_name = str(inserting_gamma)+'/images/inserting_gamma'+str(inserting_gamma)+color_type+'.pdf'
        make_histogram_plots(str(inserting_gamma)+'/data/', str(inserting_gamma)+'/images/')


    #datafile = 'data/T_final.json'
    #color_type = 'dW_dh'
    #image_name = 'gamma'+str(Ta)+color_type+'.pdf'
    #visualize_frame_jsonfile(datafile, image_name, color_type)
    #color_type = 'probing'
    #image_name = 'gamma'+str(Ta)+color_type+'.pdf'
    #visualize_frame_jsonfile(datafile, image_name, color_type)
    #color_type = 'd2W_dv2'
    #image_name = 'gamma'+str(Ta)+color_type+'.pdf'
    #visualize_frame_jsonfile(datafile, image_name, color_type)
    #visulaise_frames(n_steps, write_freq, color_type)
    return T, W_steps

sim_type = sys.argv[1]

Nxy = int(sys.argv[2])

gamma = float(sys.argv[3])

n_steps = int(sys.argv[4])

dt = float(sys.argv[5])

inserting_gamma = float(sys.argv[6])

Kc, A0c, Gc = 1.0, 1.0, 0

run_simulation(Nxy, Nxy, gamma, Kc, A0c, Gc, sim_type, n_steps, dt, inserting_gamma)
