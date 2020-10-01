#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Friday May29 2020

@email ---  aamiri@pks.mpg.de
"""

import numpy as np
from math import floor
import random
import json
from collections import namedtuple

#=============================
#------ some data types ------
class Point:
    """A class that defines points generally in 3D: Point(x,y,z)
    I use these objects to define positin, froces, and other vector quantities"""
    def __init__(self, x = 0, y = 0, z=0):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        """This allow to print the object"""
        return "({0},{1},{2})".format(self.x,self.y,self.z)

    def __add__(self,other):
        """overloading a + operator to add two point objects"""
        x = self.x + other.x
        y = self.y + other.y
        z = self.z + other.z
        return Point(x,y,z)

    def __sub__(self,other):
        """overloading a - operator to subtract two point objects"""
        x = self.x - other.x
        y = self.y - other.y
        z = self.z - other.z
        return Point(x,y,z)

    def __truediv__(self,m):
        """overloading a / operator to divide a point object by a scalr m"""
        x = self.x/m
        y = self.y/m
        z = self.z/m
        return Point(x,y,z)

    def __mul__(self,other):
        """overloading a * operator to cross product vectors: self x other"""
        x = self.y * other.z - other.y * self.z
        y = -self.x *other.z + self.z * other.x
        z = self.x * other.y - self.y * other.x
        return Point(x,y,z)

    def X(self,m):
        """A function to mulitply elements of the object by a scalar m"""
        x = self.x * m
        y = self.y * m
        z = self.z * m
        return Point(x,y,z)

    def dot(self,other):
        """A function for dot product of two objects self . other"""
        return self.x*other.x + self.y*other.y + self.z*other.z

    def Norm(self):
        """Norm of the object"""
        return(np.sqrt(self.x*self.x + self.y*self.y + self.z*self.z))

class Quadrant:
    """A class to define Quadrant objects.
    I use these objects for stablishing periodicity."""
    def __init__(self, x = 0, y = 0):
        self.x = x
        self.y = y

    def __str__(self):
        return "({0},{1})".format(self.x,self.y)

    def __add__(self,other):
        x = self.x + other.x
        y = self.y + other.y
        return Quadrant(x,y)

    def __sub__(self,other):
        x = self.x - other.x
        y = self.y - other.y
        return Quadrant(x,y)

    def Norm(self):
        return(np.sqrt(self.x*self.x + self.y*self.y))

class Corner:
    """A n fold vertex has n corners.
    the corner belongs to a cell c.
    at this corner, there is an edge before the vertex ep with the direction dp,
    and another edge after en with the direction dn."""
    def __init__(self, c = [], ep = [], dp = True, en = [], dn = True):
        self.c = c
        self.ep = ep
        self.dp = dp
        self.en = en
        self.dn = dn

    def is_corner_neighbor_with(self,other):
        """wehter two corners are neighbors, in other words do they share an edge in between them?"""
        if self.ep==other.ep or self.ep==other.en:
            return True
        elif self.en==other.ep or self.en==other.en:
            return True
        else:
            return False


#====== Some functions
def add_point_if_not_there(list_points, p):
    """adds point p to the list of points if already not in the list."""
    does_exist=False
    for op in list_points:
        delta_p = op-p
        dist = delta_p.Norm()
        if dist<1e-4:
            does_exist=True
    if does_exist==False:
        list_points.append(p)
    return list_points

def fit_in_periodic_box(p, Lx, Ly):
    """by looking at the coordinates of point p and comparing it with the size of
    the periodic box, decides how the point should be added."""
    p.x, p.y = p.x/Lx, p.y/Ly
    pqx, pqy=floor(p.x), floor(p.y)
    p.x, p.y=(p.x-pqx)*Lx, (p.y-pqy)*Ly
    return p, Quadrant(pqx,pqy)


# ==============================================================
# --------  HERE I define class for vertex, edge, triangl, Cell, etc ---------
# ==============================================================

# ----- Vertex class ------
#--------------------------
class Vertex:
    """ Vertex objects. initiated with an id: Number, type:vertexType, and it's coordinate
    Force on the vertex: dW_dv of type Point.
    List of all edges connected edges to this vertex: connectedEdges
    List all corners around the vertex: ListCorners
    second order response to out of plane smalle deformation: d2W_dh2"""
    def __init__(self, Number, vertexType, coord):
        self.Number = Number
        self.Type = vertexType
        self.coord = coord
        self.dW_dv = Point(0,0,0)
        self.storing_order = -1
        self.connectedEdges = []
        self.ListCorners=[]
        self.d2W_dh2 = 0
        self.Norm_dW_dh = 0
        self.last_tension = 0#if an edge shrinks to become 4-fold vertex, this vertex will remember the tension of that edge
        self.probing_k = 0
        self.final_opening_area = 0

    def set_force_to_zero(self):
        self.dW_dv = Point(0,0,0)

    def connected_edge_after(self, e):
        left_c = []
        if e.tailVertex==self:
            left_c = e.c1
        elif e.headVertex==self:
            left_c = e.c2

        for ite in self.connectedEdges:
            if ite==e:
                continue
            if ite.tailVertex==self and ite.c2==left_c:
                return ite
            elif ite.headVertex==self and ite.c1==left_c:
                return ite
    def order_connected_edges(self):
        ordered_list=[]
        this_e = self.connectedEdges[0]
        ordered_list.append(this_e)
        while len(ordered_list) <len(self.connectedEdges):
            this_e = self.connected_edge_after(this_e)
            ordered_list.append(this_e)
        self.connectedEdges = ordered_list

    def get_non_neighboring_corner_pairs(self):
        if len(self.ListCorners)<4:
            return
        pair1, pair2 = [], []
        cor0 = self.ListCorners[0]
        cor1 = self.ListCorners[1]
        cor2 = self.ListCorners[2]
        cor3 = self.ListCorners[3]
        if cor0.is_corner_neighbor_with(cor1)==False:
            pair1 = [cor0, cor1]
            pair2 = [cor2, cor3]
        elif cor0.is_corner_neighbor_with(cor2)==False:
            pair1 = [cor0, cor2]
            pair2 = [cor1, cor3]
        elif cor0.is_corner_neighbor_with(cor3)==False:
            pair1 = [cor0, cor3]
            pair2 = [cor1, cor2]

        return pair1, pair2

    def W_corner(self):
        Wv = 0
        for cor in self.ListCorners:
            cor.c.update_geometric_info()
            Wv = Wv + 0.5 * cor.c.Kc * (cor.c.area - cor.c.A0c)* (cor.c.area - cor.c.A0c)
            Wv = Wv + 0.5 * cor.c.Gc * cor.c.perimeter * cor.c.perimeter
        for e in self.connectedEdges:
            Wv = Wv + e.lineTension * e.length
        return Wv

    def estimate_second_derivative(self, delta):
        old_coord = self.coord
        dr = Point(0, 0, delta)
        dr2 = Point(0, 0, 2.0*delta)
        W0 = self.W_corner()
        self.coord = old_coord + dr
        Wf = self.W_corner()
        self.coord = old_coord + dr2
        Wff = self.W_corner()
        self.coord = old_coord # back to its old coordinates
        self.Norm_dW_dh = (Wf-W0)/delta
        self.d2W_dh2 = (Wff - 2.0*Wf + W0)/(delta*delta)
        return self.d2W_dh2


# ----- Edge class ------
#--------------------------
class Edge:
    """ Edge object has a tail vertex and a head vertex.
    it is oriented counter-clock-wise around cell c1 and clockwise around cell c2."""
    def __init__(self, Number, edgeType, tailVertex, headVertex, lineTension, q_tail, q_head, T):
        self.Number = Number
        self.Type = edgeType
        self.tailVertex = tailVertex
        self.headVertex = headVertex
        self.lineTension = lineTension
        self.q_tail = q_tail
        self.q_head = q_head
        self.T = T
        self.storing_order = -1
        self.crossBdry = False
        self.c1 = []
        self.c2 = []
        self.length = 0
        self.dl = self.get_dl()

    def replace_vertex_with(self, v_old, v_new, q_new):
        if self.tailVertex==v_old:
            self.tailVertex = v_new
            self.q_tail = q_new
        elif self.headVertex==v_old:
            self.headVertex = v_new
            self.q_head = q_new

    def get_dl(self):
        h = self.headVertex.coord + Point(self.q_head.x * self.T.SystemSize.x, self.q_head.y * self.T.SystemSize.y, 0)
        t = self.tailVertex.coord + Point(self.q_tail.x * self.T.SystemSize.x, self.q_tail.y * self.T.SystemSize.y, 0)

        self.dl = h-t
        self.length = self.dl.Norm()
        return self.dl

    def get_vertex_opening_along(self, v_coord, v, delta):
        e_dl = self.dl
        e_dl = e_dl/(e_dl.Norm())
        dl_opening = e_dl.X(delta)
        if self.tailVertex==v:
            return v_coord+dl_opening
        else:
            return v_coord-dl_opening

# ----- Cell class ------
#--------------------------
class Cell:
    """ListEdgesOrder makes edges either CW or CCW. Then Each face will be assigned
    to the cell with an ORDER to make it CCW around outward normal to cell"""
    def __init__(self, Number, cType, Kc, A0c, Ta, Gc):
        self.Number = Number
        self.Type = cType
        self.Kc = Kc
        self.A0c = A0c
        self.Ta = Ta
        self.Gc = Gc
        self.ListEdges = []
        self.ListVertices = []
        self.ListEdgesOrder = []
        self.center = 0
        self.area = 0
        self.perimeter = 0
        self.crossBdry = False
        self.storing_order = -1

    def remove_edge(self, e):
        #print('cell:'+str(self.Number)+'--------   Ne before removing: '+str(len(self.ListEdges)))
        for E in self.ListEdges:
            this_e = E[0]
            if this_e==e:
                self.ListEdges.remove(E)
                #print('cell:'+str(self.Number)+'--------   Ne AFTER removing: '+str(len(self.ListEdges)))
                return


    def remove_vertex(self, v):
        #print('cell:'+str(self.Number)+'--------   Nv before removing: '+str(len(self.ListVertices)))
        if v in self.ListVertices:
            self.ListVertices.remove(v)
        #print('cell:'+str(self.Number)+'--------   Nv AFTER removing: '+str(len(self.ListVertices)))

    def add_vertex(self, v):
        if v in self.ListVertices:
            return
        else:
            self.ListVertices.append(v)

    def add_edge(self, e, d):
        does_exist = False
        for E in self.ListEdges:
            E0, E1 = E[0], E[1]
            if e==E0:
                does_exist=True
        if does_exist==False:
            self.ListEdges.append((e,d))

    def getEdgeBeforeAfter(self, v):
        e_p, d_p, e_n, d_n = [], [], [], []
        for E in self.ListEdges:
            e, d = E[0], E[1]
            if e.tailVertex==v:
                if d:
                    e_n, d_n = e, True
                else:
                    e_p, d_p = e, False
            if e.headVertex==v:
                if d:
                    e_p, d_p = e, True
                else:
                    e_n, d_n = e, False

        return e_p, d_p, e_n, d_n

    def get_next_edge(self, E_given):
        e_given, dir_given = E_given[0], E_given[1]
        for E in self.ListEdges:
            if e_given==E[0]:
                continue
            if dir_given and E[1] and E[0].tailVertex==e_given.headVertex:
                next_E = E
                return next_E
            elif dir_given==False and E[1] and E[0].tailVertex==e_given.tailVertex:
                next_E = E
                return next_E
            elif dir_given and E[1]==False and E[0].headVertex==e_given.headVertex:
                next_E = E
                return next_E
            elif dir_given==False and E[1]==False and E[0].headVertex==e_given.tailVertex:
                next_E = E
                return next_E

    def order_edge_list(self):
        ordered_list = []
        E0 = []
        for E in self.ListEdges:
            if E[0].crossBdry==False:
                E0 = E
                ordered_list.append(E)
                break
        E_this = E0
        for e_id in range(len(self.ListEdges)-1):
            E_next = self.get_next_edge(E_this)
            ordered_list.append(E_next)
            E_this = E_next
        self.ListEdges = ordered_list
        #now let's also order the verticess
        ordered_verts = []
        for E in self.ListEdges:
            if E[1]:
                ordered_verts.append(E[0].tailVertex)
            else:
                ordered_verts.append(E[0].headVertex)
        self.ListVertices = ordered_verts

    def update_corners(self):
        for v in self.ListVertices:
            e_p, d_p, e_n, d_n = self.getEdgeBeforeAfter(v);
            cor = Corner(self, e_p, d_p, e_n, d_n)
            v.ListCorners.append(cor)

    def update_geometric_info(self):
        self.center = Point(0,0,0)
        self.area = 0
        areaVec = Point(0,0,0)
        self.perimeter = 0
        p0, p_this, p_next = [], [], []
        E0 = self.ListEdges[0]
        e0, d0 = E0[0], E0[1]
        if d0:
            p0=e0.tailVertex.coord
        else:
            p0 = e0.headVertex.coord
        p_this = p0

        self.center = self.center + p_this
        e_ct=0
        for E in self.ListEdges:
            e, d = E[0], E[1]
            e.get_dl()

            self.perimeter = self.perimeter + e.length
            if d:
                p_next = p_this + e.dl
            else:
                p_next = p_this - e.dl
            if e_ct<len(self.ListEdges)-1:
                self.center = self.center + p_next
            areaVec = areaVec + (p_this-p0) * (p_next-p0)
            p_this = p_next
            e_ct = e_ct + 1

        self.center = self.center/len(self.ListEdges)
        self.area = 0.5 * areaVec.dot(Point(0,0,1))


# ----- Tissue class ------
#--------------------------
class Tissue:
    def __init__(self, Ta, Kc, A0c, Gc):
        self.Ta = Ta
        self.Kc = Kc
        self.A0c = A0c
        self.Gc = Gc
        self.energy = 0
        self.isPeriodic = False
        self.ListCell=[]
        #self.ListFace = []
        self.ListVertex = []
        self.ListEdge = []
        self.delta_l = 0.1*np.sqrt(self.A0c)
        self.is_with_T1_trans = False

    def max_cell_id(self):
        max_id = 0
        for c in self.ListCell:
            if c.Number>max_id:
                max_id = c.Number
        return max_id

    def update_vertex_connections(self):
        v_storing_order = 0
        for v in self.ListVertex:
            v.storing_order = v_storing_order
            v.connectedEdges = []
            v.ListCorners = []
            v_storing_order = v_storing_order + 1

        e_storing_order=0
        for e in self.ListEdge:
            e.storing_order = e_storing_order
            e.tailVertex.connectedEdges.append(e)
            e.headVertex.connectedEdges.append(e)
            e_storing_order = e_storing_order + 1

        c_storing_order=0
        for c in self.ListCell:
            c.order_edge_list()
            c.update_corners()
            c.storing_order = c_storing_order
            c_storing_order= c_storing_order + 1

    def W(self):
        self.energy = 0
        for c in self.ListCell:
            c.update_geometric_info()
            self.energy = self.energy + 0.5 * c.Kc * (c.area - c.A0c)* (c.area - c.A0c)
            self.energy = self.energy + 0.5 * c.Gc * c.perimeter * c.perimeter

        for e in self.ListEdge:
            #e.get_dl()
            self.energy = self.energy + e.lineTension*e.length #+ 0.5*c.Gc * c.perimeter*c.perimeter

        return self.energy

    def shrink_edge_to_vertex(self, e):
        if len(e.tailVertex.connectedEdges)>3 or len(e.headVertex.connectedEdges)>3:# at this point we are avoinding 5-fold vertices
            return
        #tail vertex is being removed
        v_tail, v_head = e.tailVertex, e.headVertex
        v_head.last_tension = e.lineTension
        tail_corners = v_tail.ListCorners
        tail_connected_edges = v_tail.connectedEdges
        c1, c2 = e.c1, e.c2
        c3=[]
        for cor in tail_corners:
            if cor.c==c1 or cor.c==c2:
                continue
            else:
                c3 = cor.c
        if len(c1.ListEdges)<4 or len(c2.ListEdges)<4:#it is already a triangle, and can not lose another edge
            return
        midCoord = (v_tail.coord+v_head.coord)/2
        #e.tailVertex.coord = midCoord
        e.headVertex.coord = midCoord
        c1.remove_edge(e)
        c2.remove_edge(e)
        c1.remove_vertex(v_tail)
        c2.remove_vertex(v_tail)
        c3.remove_vertex(v_tail)
        c3.add_vertex(v_head)
        #print("at the end c3 has:"+str(len(c3.ListVertices))+" , "+str(len(c3.ListEdges)))

        edges_c1, edges_c2 = c1.ListEdges, c2.ListEdges

        for connected_e in tail_connected_edges:
            if connected_e==e:
                continue
            connected_e.replace_vertex_with(v_tail, v_head, e.q_head)
        self.ListEdge.remove(e)
        self.ListVertex.remove(v_tail)
        self.update_vertex_connections()
        for c in self.ListCell:
            c.update_geometric_info()

        return v_head
            #print(len(c.ListEdges), len(c.ListVertices), c.area)
        #print("==============================")

    def shrink_random_edge(self):
        n_e = len(self.ListEdge)
        success = False
        while success==False:
            eid = random.randrange(n_e)
            e = self.ListEdge[eid]
            if e.c1.crossBdry==False and e.c2.crossBdry==False:
                self.shrink_edge_to_vertex(e)
                success = True


    def shrink_small_edges(self):
        for e in self.ListEdge:
            if e.crossBdry:
                continue

            if e.c1.crossBdry==True or e.c2.crossBdry==True:
                continue

            if len(e.tailVertex.connectedEdges)>3 or len(e.headVertex.connectedEdges)>3:
                continue# we still don't handle 5-fold vertices

            if e.length<self.delta_l:
                self.shrink_edge_to_vertex(e)

    def force_corner_at_vertex(self, cor, v):
        p_prev, p_next = [], []
        if cor.ep.headVertex==v:
            p_prev = v.coord - cor.ep.dl
        elif cor.ep.tailVertex==v:
            p_prev = v.coord + cor.ep.dl

        if cor.en.headVertex==v:
            p_next = v.coord - cor.en.dl
        elif cor.en.tailVertex==v:
            p_next = v.coord + cor.en.dl

        fv = Point(p_next.y - p_prev.y, p_prev.x - p_next.x, 0.0)
        f_cor = fv.X( -0.5 * cor.c.Kc *(cor.c.area - cor.c.A0c) )
        dl_prev, dl_next =  p_prev-v.coord, p_next-v.coord
        dl_prev, dl_next = dl_prev/dl_prev.Norm(), dl_next/dl_prev.Norm()
        f_cor = f_cor + dl_prev.X(cor.ep.lineTension)
        f_cor = f_cor + dl_next.X(cor.en.lineTension)

        return f_cor

    def check_opening_forces_four_fold(self, v):
        #print(v.coord)
        pair1, pair2 = v.get_non_neighboring_corner_pairs()
        forces1 = [self.force_corner_at_vertex(pair1[0], v), self.force_corner_at_vertex(pair1[1], v)]
        forces2 = [self.force_corner_at_vertex(pair2[0], v), self.force_corner_at_vertex(pair2[1], v)]
        #print('at v:'+str(v.Number))
        #print(pair1[0].dp, pair1[0].ep.c1.Number,pair1[0].ep.c2.Number, pair1[0].dn, pair1[0].en.c1.Number,pair1[0].en.c2.Number)
        #print(pair1[1].dp, pair1[1].ep.c1.Number,pair1[1].ep.c2.Number, pair1[1].dn, pair1[1].en.c1.Number,pair1[1].en.c2.Number)
        #print(pair2[0].dp, pair2[0].ep.c1.Number,pair2[0].ep.c2.Number, pair2[0].dn, pair2[0].en.c1.Number,pair2[0].en.c2.Number)
        #print(pair2[1].dp, pair2[1].ep.c1.Number,pair2[1].ep.c2.Number, pair2[1].dn, pair2[1].en.c1.Number,pair2[1].en.c2.Number)
        f1, f2 = forces1[1]-forces1[0], forces2[1]-forces2[0]
        #print('f1: '+str(f1.Norm())+ '  f2: '+str(f2.Norm()))

        return pair1, forces1, pair2, forces2

    def prospect_edge_direction_at_corner(self, e, cor):
        #print(e.tailVertex.Number, e.headVertex.Number)
        #print(cor.ep.tailVertex.Number, cor.ep.headVertex.Number)
        #print(cor.en.tailVertex.Number, cor.en.headVertex.Number)
        #print("--------")
        if cor.dp:
            if cor.ep.headVertex==e.tailVertex:
                return True
            elif cor.ep.headVertex==e.headVertex:
                return False

        if cor.dp==False:
            if cor.ep.tailVertex==e.tailVertex:
                return True
            elif cor.ep.tailVertex==e.headVertex:
                return False

        return True

    def open_vertex_to_edge(self, v, separating_corners, expanding_corners, extra_l):
        """ Separating corners move apart from each other along the new edge.
        expanding corners gain the new edge.
        the vertex v will serve as tail vertex of the new edge.
        extra_l is usually a small addition to avoid closing an edge right after its opening."""
        delta = self.delta_l + extra_l
        opening_tenion = v.last_tension
        vcoord = v.coord
        dl_new = separating_corners[0].c.center - separating_corners[1].c.center
        dl_new = dl_new.X(delta/(2.0*dl_new.Norm()))
        tail_coord = vcoord - dl_new
        head_coord = vcoord + dl_new
        v.coord = tail_coord
        #now let's create the new vertex:
        new_vertex = self.similar_or_new_vertex(head_coord, 1)
        separating_corners[0].c.add_vertex(new_vertex)
        separating_corners[0].c.remove_vertex(v)
        expanding_corners[0].c.add_vertex(new_vertex)
        expanding_corners[1].c.add_vertex(new_vertex)
        separating_corners[0].ep.replace_vertex_with(v, new_vertex, Quadrant(0,0))
        separating_corners[0].en.replace_vertex_with(v, new_vertex, Quadrant(0,0))
        new_edge = self.similar_or_new_edge(v, new_vertex, opening_tenion, Quadrant(0,0), Quadrant(0,0))
        #new_edge = self.similar_or_new_edge(v, new_vertex, self.Ta, Quadrant(0,0), Quadrant(0,0))
        #print('length of the new edge is equal to: '+str(new_edge[0].length))
        new_edge[0].crossBdry = False
        d_e = self.prospect_edge_direction_at_corner(new_edge[0], expanding_corners[0])
        #print(d_e)
        if d_e:
            new_edge[0].c1=expanding_corners[0].c
            new_edge[0].c2=expanding_corners[1].c
        else:
            new_edge[0].c1=expanding_corners[1].c
            new_edge[0].c2=expanding_corners[0].c
        expanding_corners[0].c.add_edge(new_edge[0], d_e)
        d_e = self.prospect_edge_direction_at_corner(new_edge[0], expanding_corners[1])
        #print(d_e)
        if d_e:
            new_edge[0].c1=expanding_corners[1].c
            new_edge[0].c2=expanding_corners[0].c
        else:
            new_edge[0].c1=expanding_corners[0].c
            new_edge[0].c2=expanding_corners[1].c
        expanding_corners[1].c.add_edge(new_edge[0], d_e)


    def open_four_fold_vertices(self):
        """
        from a set of all 4-fold vertices, we pick one vertex at random.
        Then we find the preferred opening direction, and open it into
        a new edge.
        A 4-fold vertex has been generated as a previous edge shrinking.
        The vertex rememebrs the line tension value of that edge v.last_tension
        We use this last_tension as tension value for the opening new edge."""
        v_ct = 0
        vlist_4fold = []
        for v in self.ListVertex:
            if len(v.connectedEdges)>3:
                vlist_4fold.append(v_ct)
            v_ct = v_ct + 1

        n_v = len(vlist_4fold)
        if n_v<1: # we have not found any 4-fold vertex to open it into edge
            return

        vid_l = random.randrange(n_v)# we pich
        opening_id = vlist_4fold[vid_l]
        opening_v = self.ListVertex[opening_id]

        pair1, forces1, pair2, forces2 = self.check_opening_forces_four_fold(opening_v)

        df1 = forces1[0] - forces1[1]
        df2 = forces2[0] - forces2[1]

        if df1.Norm()<0.25*opening_v.last_tension and df2.Norm()<0.25*opening_v.last_tension:
            return

        extra_l = 0.01*np.sqrt(self.A0c)
        if random.uniform(0,1)>0.15 and df1.Norm()>df2.Norm():# the first condition is to add some stochasticity
            self.open_vertex_to_edge(opening_v, pair1, pair2, extra_l)
            self.update_vertex_connections()
            for c in self.ListCell:
                c.update_geometric_info()
            return
        else:
            self.open_vertex_to_edge(opening_v, pair2, pair1, extra_l)
            self.update_vertex_connections()
            for c in self.ListCell:
                c.update_geometric_info()
            return


    def flip_edge(self, e):
        v = self.shrink_edge_to_vertex(e)
        if len(v.connectedEdges)<4:
            print("Warning: openeing a vertex that is not 4-fold")
        pair1, forces1, pair2, forces2 = self.check_opening_forces_four_fold(v)
        df1 = forces1[0] - forces1[1]
        df2 = forces2[0] - forces2[1]
        extra_l = 0.3*np.sqrt(self.A0c)
        if df1.Norm()<df2.Norm():
            self.open_vertex_to_edge(v, pair1, pair2, extra_l)
        else:
            self.open_vertex_to_edge(v, pair2, pair1, extra_l)

        self.update_vertex_connections()
        for c in self.ListCell:
            c.update_geometric_info()

    def divide_cell(self, c):
        """ Not completed yet!
        dividing cell c with a random axis """
        e1_i = int((c.ListEdges)/2)-1
        e2_i = e1_i + int((c.ListEdges)/2)
        E1, E2 = c.ListEdges[e1_i], c.ListEdges[e2_i]

    def probe_vertex(self, v0, n_steps, old_coords):
        """here we want to apply out of plane force to vertex v.
        only vertex v can move out of plane, the rest will only
        have in-plane movements."""
        self.is_with_T1_trans = False
        eps_val = 1e-20
        f_out = Point(0,0,-0.01)
        dt, dampCoef = 0.1, 1
        for s in range(n_steps+1):
            self.update_derivatives_analytically()
            v0.dW_dv = v0.dW_dv + f_out
            for v in self.ListVertex:
                velocity = v.dW_dv/(-dampCoef)
                dr = velocity.X(dt)
                v.coord = v.coord + dr
        v0.probing_k = abs(f_out.z)/(abs(v0.coord.z)+eps_val)

        for i in range(len(old_coords)):#now we move back all vertices to their oroginal position
            self.ListVertex[i].coord = old_coords[i]

    def check_vertex_opening(self, v0, n_steps, inserting_gamma):
        """to check if a vertex can sunccesfully open to a cell"""
        #for e in v0.connectedEdges:
        #    if e.crossBdry:
        #        return
        for cor in v0.ListCorners:
            if cor.c.crossBdry:
                return

        self.open_vertex_to_cell(v0, inserting_gamma)

        self.is_with_T1_trans = False
        dt, dampCoef = 0.1, 1.0
        for s in range(n_steps+1):
            self.update_derivatives_analytically()
            v0.dW_dv = v0.dW_dv
            for v in self.ListVertex:
                velocity = v.dW_dv/(-dampCoef)
                dr = velocity.X(dt)
                v.coord = v.coord + dr
        v0.final_opening_area = self.ListCell[-1].area



    def similar_or_new_vertex(self, vcoord, vtype):
        cut = 1e-6
        for v in self.ListVertex:
            if v.Type != vtype:
                continue
            dp = v.coord - vcoord
            dist = dp.Norm()
            if dist<cut:
                return v
        #if the there doesnt exist a vertex with similar coordinates them create it
        nv = len(self.ListVertex)
        self.ListVertex.append(Vertex(nv, vtype, vcoord))
        return self.ListVertex[-1]

    def similar_or_new_edge(self, v_tail, v_head, ltens, qt, qh):
        for e in self.ListEdge:
            if e.tailVertex== v_tail and e.headVertex== v_head:
                return (e,True)
            elif e.tailVertex== v_head and e.headVertex== v_tail:
                return (e,False)
        # --> if there doesnt exist an edge with similar vertices then create it
        ne = len(self.ListEdge)
        self.ListEdge.append(Edge(ne, 1, v_tail, v_head, ltens, qt, qh, self))
        return (self.ListEdge[-1],True)

    def cell_between_edges(self, e1, e2):
        if e1.c1==e2.c1 or e1.c1==e2.c2:
            return e1.c1
        elif e1.c2==e2.c1 or e1.c2==e2.c2:
            return e1.c2

    def open_vertex_to_cell(self, v, inserting_gamma):
        l_i = 0.015
        if len(v.connectedEdges)==3:
            l_i=0.015
        elif len(v.connectedEdges)==4:
            l_i = 0.00987#0.012
        v_coord = v.coord
        self.ListVertex.remove(v)
        v.order_connected_edges()
        new_vert_coords = []
        #self.delta_l
        for e in v.connectedEdges:
            nv_coord = e.get_vertex_opening_along(v_coord, v, l_i)
            new_vert_coords.append(nv_coord)
        new_cell_listEdges = []
        new_cell = Cell(self.max_cell_id()+1, 1, self.Kc, self.A0c, inserting_gamma, self.Gc)

        list_new_verts = []
        for vn, ec in zip(new_vert_coords, v.connectedEdges):
            v_new = self.similar_or_new_vertex(vn, 1)
            list_new_verts.append(v_new)
            ec.replace_vertex_with(v, v_new, Quadrant(0,0))

        list_new_edges = []
        for i in range(len(list_new_verts)):
            vp_idx, vn_idx = i, i+1
            if i==len(list_new_verts)-1:
                vn_idx = 0
            vp = list_new_verts[vp_idx]
            vn = list_new_verts[vn_idx]

            new_E = self.similar_or_new_edge(vp, vn, inserting_gamma, Quadrant(0,0), Quadrant(0,0))
            new_E[0].crossBdry, new_E[0].c1 = False, new_cell
            new_E[0].c2 = self.cell_between_edges(v.connectedEdges[vp_idx], v.connectedEdges[vn_idx])
            new_E[0].c2.add_vertex(vp)
            new_E[0].c2.remove_vertex(v)
            new_E[0].c2.add_edge(new_E[0], False )
            list_new_edges.append(new_E)

        new_cell.ListVertices =  list_new_verts
        new_cell.ListEdges =    list_new_edges

        for E in new_cell.ListEdges:
            E[0].lineTension = inserting_gamma

        print(len(new_cell.ListEdges), len(new_cell.ListVertices))

        self.ListCell.append(new_cell)
        #for v in new_cell.ListVertices:
        #    print(v.Number)

        for c in self.ListCell:
            c.update_geometric_info()
        self.update_vertex_connections()

        #In case you want to scale the opening cell area
        #A_n = 0.015*0.015
        #new_cell.update_geometric_info()
        #A_c = new_cell.area
        #c_cent = Point(0,0,0)
        #for v in new_cell.ListVertices:
        #    c_cent = c_cent + v.coord
        #c_cent = c_cent/len(new_cell.ListVertices)
        #print("init area: ", A_c)
        #for v in new_cell.ListVertices:
        #    vi = v.coord - c_cent
        #    sf = np.sqrt(A_n/A_c)
            #v.coord = v.coord/sf
        #    v.coord = c_cent + Point(vi.x/sf, vi.y/sf, 0)
        #new_cell.update_geometric_info()
        #print("the are of new cell: ", new_cell.area)


    def update_derivatives_analytically(self):
        [v.set_force_to_zero() for v in self.ListVertex]
        [c.update_geometric_info() for c in self.ListCell]

        for e in self.ListEdge:
            fe = e.dl.X(e.lineTension/e.length)
            e.tailVertex.dW_dv = e.tailVertex.dW_dv - fe
            e.headVertex.dW_dv = e.headVertex.dW_dv + fe

        # -- I add the forces from cells that contain each vertex
        for v in self.ListVertex:
            p_prev, p_next = [], []
            for cor in v.ListCorners:
                if cor.ep.headVertex==v:
                    p_prev = v.coord - cor.ep.dl
                elif cor.ep.tailVertex==v:
                    p_prev = v.coord + cor.ep.dl

                if cor.en.headVertex==v:
                    p_next = v.coord - cor.en.dl
                elif cor.en.tailVertex==v:
                    p_next = v.coord + cor.en.dl

                fv = Point(p_next.y - p_prev.y, p_prev.x - p_next.x, 0.0)
                v.dW_dv = v.dW_dv + fv.X( 0.5 * cor.c.Kc *(cor.c.area - cor.c.A0c) )

    def find_vertex_at_center(self):
        x_cent, y_cent = self.SystemSize.x/2.0, self.SystemSize.y/2.0
        t_cent = Point(x_cent, y_cent, 0)
        dist = 1000
        for v in self.ListVertex:
            d = (t_cent - v.coord).Norm()
            if d<dist:
                dist = d
                v_cent = v
        return v_cent

    def minimize_dynamically(self, n_steps, dt=0.1, write_freq=10, is_writing=True, is_with_T1_trans=True):
        """to minimize ..."""
        #print(dt, write_freq)
        self.is_with_T1_trans = is_with_T1_trans
        W_steps = []
        dampCoef = 1.0
        for s in range(n_steps+1):
            self.update_derivatives_analytically()
            for v in self.ListVertex:
                velocity = v.dW_dv/(-dampCoef)
                dr = velocity.X(dt)
                v.coord = v.coord + dr
            if is_writing and s%write_freq==0:
                if s>n_steps-2:
                    delta_v = 0.00001*np.sqrt(self.A0c)
                    for v in self.ListVertex:
                        v.estimate_second_derivative(delta_v)
                self.write_tissue_to_file("data/T"+str(int(s))+".json", True)
                W_steps.append(self.W())
            if self.is_with_T1_trans and s%2:
                self.shrink_small_edges()
            if self.is_with_T1_trans and s%4==0:
                self.open_four_fold_vertices()

        self.shrink_small_edges()
        self.write_tissue_to_file("data/T"+str(int(n_steps))+".json", True)
        return W_steps
            #for c in self.ListCell:
            #    c.update_geometric_info()
            #self.update_vertex_connections()

    def create_hexagonal_tissue(self, Nx, Ny, Kc, A0c, Ta, Gc):
        ctype = 1
        r = np.sqrt(2.0*self.A0c/(3.0*np.sqrt(3)))
        lx, ly = np.sqrt(3.0)*r, 1.5*r
        dx = Point(lx,0,0)
        dy = Point(0,ly,0)
        Lx, Ly = Nx*lx,  Ny*ly
        self.SystemSize = Quadrant(Lx, Ly)
        self.isPeriodic = True


        P_o = Point(0,0,0)
        cid = 0
        for nx in range(Nx):
            for ny in range(Ny):
                cell = Cell(cid, ctype, Kc, A0c, Ta, Gc)
                self.ListCell.append(cell)

                displacement = dx.X(nx) + dy.X(ny)
                c_center = P_o + displacement
                if ny%2==1:
                    c_center = c_center + Point(lx/2.0, 0, 0)

                rcos30 = r*np.sqrt(3.0)/2.0
                rsin30 = r/2.0
                #--- p0
                p0= c_center+Point(rcos30, rsin30, 0)
                p0, q0 = fit_in_periodic_box(p0,Lx,Ly)
                #---- p1
                p1 = c_center+Point(0, r, 0)
                p1, q1 = fit_in_periodic_box(p1,Lx,Ly)
                #---- p2
                p2= c_center+Point(-rcos30, rsin30, 0)
                p2, q2 = fit_in_periodic_box(p2,Lx,Ly)
                #---- p3
                p3 = c_center+Point(-rcos30, -rsin30, 0)
                p3, q3 = fit_in_periodic_box(p3,Lx,Ly)
                #---- p4
                p4 = c_center+Point(0, -r, 0)
                p4, q4 = fit_in_periodic_box(p4,Lx,Ly)
                #---- p5
                p5 = c_center+Point(rcos30, -rsin30, 0)
                p5, q5 = fit_in_periodic_box(p5,Lx,Ly)

                totQuad = q0 + q1 + q2 + q3 + q4 + q5
                if totQuad.Norm() >0:
                    cell.crossBdry = True

                v0, v1 = self.similar_or_new_vertex(p0, 1), self.similar_or_new_vertex(p1, 1)
                v2, v3 = self.similar_or_new_vertex(p2, 1), self.similar_or_new_vertex(p3, 1)
                v4, v5 = self.similar_or_new_vertex(p4, 1), self.similar_or_new_vertex(p5, 1)

                e0 = self.similar_or_new_edge(v0, v1, self.Ta, q0, q1)
                if e0[1]:
                    e0[0].c1 = cell
                else:
                    e0[0].c2 = cell

                e1 = self.similar_or_new_edge(v1, v2, self.Ta, q1, q2)
                if e1[1]:
                    e1[0].c1 = cell
                else:
                    e1[0].c2 = cell

                e2 = self.similar_or_new_edge(v2, v3, self.Ta, q2, q3)
                if e2[1]:
                    e2[0].c1 = cell
                else:
                    e2[0].c2 = cell

                e3 = self.similar_or_new_edge(v3, v4, self.Ta, q3, q4)
                if e3[1]:
                    e3[0].c1 = cell
                else:
                    e3[0].c2 = cell

                e4 = self.similar_or_new_edge(v4, v5, self.Ta, q4, q5)
                if e4[1]:
                    e4[0].c1 = cell
                else:
                    e4[0].c2 = cell

                e5 = self.similar_or_new_edge(v5, v0, self.Ta, q5, q0)
                if e5[1]:
                    e5[0].c1 = cell
                else:
                    e5[0].c2 = cell

                cell.ListEdges=[e0,e1,e2,e3,e4,e5]
                cell.ListVertices = [v0,v1,v2,v3,v4,v5]

                cell.update_geometric_info()
                cid = cid +1
        self.update_vertex_connections()

        #for v in self.ListVertex:
        #    dr = Point(random.uniform(-1.0,1.0), random.uniform(-1.0,1.0), 0)
        #    v.coord = v.coord +dr.X(0.2)
    def build_tissue_from_file(self, datafile):
        self.ListCell, self.ListVertex, self.ListEdge = [], [], []
        with open(datafile, 'r') as fp:
            data = json.load(fp)

        self.SystemSize = Quadrant(data["box_size"][0], data["box_size"][1])
        self.isPeriodic = data["is_periodic"]

        ver_coords = data["vertices"]
        for v_coord in ver_coords:
            v = self.similar_or_new_vertex(Point(v_coord[0], v_coord[1], 0), 1)

        self.Kc, self.A0c, self.Ta, self.Gc = 1.0, 1.0, data["lineTension"], 0

        for cid in range(data["n_cells"]):
            cell = Cell(cid, 1, self.Kc, self.A0c, self.Ta, self.Gc)
            self.ListCell.append(cell)
            cell.storing_order = cid

        tissue_edges = data["edges"]
        for te in tissue_edges:
            crosses, t_id, h_id, t_qx, t_qy, h_qx, h_qy, tension, c1_id, c2_id = te[0], te[1], te[2], te[3], te[4], te[5], te[6], te[7], te[8], te[9]
            e = self.similar_or_new_edge(self.ListVertex[t_id], self.ListVertex[h_id], tension, Quadrant(t_qx, t_qy), Quadrant(h_qx, h_qy))
            e[0].crossBdry = crosses
            e[0].c1 = self.ListCell[c1_id]
            e[0].c2 = self.ListCell[c2_id]

        cells_ListVertex_ids = data["cells_ListVertex_ids"]
        cells_ListEdge_ids = data["cells_ListEdge_ids"]
        cells_bdry_crossing = data["cells_bdry_crossing"]

        cid = 0
        for cell_verts, cell_edges, cell_crossing in zip(cells_ListVertex_ids, cells_ListEdge_ids, cells_bdry_crossing):
            cell = self.ListCell[cid]
            for iv in range(len(cell_verts)):
                cell.ListVertices.append(self.ListVertex[cell_verts[iv]])
            #cell.ListVertices = self.ListVertex[np.array(cell_verts)]
            cell.crossBdry = cell_crossing
            cell.ListEdges = []
            for E in cell_edges:
                cell.ListEdges.append((self.ListEdge[E[0]], E[1]))
            cell.update_geometric_info()
            cid = cid + 1

        self.update_vertex_connections()

#        eid = 0
#        for te in tissue_edges:
#            crosses, t_id, h_id, t_qx, t_qy, h_qx, h_qy, tension, c1_id, c2_id = te[0], te[1], te[2], te[3], te[4], te[5], te[6], te[7], te[8], te[9]
#            self.ListEdge[eid].c1 = self.ListCell[c1_id]
#            self.ListEdge[eid].c2 = self.ListCell[c2_id]
#            eid = eid+1

    def write_tissue_to_file(self, filename, write_d2W_dv2):
        tissue_dict = {}
        tissue_dict["is_periodic"]=self.isPeriodic
        tissue_dict["n_cells"] = len(self.ListCell)
        tissue_dict["n_vertices"] = len(self.ListVertex)
        tissue_dict["n_edges"] = len(self.ListEdge)
        tissue_dict["box_size"]=[self.SystemSize.x, self.SystemSize.y]
        tissue_dict["lineTension"] = self.Ta
        vertex_list=[]
        [vertex_list.append([v.coord.x, v.coord.y, v.coord.z]) for v in self.ListVertex]
        tissue_dict["vertices"]=vertex_list
        if write_d2W_dv2:
            d2W = [v.d2W_dh2 for v in self.ListVertex]
            tissue_dict["d2W_dv2"] = d2W
        dWdh = [v.Norm_dW_dh for v in self.ListVertex]
        tissue_dict["dW_dh"] = dWdh
        k_probe = [v.probing_k for v in self.ListVertex]
        tissue_dict["probing_k"] = k_probe
        area_probe = [abs(v.final_opening_area)+1e-6 for v in self.ListVertex]
        tissue_dict["inserting_area"] = area_probe
        edge_list = []
        [edge_list.append([e.crossBdry, e.tailVertex.storing_order, e.headVertex.storing_order, e.q_tail.x, e.q_tail.y, e.q_head.x, e.q_head.y, e.lineTension, e.c1.storing_order, e.c2.storing_order]) for e in self.ListEdge]
        tissue_dict["edges"] = edge_list

        cell_list = []
        c_area = []
        cells_ListVertex_ids = []
        cells_ListEdge_ids = []
        cells_bdry_crossing = []
        for c in self.ListCell:
            if c.crossBdry:
                cells_bdry_crossing.append(True)
            else:
                cells_bdry_crossing.append(False)
            #    continue
            c_area.append(c.area)
            #c_edges = np.zeros((len(c.ListEdges)+1, 2))
            c_edges = []
            #c_edges[0,:]=np.array([c.ListEdges[0][0].tailVertex.coord.x, c.ListEdges[0][0].tailVertex.coord.y])
            #c_edges.append([c.ListEdges[0][0].tailVertex.coord.x, c.ListEdges[0][0].tailVertex.coord.y])
            cell_verts_ids = []
            cell_edges_ids = []
            for v in c.ListVertices:
                cell_verts_ids.append(v.storing_order)
            eid=0
            for E in c.ListEdges:
                e, d = E[0], E[1]
                cell_edges_ids.append([e.storing_order, d])
                #e.get_dl()
                #c_edges[eid+1,:] = np.array([e.dl.x, e.dl.y])
                if d:
                    c_edges.append([e.tailVertex.coord.x, e.tailVertex.coord.y])
                else:
                    c_edges.append([e.headVertex.coord.x, e.headVertex.coord.y])
                eid = eid+ 1
            #print(c_edges)
            cell_list.append(c_edges)
            cells_ListVertex_ids.append(cell_verts_ids)
            cells_ListEdge_ids.append(cell_edges_ids)
        tissue_dict["cells"] = cell_list
        tissue_dict["cells_ListVertex_ids"] = cells_ListVertex_ids
        tissue_dict["cells_ListEdge_ids"] = cells_ListEdge_ids
        tissue_dict["cells_bdry_crossing"] = cells_bdry_crossing


        tissue_dict["cell_area"]= c_area
        if filename.endswith(".json"):
            with open(filename, 'w') as fp:
                json.dump(tissue_dict, fp)
        else:
            with open(filename+".json", 'w') as fp:
                json.dump(tissue_dict, fp)
