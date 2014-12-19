#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       ortho2.py
#           stdin: {Objects}
#       
#       Script perform orthogonalization of selected ways, algorithm is different from JOSM's one
#       If some nodes are selected, they are not moved
#       If 2-node way having tag fixme=ortho2_py_base_vector is selected it will be used for base vector for orthogonalization
#       If two nodes are selected, they will be used for base vector (just like in JOSM)
#       Or base vector will be autocalculated from selected geometry
#       Note: nodes of ortho2_py_base_vector could be moved in the result, if connected to other way undergoing orthogonalization
#       
#       Copyright 2013-2014 OverQuantum
#       2013-03-13/19 main work on path_ortho, till RC3
#       2013-03-23 path_ortho project finalized, see https://github.com/OverQuantum/path_ortho
#       2014-12-13/16 converting code to python
#       2014-12-16 first working ortho on python, several objects support, no collapse
#       2014-12-17 added fixed nodes and base vector from them
#       2014-12-18 added support of ortho2_py_base_vector
#       2014-12-19 edged with both fixed nodes excluded from base vector and dirgroup calculations
#       
#       CONSIDER: introduce distance threshold for not-moved nodes
#       CONSIDER: collapsing and cleaning
#       
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

import sys
import math
import projections
import copy
from OsmData import OsmData, Map, LON, LAT, ACTION, MODIFY, REF, NODES, TAG

if sys.version_info[0] < 3:
    reload(sys)
    sys.setdefaultencoding("utf-8")                    # a hack to support UTF-8    

def main():
    data1 = OsmData()  # Dependencies
    data2 = OsmData()  # Selected objects
    data1.read(sys.stdin)
    data2.read(sys.stdin)
    
    if len(data2.ways)==0:
        data1.addcomment("No geometry for orthogonalization")
    else:
        path_ortho(data1, data2)
    
    #data.addcomment("Done.")
    data1.write(sys.stdout) #only nodes in dependencies could be modified
    return 0

def path_ortho(nodedata, selected):
    edges = []                              # all edges of all ways
    nodes = dict()                          # all used nodes

    fixnodes = len(selected.nodes)
    vbase1 = [0,0]
    
    #load geometry into our arrays with additional data
    n=0
    num = 0
    anchor = [0,0]                          # base point, all coords will be relative to it, to remove loss of float accuracy due to big absolute value
    basevectorflag = 0
    for wayid in selected.ways.keys():
        prevnode = None
        prevfix = 0
        if len(selected.ways[wayid][REF])==2 and basevectorflag==0:
            if 'fixme' in selected.ways[wayid][TAG]:
                if selected.ways[wayid][TAG]['fixme']=='ortho2_py_base_vector':
                    basevectorflag = 1      # manual base vector found
        for nodeid in selected.ways[wayid][REF]:    
            node = nodedata.nodes[nodeid]
            proj = projections.from4326((node[LON],node[LAT]), "EPSG:3857")
            if num==0:
                anchor = copy.deepcopy(proj)   #first found node
            proj[0] -= anchor[0]
            proj[1] -= anchor[1]
            if basevectorflag==1:
                vbase1[0] = proj[0]
                vbase1[1] = proj[1]
                basevectorflag=2
                continue
            elif basevectorflag==2:
                vbase1[0] -= proj[0]   #second node of manual base vector
                vbase1[1] -= proj[1]
                basevectorflag=3
                continue

            if nodeid in selected.nodes:
                fix = 1
            else:
                n+=1
                fix = 0

            if not nodeid in nodes:
                nodes[nodeid] = [proj,copy.deepcopy(proj), [None,None], [None, None], fix] #prev coords, new coords, C-params, dirgroups, fixed flag

            if prevnode!= None:
                if fix==0 or prevfix==0: #skip edges with both fixed nodes
                    edges.append([prevnode,nodeid,0])   #edge members: node1, node2, direction
                
            num+=1
            prevnode = nodeid
            prevfix = fix

    if n==0:
        nodedata.addcomment("All nodes fixed, nothing to move")
        return 0

    if basevectorflag==3:
        #use manual base vector
        vbase1 = normalize(vbase1)
    elif fixnodes==2:
        #use two selected nodes as a base vector
        num = 0
        for nodeid in selected.nodes.keys():
            if nodeid in nodes:
                proj = nodes[nodeid][0] #already calculated
            else:
                node = selected.nodes[nodeid]
                proj = projections.from4326((node[LON],node[LAT]), "EPSG:3857")
                proj[0] -= anchor[0]
                proj[1] -= anchor[1]
            if num==0:
                vbase1[0] = proj[0]
                vbase1[1] = proj[1]
            elif num==1: #just in case num became >1
                vbase1[0] -= proj[0]
                vbase1[1] -= proj[1]
            num +=1
        vbase1 = normalize(vbase1)
    else:
        #calculate base vector from geometry

        vbase2 = [0,0]
        
        #calculate two versions of base direction vector
        #1st on vectors in range (0;90) degrees
        #2nd on vectors in range (-45;45) degrees    
        for edge in edges:
            n1 = edge[0]
            n2 = edge[1]
            x1 = nodes[n2][0][0]-nodes[n1][0][0]
            y1 = nodes[n2][0][1]-nodes[n1][0][1]
            d = x1*x1 + y1*y1
            x1 *= d  #increase length for increasing effect of longer vectors
            y1 *= d
            if x1<0:
                x1=-x1
                y1=-y1
            #(x1;y1) now in range (-90;90) degr
            
            if y1<0:
                x2=-y1
                y2=x1
            else:
                x2=x1
                y2=y1
            #(x2,y2) now in range (0;90) degr

            vbase1[0]+=x2 #data goes into 1st base vector
            vbase1[1]+=y2
            
            if x1>abs(y1):
                x4=x1
                y4=y1  #in range (-45;45) degr
            else:
                if y1<0:
                    x4=-y1
                    y4=x1 #was in range (-90;-45) degr
                else:
                    x4=y1
                    y4=-x1 #was in range (45;90) degr
            #(x4,y4) now in range (-45;45) degr

            vbase2[0]+=x4 #data goes into 2nd base vector
            vbase2[1]+=y4

        #normalize both base vectors
        vbase1 = normalize(vbase1)
        vbase2 = normalize(vbase2)

        #calculate for both base vector square error
        sumv2=0
        sumv4=0
        for edge in edges:
            #path vector from one node to another
            n1 = edge[0]
            n2 = edge[1]
            x1 = nodes[n2][0][0]-nodes[n1][0][0]
            y1 = nodes[n2][0][1]-nodes[n1][0][1]
        
            #for 1st base vector
            v2 = x1 * vbase1[0] + y1 * vbase1[1] #length of path vector along base vector
            v4 = x1 * vbase1[1] - y1 * vbase1[0] #length of path vector perpendicular to base vector
            v2*=v2 #square
            v4*=v4
            if v2>v4:
                sumv2+=v4 #path vector is along base vector, square error is defined by perpendicular
            else:
                sumv2+=v2 #path vector is perpendicular to base vector, square error is defined by along length
            
            #for 2nd base vector
            v2 = x1 * vbase2[0] + y1 * vbase2[1] #length of path vector along base vector
            v4 = x1 * vbase2[1] - y1 * vbase2[0] #length of path vector perpendicular to base vector
            v2*=v2 #square
            v4*=v4
            if  v2>v4:
                sumv4+=v4 #path vector is along base vector, square error is defined by perpendicular
            else:
                sumv4+=v2 #path vector is perpendicular to base vector, square error is defined by along length
        
        if sumv2>sumv4:    #square error of 1st base vector is larger, so we will use 2nd vector as base
            vbase1=vbase2
        
        #for now on vbase1 is a base vector
        
    #set directions
    for edge in edges:
    #path vector from one node to another
        n1 = edge[0]
        n2 = edge[1]
        x1 = nodes[n2][0][0]-nodes[n1][0][0]
        y1 = nodes[n2][0][1]-nodes[n1][0][1]
        v2 = abs(x1*vbase1[0]+y1*vbase1[1]) #length of path vector along base vector
        v4 = abs(x1*vbase1[1]-y1*vbase1[0]) #length of path vector perpendicular to base vector
        if v2>v4:
            edge[2]=0 #path vector is along base vector
        else:
            edge[2]=1 #path vector is perpendicular to base vector

    #set dirgroups
    #each dirgroup contains all edges having same direction (0/1, see above) and connected by nodes
    #all nodes from one dirgroup are projected to single line during orthogonalization
    grouplist=[]   #list of dirgroups
    for edge in edges:
        n1 = edge[0]
        n2 = edge[1]
        dir = edge[2]
        if nodes[n1][3][dir]==None:
            if nodes[n2][3][dir]==None:
                #new group
                dirgr1 = [dir,n1,n2]         #group members: dirction, node1, node2, node3....
                nodes[n1][3][dir] = dirgr1
                nodes[n2][3][dir] = dirgr1
                grouplist.append(dirgr1)
            else:
                #add n1 to dirgroup
                dirgr1 = nodes[n2][3][dir]
                nodes[n1][3][dir]=dirgr1
                dirgr1.append(n1)
        else:
            if nodes[n2][3][dir]==None:
                #add n2 to dirgroup
                dirgr1 = nodes[n1][3][dir]
                nodes[n2][3][dir]=dirgr1
                dirgr1.append(n2)
            else:
                dirgr1 = nodes[n1][3][dir]
                dirgr2 = nodes[n2][3][dir]
                if dirgr1!=dirgr2:       
                    #combining different groups
                    len2 = len(dirgr2)
                    for i in range(1, len2):
                        nodeid=dirgr2[i]
                        nodes[nodeid][3][dir] = dirgr1  #update all nodes to dirgr1
                        dirgr1.append(nodeid)
                    for i in range(0,len(grouplist)-1):   #delete dirgr2
                        if grouplist[i]==dirgr2:
                            grouplist.pop(i)
                            break
                    del dirgr2 #to free memory
    
    #calculate nodes c1/c2 parameters
    #C-param defines line equation A*x+B*y = C, where (A,B) - vector perpendicular to line
    #(A,B) is our base vector for dir=1
    for dirgroup in grouplist:
        c = 0.0
        n = 0
        nfix = 0
        len1 = len(dirgroup)
        dir = dirgroup[0]
        cfix = 0.0
        for i in range(1, len1):
            nodeid = dirgroup[i]
            n+=1
            if dir==0:
                c0 = nodes[nodeid][0][0]*vbase1[1]-nodes[nodeid][0][1]*vbase1[0] #perpendicular to base vector
            else:
                c0 = nodes[nodeid][0][0]*vbase1[0]+nodes[nodeid][0][1]*vbase1[1] #along base vector
            c += c0
            if nodes[nodeid][4] == 1: #fixed node
                cfix += c0
                nfix += 1
        if n==0:
            c = None    #protection from division by zero
        else:
            c/=n         #average c-param
        if nfix > 0:
            c = cfix/nfix  #fixed nodes in group - use only their c-param average
        for i in range(1, len1):
            nodeid = dirgroup[i]
            if nodes[nodeid][4] == 1:
                nodes[nodeid][2][dir] = None  #fixed nodes are not moved
            else:
                nodes[nodeid][2][dir] = c
    
    #calculate new nodes positions from c1/c2 and update geometry
    #A*x+B*y=C1 - ortho-direction
    #B*x-A*y=C2 - along-direction
    #x=(A*C1+B*C2)/(A*A+B*B)
    #y=(-A*C2+B*C1)/(A*A+B*B)
    #where (A*A+B*B) = 1
    n=0
    for nodeid in nodes.keys():
        if nodes[nodeid][4] == 0: #not fixed
            c1 = nodes[nodeid][2][0]
            c2 = nodes[nodeid][2][1]
            
            #calc missing cl/c2 - for ending nodes, intermediates in chains and in case of 
            if c1==None:
                c1 = nodes[nodeid][0][0]*vbase1[1]-nodes[nodeid][0][1]*vbase1[0] #perpendicular to base vector
            if c2==None:
                c2 = nodes[nodeid][0][0]*vbase1[0]+nodes[nodeid][0][1]*vbase1[1] #along base vector
                
            #calc new positions
            nodes[nodeid][1][0]=vbase1[0]*c2+vbase1[1]*c1   #x
            nodes[nodeid][1][1]=-vbase1[0]*c1+vbase1[1]*c2  #y
    
            #update geometry back
            node = nodedata.nodes[nodeid]
            node_move = nodes[nodeid]
            if (node_move[1][0] != node_move[0][0]) or (node_move[1][1] != node_move[0][1]):
                node_lon,node_lat = projections.to4326((node_move[1][0]+anchor[0],node_move[1][1]+anchor[1]), "EPSG:3857")
                if (node[LON]!=node_lon) or (node[LAT]!=node_lat):
                    #both checks are mainly useless, as small changes present each time
                    node[LON]=node_lon
                    node[LAT]=node_lat
                    node[ACTION]=MODIFY    
                    n+=1
    
    nodedata.addcomment("Nodes moved: "+str(n)+", fixed: "+str(fixnodes))

def normalize(vector):
    l = 1.0/math.sqrt(vector[0]*vector[0] + vector[1]*vector[1])
    return (vector[0] * l, vector[1] * l)    
    
if __name__ == '__main__':
    main()
