from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import numpy as np

from colorspacious import cspace_convert

from cmpuc import (color_triangle_mesh, default_triangle_style, cvd_space, default_color_space)






def sRGB_mixing_in_UCS_3D(ax,
                color_points,
                nsegs = 20,
                labels = ['Part 1','Part 2','Part 3' ],
                order = [0,1,2],
                norm = 1.0,
                font_options = {},
                label_offset = 1.0,
                use_deuteranomaly = False, verbose = False,
                undisplayable_action = 'clip',
                color_for_undisplayable = (1.0, 1.0, 1.0, 0.0),
                triangle_style = default_triangle_style,
                color_space = default_color_space,
                alpha = 1.0 ):

    '''color_points are in sRGB1'''
    assert undisplayable_action in ['remove', 'clip', 'replace']
    #from matplotlib.patches import Polygon

    def fix_vert_order(verts):
        order = [1,2,0]
        verts_out  = np.zeros((3,3))
        for i in range(3):
            verts_out[i] = verts[i,order]
        return verts_out




    cp = np.asarray(color_points)
    z = np.zeros(3)

    super_vert_sets = [
        # primaries
        [z,     cp[0],       cp[1]],
        [cp[0], cp[1], cp[0]+cp[1]],
        [z,     cp[1],       cp[2]],
        [cp[1], cp[2], cp[1]+cp[2]],
        [z,     cp[2],       cp[0]],
        [cp[2], cp[0], cp[2]+cp[0]],
        # secondaries
        [cp[0],             cp[0]+cp[1], cp[0]+cp[2]],
        [cp[0]+cp[1]+cp[2], cp[0]+cp[1], cp[0]+cp[2]],
        [cp[1],             cp[1]+cp[0], cp[1]+cp[2]],
        [cp[0]+cp[1]+cp[2], cp[1]+cp[0], cp[1]+cp[2]],
        [cp[2],             cp[2]+cp[0], cp[2]+cp[1]],
        [cp[0]+cp[1]+cp[2], cp[2]+cp[0], cp[2]+cp[1]]
        ]

    vert_list = []
    color_list = []

    for sv in super_vert_sets:
        vl, cl = color_triangle_mesh(
                color_points = sv,
                mode = '3D',
                nsegs = nsegs,
                order = order,
                norm = norm,
                verbose = verbose)
        vert_list = vert_list + vl
        color_list = color_list + cl
    #print(norm)

    #print(color_list[-1])
    ## now we filter them
    vert_list_out = []
    color_list_out = []
    for vertices, sRGB1_color in zip(vert_list, color_list):

        #print(vertices)
        lab_vertices = cspace_convert(vertices,  "sRGB1", color_space,)
        #cspace_convert(color, chemical_map.color_space, "sRGB1")

        if use_deuteranomaly: sRGB1_color = cspace_convert(sRGB1_color, cvd_space, "sRGB1")
        if  np.any(sRGB1_color<0,0) or np.any(1.0<sRGB1_color):
            if undisplayable_action != 'remove':
                if undisplayable_action == 'replace':
                    sRGB1_color = color_for_undisplayable
                else:
                    sRGB1_color = np.clip(sRGB1_color,0,1)
                vert_list_out.append(fix_vert_order(lab_vertices))
                color_list_out.append(sRGB1_color)

        else:
#            ax.add_patch(  Polygon(vertices, color =np.clip(sRGB1_color,0,1), **triangle_style))
            vert_list_out.append(fix_vert_order(lab_vertices))
            color_list_out.append(sRGB1_color)

    poly = Poly3DCollection(vert_list_out, facecolors=color_list_out, alpha=alpha, **triangle_style)
    ax.add_collection3d(poly, zdir='z')

    if False: # for some, reason the polygons keep overriding the lines.
        cp = my_map.get_color_points()*1.01
        cp = np.vstack([cp, (-0.5,0,0)])

        point_pairs = [(0,1), (1,2), (2,0), (3,0), (3,1), (3,2)]
        for pair in point_pairs:
            xs = [ cp[pair[0]][1],  cp[pair[1]][1] ]
            ys = [ cp[pair[0]][2],  cp[pair[1]][2] ]
            zs = [ cp[pair[0]][0],  cp[pair[1]][0] ]
            ax.plot(xs, ys, zs, color = 'k', linewidth = 2.0)


    if len(labels) ==3:
        #cp = chemical_map.get_color_points()
        lab_cp = cspace_convert(cp,  "sRGB1", color_space)
        for i in range(3):
            ax.text(lab_cp[i,1], lab_cp[i,2], lab_cp[i,0]+ label_offset,labels[i] , ha = 'center', va = 'bottom', **font_options)
        #ax.text(   0,             -label_offset, labels[0], ha = 'center', va = 'top',    **font_options)
        #ax.text(  width,         -label_offset, labels[1], ha = 'center', va = 'top',    **font_options)
        #ax.text( width/2.0, height+label_offset, labels[2], ha = 'center', va = 'bottom', **font_options)

    ax.set_xlabel('a*')
    ax.set_ylabel('b*')
    ax.set_zlabel('L*')

    vlo = np.asarray(vert_list_out)
    # they are reordered to be like xyz so a,b,L
    a_lims = (vlo[:,:,0].min(), vlo[:,:,0].max())
    b_lims = (vlo[:,:,1].min(), vlo[:,:,1].max())
    L_lims = (vlo[:,:,2].min(), vlo[:,:,2].max())

    ax.set_xlim(a_lims[0], a_lims[1])
    ax.set_ylim(b_lims[0], b_lims[1])
    ax.set_zlim(L_lims[0], L_lims[1]+ label_offset)




if __name__ == "__main__":

    #################
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    nsegs = 8

    color_points = [(1,0,0),(0,1,0),(0,0,1)]

    sRGB_mixing_in_UCS_3D(ax,  color_points, nsegs = nsegs)




    plt.show()
