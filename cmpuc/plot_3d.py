from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import numpy as np

from cmpuc.core import uniform_chemical_map, color_triangle_mesh, default_triangle_style,  cvd_space

from colorspacious import cspace_convert

#default_color_space = "CIELab"
#cvd_space = {"name": "sRGB1+CVD", "cvd_type":"deuteranomaly", "severity": 50}



def UCS_pyramid_3D(ax,
                chemical_map,
                nsegs = 20,
                labels = ['Part 1','Part 2','Part 3' ],
                order = [0,1,2],
                norm = 1.0,
                font_options = {},
                label_offset = 1.0,
                use_deuteranomaly = False, verbose = False,
                undisplayable_action = 'remove',
                color_for_undisplayable = (1.0, 1.0, 1.0, 0.0),
                triangle_style = default_triangle_style,
                alpha = 1.0 ):

    assert undisplayable_action in ['remove', 'clip', 'replace']
    #from matplotlib.patches import Polygon

    def fix_vert_order(verts):
        order = [1,2,0]
        verts_out  = np.zeros((3,3))
        for i in range(3):
            verts_out[i] = verts[i,order]
        return verts_out


    vert_list, color_list = color_triangle_mesh(
                color_points = chemical_map.get_color_points(),
                mode = '3D',
                nsegs = nsegs,
                order = order,
                norm = norm,
                verbose = verbose)

    for i in range(3):
        cp = chemical_map.get_color_points()
        cp[i,:] = 0.0
        vl, cl = color_triangle_mesh(
                color_points = cp,
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
    for vertices, color in zip(vert_list, color_list):
        sRGB1_color = cspace_convert(color, chemical_map.color_space, "sRGB1")

        if use_deuteranomaly: sRGB1_color = cspace_convert(sRGB1_color, cvd_space, "sRGB1")
        if  np.any(sRGB1_color<0,0) or np.any(1.0<sRGB1_color):
            if undisplayable_action != 'remove':
                if undisplayable_action == 'replace':
                    sRGB1_color = color_for_undisplayable
                else:
                    sRGB1_color = np.clip(sRGB1_color,0,1)
                vert_list_out.append(fix_vert_order(vertices))
                color_list_out.append(sRGB1_color)

        else:
#            ax.add_patch(  Polygon(vertices, color =np.clip(sRGB1_color,0,1), **triangle_style))
            vert_list_out.append(fix_vert_order(vertices))
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
        cp = chemical_map.get_color_points()
        for i in range(3):
            ax.text(cp[i,1], cp[i,2], cp[i,0]+ label_offset,labels[i] , ha = 'center', va = 'bottom', **font_options)
        #ax.text(   0,             -label_offset, labels[0], ha = 'center', va = 'top',    **font_options)
        #ax.text(  width,         -label_offset, labels[1], ha = 'center', va = 'top',    **font_options)
        #ax.text( width/2.0, height+label_offset, labels[2], ha = 'center', va = 'bottom', **font_options)

    ax.set_xlabel('a*')
    ax.set_ylabel('b*')
    ax.set_zlabel('L*')

    ax.set_xlim(-chemical_map.radius, chemical_map.radius)
    ax.set_ylim(-chemical_map.radius, chemical_map.radius)
    ax.set_zlim(0, chemical_map.L_plane+ label_offset)


    #ax.set_aspect('equal','box')
    #ax.set_xlim(0,1)
    #ax.set_ylim(0,sin(60*pi/180))
    #ax.set_axis_off()

#############################

if __name__ == "__main__":
    L_plane = 61.5
    # First hue angle of triangle vertices
    angle_0 = -86.0
    angle_0 = 0.0
    optimize_angle = False
    ## a radius for the triangle
    radius = 30


    nsegs = 5
    my_map =   uniform_chemical_map(
                                L_plane = L_plane,
                                angle_0 = angle_0,
                                radius = radius, )

    #my_map.maximize_triangle_radius_and_angle_0()
    my_map.maximize_triangle_radius()
    #################
    fig = plt.figure()
    ax = fig.gca(projection='3d')


    UCS_pyramid_3D(ax,  chemical_map = my_map, nsegs = nsegs)




    plt.show()
