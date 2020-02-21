from cmpuc import uniform_chemical_map,  isoluminant_triangle, triangle_on_isoluminant_slice
from matplotlib import pyplot as plt
import numpy as np

# maximum lightness L* of 61.5 gives the widest triangle
L_plane = 61.5
# First hue angle of triangle vertices
angle_0 = 60.0
## a radius for the triangle
radius = 30
###### these you'll likely never change
## Center of Triangle
center = (0,0)

my_map =   uniform_chemical_map(
							L_plane = L_plane,
							angle_0 = angle_0,
							radius = radius, center = center)
my_map.maximize_triangle_radius_and_angle_0(verbose = True)


### this will raise the L_plane until the data is at the edge of displayble colors
auto_tune_L_plane_to_data = False
# using the contrast boost allows you to have some points that go outside the
# displayble color range so that you can better see the majority
contrast_boost = 1.0

from cmpuc.analytic_data import test_compositions
sRGB1_map  = my_map(
				data = test_compositions,
				contrast_boost = contrast_boost,
				auto_tune_L_plane_to_data = auto_tune_L_plane_to_data,
				use_deuteranomaly = False)


################## now we can make plot and  accessory plots

fig, ax = plt.subplots()
ax.imshow(np.clip(sRGB1_map,0,1), origin = 'lower')


## when looking at the LAB slice, the colors are computed on a grid, you can control how fine that grid is.
fig, ax = plt.subplots()
triangle_on_isoluminant_slice(ax, my_map, ab_step_size = 1.0)


## This is is the number of triangles along the side of the mixing diagram
number_of_triangles = 6
fig, ax = plt.subplots()
isoluminant_triangle(ax, my_map,
        nsegs = number_of_triangles,
        labels = ['Part 1','Part 2','Part 3' ])

plt.show()
