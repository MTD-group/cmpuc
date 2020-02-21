# Composition Mapping within Perceptually Uniform Colorspaces (cmpuc) #
cmpuc is a small tool for performing perceptually uniform color mappings of three component composition fields. A common application of is EDS mapping of elements where traditionally red, green, and blue are used. cmpuc is designed to integrate into the Matplotlib/NumPy ecosystem and relies on [colorspacious](https://colorspacious.readthedocs.io/en/latest/) for colorspace transformations. The code is still at the prototype stage but should work for most cases. 

## Why? ##
For chemical mapping, the red, green, and blue primaries are typically used for three species. This method has serious interpretability mostly stemming from the large brightness differences between the sRGB primaries.

## Citing Cmpuc ##
Our paper exploring the concepts, the improvements, and a case study are [here](Will_be_on_ArXiv_soon).

## Usage ##
The basic idea is to map composition onto an inverted triangular pyramid in a uniform color space. The wider the pyramid base, the more chemical constrast there is. Likewise, the taller the pyramid, the more consentration contrast there is. In this example (*simple_example.py*), we put the base triangle plane at L = 61.5, and use a feature of the code to rotate and stretch the triange to be as wide as possibe while still fitting in the sRGB color space. From here, we plot some test data, show the base triange in the unform color space, and make mixing diagram to aid in reading the color mapping. For comparison, we also do the same with sRGB primaries.

```python
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
image  = my_map(
				data = test_compositions,
				contrast_boost = contrast_boost,
				auto_tune_L_plane_to_data = auto_tune_L_plane_to_data,
				use_deuteranomaly = False)

################## now we can make plot and  accessory plots

fig, ax = plt.subplots()
ax.imshow(np.clip(image,0,1), origin = 'lower')
fig.savefig('cmpuc_mapping.png',transparent = True, dpi = 150 )

## when looking at the LAB slice, the colors are computed on a grid, you can control how fine that grid is.
fig, ax = plt.subplots()
triangle_on_isoluminant_slice(ax, my_map, ab_step_size = 1.0)
fig.savefig('isoslice.png',transparent = True, dpi = 150 )

## This is is the number of triangles along the side of the mixing diagram
number_of_triangles = 20
fig, ax = plt.subplots()
isoluminant_triangle(ax, my_map,
        nsegs = number_of_triangles)
fig.savefig('isoluminant_triangle.png',transparent = True, dpi = 150 )

#### compare to RGB method
from cmpuc import sRGB1_colormap, sRGB1_triangle
color_points = [(1,0,0),(0,1,0),(0,0,1)]
fig, ax = plt.subplots()
image_sRGB1 = sRGB1_colormap(data = test_compositions, color_points = color_points)
ax.imshow(np.clip(image_sRGB1,0,1), origin = 'lower')
fig.savefig('sRGB1_mapping.png',transparent = True, dpi = 150 )

#####
fig, ax = plt.subplots()
sRGB1_triangle(ax, color_points = color_points, order = [0,2,1], #sRGB1 has opposite handednes from lab
        nsegs = number_of_triangles)
fig.savefig('sRGB1_triangle.png',transparent = True, dpi = 150 )

plt.show()
```
Here are the results of this example script:

![cmpuc_mapping](/cmpuc/examples/cmpuc_mapping.png) ![sRGB1_mapping](/cmpuc/examples/sRGB1_mapping.png)

1. The perceptually uniform chemical mapping.

![isoslice](/cmpuc/examples/isoslice.png)

2. The pyramid base plane in an isoluminant slice of the perceptually uniform color space. 

![isoluminant_triangle](/cmpuc/examples/isoluminant_triangle.png)

3. The chemical mixing triangle for understandin the colormapping, lower signal shrinks the triangle towards the black point.




