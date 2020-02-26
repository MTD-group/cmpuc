# Composition Mapping within Perceptually Uniform Colorspaces (cmpuc) #
cmpuc is a tool for performing perceptually uniform color mappings of three component composition fields. A common application occurs in energy-dispersive x-ray spectroscopy (EDS) mapping of elements in multicomponent materials where traditionally red, green, and blue are used to represent different comprising elements. cmpuc provides better visualization fidelity by replacing the traditional primary colors

cmpuc integrates into the Matplotlib/NumPy ecosystem and relies on [colorspacious](https://colorspacious.readthedocs.io/en/latest/) for colorspace transformations. The code is still at the prototype stage but should work for most cases. By default, cmpuc uses the CIELAB color space; however, CAM02-UCS is implemented although untested. 

## Why? ##
For chemical mapping, the red, green, and blue primaries are typically used for three species. This method has serious interpretability mostly stemming from the large brightness differences between the sRGB primaries.

## Citing Cmpuc ##
Our paper exploring the concepts, the improvements, and a case study are [here](Will_be_on_ArXiv_soon).
[![DOI](https://zenodo.org/badge/181064496.svg)](https://zenodo.org/badge/latestdoi/181064496)

## Usage ##
The basic idea is to map composition onto an inverted triangular pyramid in a uniform color space (UCS). The wider the pyramid base, the more chemical constrast there is. Likewise, the taller the pyramid, the more concentration contrast there is. In this example (*simple_example.py*), we put the base triangle plane at L = 61.5, and use a feature of the code to rotate and stretch the triange to be as wide as possibe while still fitting in the sRGB color space. From here, we plot some test data, show the base triange in the unform color space, and construct a mixing diagram to aid in reading the color mapping. For comparison, we also do the same with sRGB primaries.

```python
from cmpuc import uniform_chemical_map,  isoluminant_triangle, triangle_on_isoluminant_slice
from matplotlib import pyplot as plt
import numpy as np

## maximum lightness L* of 61.5 gives the widest triangle in CIELab
L_plane = 61.5
# First hue angle of triangle vertices
angle_0 = 60.0
## a radius for the triangle
radius = 30
## Center of Triangle
center = (0,0)

my_map =   uniform_chemical_map(
               L_plane = L_plane,
               angle_0 = angle_0,
               radius = radius, center = center)
my_map.maximize_triangle_radius_and_angle_0(verbose = True)

```
Now we can import some test data. The ```my_map``` expects three 2D NumPy arrays in a list or array in the form of:  ```my_map([chemical_1_2D_array, chemical_2_2D_array, chemical_3_2D_array])``` and it returns an image in a NumPy array compatible with Matplotlib's imshow and imsave *i.e.* ```image[x_index,y_index,rgb_index]```. 
```python
### This will raise the L_plane until the data is at the edge of displayble colors
auto_tune_L_plane_to_data = False
# Using the contrast boost allows you to have some points that go outside the
# displayble color range so that you can better see the majority, use with great caution
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

## when looking at the LAB slice, the colors are computed on a grid,
## you can control how fine that grid is.
fig, ax = plt.subplots()
triangle_on_isoluminant_slice(ax, my_map, ab_step_size = 1.0)
fig.savefig('isoslice.png',transparent = True, dpi = 150 )

## This is is the number of triangles along the side of the mixing diagram
number_of_triangles = 20
fig, ax = plt.subplots()
isoluminant_triangle(ax, my_map,
        nsegs = number_of_triangles)
fig.savefig('isoluminant_triangle.png',transparent = True, dpi = 150 )

## compare to old RGB method
from cmpuc import sRGB1_colormap, sRGB1_triangle
color_points = [(1,0,0),(0,1,0),(0,0,1)]
fig, ax = plt.subplots()
image_sRGB1 = sRGB1_colormap(data = test_compositions, color_points = color_points)
ax.imshow(np.clip(image_sRGB1,0,1), origin = 'lower')
fig.savefig('sRGB1_mapping.png',transparent = True, dpi = 150 )

## The corresponding sRGB triangle
fig, ax = plt.subplots()
sRGB1_triangle(ax, color_points = color_points, order = [0,2,1], #sRGB1 has opposite handednes from lab
        nsegs = number_of_triangles)
fig.savefig('sRGB1_triangle.png',transparent = True, dpi = 150 )

plt.show()
```
Here are the results of this example script:

<img src="/cmpuc/examples/cmpuc_mapping.png" width="400" ><img src="/cmpuc/examples/sRGB1_mapping.png" width="400" >

1. (left) The perceptually uniform chemical mapping and (right) the traditional RGB primary mapping. 

<img src="/cmpuc/examples/isoluminant_triangle.png" width="400" ><img src="/cmpuc/examples/sRGB1_triangle.png" width="400" >

2. (left) The chemical mixing triangle for the cmpuc color mapping and (right) the traditional RGB chemical mixing triangle.

<img src="/cmpuc/examples/isoslice.png" width="500" >

3. The pyramid base plane is an isoluminant slice of the perceptually uniform color space. Lower signal shrinks the triangle towards the black point, *i.e.* the point of the pyramid.


## Using Real Data ##
In our EDS example (*EDS_example.py*) from our [paper](Will_be_on_ArXiv_soon), we use the same isoluminant triangle base as in the previus example. We now see how the chemical mapping performs with real data.

<img src="/cmpuc/examples/EDS_example/EDS_map.png" width="400" ><img src="/cmpuc/examples/EDS_example/mixing_triangle.png" width="400" >

<img src="/cmpuc/examples/EDS_example/line_scan.png" width="400" ><img src="/cmpuc/examples/EDS_example/CIELab_slice_with_triangle.png" width="400" >

(Top left) The chemical mapping of the EDS data with (top right) the chemical mixing triangle. (Bottom left) an average of the bottom 20 line scans from the chemical mapping. (Bottom right) The triangle in the perceptually uniform color space.

### Enhancing Signal Brightness ###
The mapping is quite dark since the signals have little overlap. We have implemented a feature to raise the pyramid base beyond the bounds of the sRGB color space while both maintaing the triangle radius and orientation and ensuring all of the mapped data is still in the bounds of the sRGB color space. By enabling this feature with ```auto_tune_L_plane_to_data=True```, there can be more brightness and signal constrast.

<img src="/cmpuc/examples/EDS_example/EDS_map_autotuned.png" width="400" >

Now the mixing triangle really must be shown as series of slices of the inverted pyramid like so:

<img src="/cmpuc/examples/EDS_example/mixing_triangle_autotuned.png" width="1000" >

In this case, we show slices at 0.5x, 0.8x, and 1.0x of the maximum lightness plane. The corresponding slices in the UCS:

<img src="/cmpuc/examples/EDS_example/CIELab_slice_with_triangle_autotuned.png" width="1000" >


