

def Jess_data_file(fname):
	fid = open(fname, 'r')
	lines = fid.readlines()
	fid.close()
	from numpy import array

	data = []
	for line in lines:
		data_line = []
		sline = line.split(',')
		for a in sline[0:-1]:
			data_line.append(float(a))
		data.append(data_line)
	data = array(data)

	return data



elements = ['Pt', 'Ca','P' ]


# maximum lightness L* of 61.5 gives the widest triangle
Lp = 61.5
# using the contrast boost allows you to have some points that go outside the
# displayble color range so that you can better see the majority
contrast_boost = 1.0
# First hue angle of triangle vertices
angle_0 = -180.0
## This is is the number of triangles along the side of the mixing diagram
number_of_triangles = 6
## when looking at the LAB slice, the colors are computed on a grid, you can control how fine that grid is.
ab_step_size = 8.0

###### these you'll likely never change
## Center of Triangle
center = (0,0)
## rotation optimization range, no need to ever go larger than 120 degrees due to symmetry
angle_range = 120.0


data_list = []
label_list = []
for el in elements:
	data = Jess_data_file( el+' At%.csv').T
	data_list.append( data)
	### the contrast_boost here is important so that the labels are correctly scaled
	label =  '%.1f%%'%(data.max()/contrast_boost) + ' '+ el
	label_list.append(label)

############################################################################
##################################################################

from numpy import loadtxt, array, zeros, ones, arctan2, mean
from ucs_chemmap.tools import *
from matplotlib.pyplot import *
from matplotlib.ticker import MultipleLocator
rcParams.update({'font.size': 12})


angle_0, radius = maximize_triangle_radius_and_angle_0( Lp = Lp, center = center, angle_0 = angle_0, angle_range = angle_range, angle_step = 0.2)
print('L*', Lp, 'First Optimal Hue Angle', angle_0, 'Optimal Triangle Radius', radius)





fig_EDS, ax_EDS = subplots(figsize = (3.2, 2.5), dpi = 440/3.0)

sRGB1_map, sRGB1_color_points, L_scale = uniform_colormap(data = data_list,
							ax = ax_EDS,
							Lp = Lp,
							angle_0 = angle_0,
							pixel_size = 1.0,
							radius = radius, center = center,
							contrast_boost = contrast_boost,
							use_deuteranomaly = False)
ax_EDS.minorticks_on()
ax_EDS.xaxis.set_major_locator(MultipleLocator(base=20))
ax_EDS.xaxis.set_minor_locator(MultipleLocator(base=10))
ax_EDS.yaxis.set_major_locator(MultipleLocator(base=20))
ax_EDS.yaxis.set_minor_locator(MultipleLocator(base=10))
fig_EDS.tight_layout(pad=0.1)
fig_EDS.savefig('EDS_map.pdf', transparent = True, dpi =1200)
fig_EDS.savefig('EDS_map.svg', transparent = True, dpi =1200)


####################################################




fig_Tri, ax_Tri = subplots(figsize = (3.2, 2.5), dpi = 440/3.0)
CIELAB_color_triangle(ax = ax_Tri,
				Lp = Lp,
				radius = radius, center = center,
				angle_0 = angle_0,
				nsegs = number_of_triangles,
				labels = label_list,
				font_options = {'fontsize': 8},
				use_deuteranomaly = False )
fig_Tri.tight_layout(pad=0.1)
fig_Tri.savefig('mixing_triangle.pdf', transparent = True, dpi =1200)
fig_Tri.savefig('mixing_triangle.svg', transparent = True, dpi =1200)



####################################################

fig_LAB_Slice, ax_LAB_Slice = subplots(figsize = (3.2, 2.5), dpi = 440/3.0)
CIELAB_triangle_with_LAB_slice(ax = ax_LAB_Slice,
				Lp = Lp,
				radius = radius, center = center,
				angle_0 = angle_0,
				ab_step_size = ab_step_size,
				use_deuteranomaly = False)
fig_LAB_Slice.tight_layout(pad=0.1)
fig_LAB_Slice.savefig('CIELab_slice_with_triangle.pdf', transparent = True, dpi =1200)
fig_LAB_Slice.savefig('CIELab_slice_with_triangle.svg', transparent = True, dpi =1200)








########## scan_profile ####

scans_to_average = 20


fig_Scan, ax_Scan = subplots(figsize = (3.2, 2.5), dpi = 440/3.0)


sum_intensity = zeros(data_list[0].shape[1])

# get the overall max signal
normed_data = []
for data_index in range(3):
	normed_data.append(data_list[data_index]/data_list[data_index].max())
signal_sum = normed_data[0]+normed_data[1]+normed_data[2]
signal_sum_max = signal_sum.max()
overall_norm = signal_sum/signal_sum_max # must be less than 3

##
for i in range(3):
	mean_data = mean(normed_data[i][0: scans_to_average],axis =0)
	#line_data = data[i][0]
	sum_intensity += mean_data
	ax_Scan.plot(mean_data/signal_sum_max , label = elements[i], color = clip(sRGB1_color_points[i],0,1))
ax_Scan.plot(sum_intensity/signal_sum_max , label = 'Total',color = 'k')
ax_Scan.set_ylabel("Normalized Intensity")
ax_Scan.set_ylim(0,0.65)
ax_Scan.set_xlim(0,data_list[0].shape[1]-1)


ax_Scan.xaxis.set_major_locator(MultipleLocator(base=20))
ax_Scan.xaxis.set_minor_locator(MultipleLocator(base=10))


ax_Scan.minorticks_on()

ax_Scan.legend(loc= 'upper right')

fig_Scan.tight_layout(pad = 0.1)
fig_Scan.subplots_adjust(right = 0.99, top = 0.99)
fig_Scan.savefig('line_scan.pdf' ,transparent = True, dpi =1200)
fig_Scan.savefig('line_scan.svg',transparent = True, dpi =1200)


###################################################


show()
