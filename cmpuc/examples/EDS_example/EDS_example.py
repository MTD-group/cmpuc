







# maximum lightness L* of 61.5 gives the widest triangle
L_plane = 61.5

# First hue angle of triangle vertices
angle_0 = -180.0

radius = 30



### this will raise the L_plane until the data is at the edge of displayble colors
auto_tune_L_plane_to_data = False


# using the contrast boost allows you to have some points that go outside the
# displayble color range so that you can better see the majority
contrast_boost = 1.0


# basic red-green color blindness simulation
use_deuteranomaly = False

###### these you'll likely never change
## Center of Triangle
center = (0,0)

elements = ['Pt', 'Ca','P' ]

def data_file(fname):
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

data_list = []
for el in elements:
	data = data_file( el+' At%.csv').T
	data_list.append( data)


def make_labels(norm=1.0):
	label_list = []
	for i, el in enumerate(elements):
		### the contrast_boost here is important so that the labels are correctly scaled
		label =  '%.1f%%'%(norm*data_list[i].max()/contrast_boost) + ' '+ el
		label_list.append(label)
	return label_list
############################################################################
##################################################################

from numpy import loadtxt, array, zeros, ones, arctan2, mean, clip
from cmpuc.core import *
from matplotlib.pyplot import *
from matplotlib.ticker import MultipleLocator
rcParams.update({'font.size': 12})



my_map =   uniform_chemical_map(
							L_plane = L_plane,
							angle_0 = angle_0,
							radius = radius, center = center)
my_map.maximize_triangle_radius_and_angle_0(verbose = True)


print('L*', my_map.L_plane, 'First Optimal Hue Angle', my_map.angle_0, 'Optimal Triangle Radius', my_map.radius)





fig_EDS, ax_EDS = subplots(figsize = (3.2, 2.5), dpi = 440/3.0)


sRGB1_map  = my_map(
				data = data_list,
				contrast_boost = contrast_boost,
				auto_tune_L_plane_to_data = auto_tune_L_plane_to_data,
				use_deuteranomaly = use_deuteranomaly)
ax_EDS.imshow(np.clip(sRGB1_map,0,1), origin = 'lower')
ax_EDS.minorticks_on()
ax_EDS.xaxis.set_major_locator(MultipleLocator(base=20))
ax_EDS.xaxis.set_minor_locator(MultipleLocator(base=10))
ax_EDS.yaxis.set_major_locator(MultipleLocator(base=20))
ax_EDS.yaxis.set_minor_locator(MultipleLocator(base=10))
fig_EDS.tight_layout(pad=0.1)
fig_EDS.savefig('EDS_map.pdf', transparent = True, dpi =1200)
fig_EDS.savefig('EDS_map.svg', transparent = True, dpi =1200)
fig_EDS.savefig('EDS_map.png', transparent = True, dpi =300)

####################################################


number_of_triangles = 32

if auto_tune_L_plane_to_data==False and contrast_boost <= 1.0:

	fig_Tri, ax_Tri = subplots(ncols = 1, figsize = (3.2, 2.5), dpi = 440/3.0)
	isoluminant_triangle(ax_Tri, my_map,
					nsegs = number_of_triangles,
					labels = make_labels(),
					font_options = {'fontsize': 8},
					use_deuteranomaly = use_deuteranomaly)


else:
	norms = [ 0.5, 0.8, 1.0]
	fig_Tri, ax_Tri = subplots(ncols = len(norms), figsize = (6.5, 2.5), dpi = 440/3.0)

	for i in range(len(norms)):
		isoluminant_triangle(ax_Tri[i], my_map, norm =norms[i],
						nsegs = int(number_of_triangles*norms[i]),
						labels = make_labels(norm=norms[i]),
						font_options = {'fontsize': 8},
						use_deuteranomaly = use_deuteranomaly,
						remove_undisplayable = True )
		ax_Tri[i].set_title('Norm = %.1f'%norms[i], fontsize = 8)


fig_Tri.tight_layout(pad=0.1)
fig_Tri.savefig('mixing_triangle.pdf', transparent = True, dpi =1200)
fig_Tri.savefig('mixing_triangle.svg', transparent = True, dpi =1200)
fig_Tri.savefig('mixing_triangle.png', transparent = True, dpi =300)


####################################################

if auto_tune_L_plane_to_data==False and contrast_boost <= 1.0:

	fig_LAB_Slice, ax_LAB_Slice = subplots(figsize = (3.2, 2.5), dpi = 440/3.0)
	triangle_on_isoluminant_slice(ax_LAB_Slice, my_map,
		ab_step_size = 1.0, use_deuteranomaly= use_deuteranomaly)

else:
	fig_LAB_Slice, ax_LAB_Slice  = subplots(ncols = len(norms), figsize = (6.5, 2.5), dpi = 440/3.0)

	for i in range(len(norms)):
		triangle_on_isoluminant_slice(ax_LAB_Slice[i], my_map, norm = norms[i],
		ab_step_size = 1.0, use_deuteranomaly= use_deuteranomaly)
		ax_LAB_Slice[i].set_title('Norm = %.1f'%norms[i],fontsize = 8)

	fig_LAB_Slice.suptitle(my_map.color_space, fontsize=10)

fig_LAB_Slice.tight_layout(pad=0.1)
fig_LAB_Slice.savefig('CIELab_slice_with_triangle.pdf', transparent = True, dpi =1200)
fig_LAB_Slice.savefig('CIELab_slice_with_triangle.svg', transparent = True, dpi =1200)
fig_LAB_Slice.savefig('CIELab_slice_with_triangle.png', transparent = True, dpi =300)


########## 3D shape 

from plot_3d import UCS_pyramid_3D
fig = plt.figure()
ax = fig.gca(projection='3d')
UCS_pyramid_3D(ax,  chemical_map = my_map, nsegs = number_of_triangles)



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

if auto_tune_L_plane_to_data==False and contrast_boost <= 1.0:
	norm = 1.0
else:
	norm = norms[1]

sRGB1_color_points = my_map.get_sRGB1_color_points(norm)
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
fig_Scan.savefig('line_scan.png',transparent = True, dpi =300)

###################################################


show()
