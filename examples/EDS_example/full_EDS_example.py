####### a handy function for dealing with xlrd's interpretation of excel coloumns


from numpy import loadtxt, array, zeros, ones, arctan2, mean


def Jess_data_file(fname):
	fid = open(fname, 'r')
	lines = fid.readlines()
	fid.close()

	data = []
	for line in lines:
		data_line = []
		sline = line.split(',')
		for a in sline[0:-1]:
			data_line.append(float(a))
		data.append(data_line)
	data = array(data)

	return data




subdir =  ''
from matplotlib.pylab import *
rcParams.update({'font.size': 12})



from ucs_chemmap.tools import *
use_separate_line_scan = True

use_fancy_triangle = True
nsegs = number_of_triangles = 40
ab_step_size = 1.0
#################################################################################
#################################################################################

##########  row 0   ##########
## a nice colorset for RGB parity
if True:
	Lp = 61.5
	color_points = [[0.0, 1.0, 0.0],[0.0, 0.0, 1.0], [1.0, 0.0, 0.0] ]
	angle_0 = -180
	contrast_boost = 1.2
	angle_range = 120
	center = (0,0)
	style = 1

if False:
	Lp = 83
	color_points = [[0.0, 1.0, 1.0],[1.0, 0.0, 1.0], [1.0, 1.0, 0.0] ]
	angle_0 = -178
	contrast_boost = 1.2
	angle_range = 120
	center = (0,0)
	style = 2

if False:
	Lp = 74
	color_points = [[0.0, 1.0, 1.0],[1.0, 0.0, 1.0], [1.0, 1.0, 0.0] ]
	angle_0 = -178
	contrast_boost = 1.4
	angle_range = 120
	center = (0,0)
	style = 3

if False:
	Lp = 75
	color_points = [[0.0, 1.0, 1.0],[1.0, 0.0, 1.0], [1.0, 1.0, 0.0] ]
	angle_0 = -110
	contrast_boost = 1.3
	angle_range = 60
	center = (0,0)
	style = 4

if False:
	Lp = 70
	color_points = [[0.0, 1.0, 1.0],[1.0, 0.0, 1.0], [1.0, 1.0, 0.0] ]
	angle_0 = -180
	contrast_boost = 1.3
	angle_range = 120
	center = (0,0)
	style = 5

if False: # off center not implemented fully
	Lp = 74
	color_points = [[0.0, 1.0, 1.0],[1.0, 0.0, 1.0], [1.0, 1.0, 0.0] ]
	angle_0 = -90
	radius = 40
	contrast_boost = 1.0
	angle_range = 30
	center = (-10,20)
	style = 6
	nsegs = 10
	ab_step_size = 1



################################################################################
################################################################################

if True:
	fig, axes = subplots(figsize = (89/25.4, 8), ncols = 2, nrows = 6, dpi = 128)

	ax_0_0 = axes[0, 0]
	ax_0_1 = axes[0, 1]
	ax_1_0 = axes[1, 0]
	ax_1_1 = axes[1, 1]
	ax_3_0 = axes[3, 0]
	ax_3_1 = axes[3, 1]
	ax_4_0 = axes[4, 0]
	ax_4_1 = axes[4, 1]
	ax_5_0 = axes[5, 0]
	ax_5_1 = axes[5, 1]


	axes[2, 1].set_axis_off()
	if use_separate_line_scan:
		axes[2, 0].set_axis_off()
		scan_fig = figure(figsize = (89/25.4*0.8, 2.0),constrained_layout=True, dpi = 128)
		scan_ax = gca()
	else:
		scan_ax  = axes[2, 0]
else:
	fig = figure(figsize = (89/25.4, 8),constrained_layout=True, dpi = 128)
	import matplotlib.gridspec as gridspec
	gs = gridspec.GridSpec(ncols= 2, nrows = 6, figure = fig)

	ax_0_0 = fig.add_subplot(gs[0, 0])
	ax_0_1 = fig.add_subplot(gs[0, 1])
	ax_1_0 = fig.add_subplot(gs[1, 0])
	ax_1_1 = fig.add_subplot(gs[1, 1])
	ax_3_0 = fig.add_subplot(gs[3, 0])
	ax_3_1 = fig.add_subplot(gs[3, 1])
	ax_4_0 = fig.add_subplot(gs[4, 0])
	ax_4_1 = fig.add_subplot(gs[4, 1])
	ax_5_0 = fig.add_subplot(gs[5, 0])
	ax_5_1 = fig.add_subplot(gs[5, 1])

	scan_ax  = fig.add_subplot(gs[2, 0:2])




angle_0, radius = maximize_triangle_radius_and_angle_0( Lp = Lp, center = center, angle_0 = angle_0, angle_range = angle_range, angle_step = 0.2)
print(Lp, angle_0, radius)



elements = ['Pt', 'Ca','P' ]

data = []
for el in elements:
	data.append( Jess_data_file(subdir + el+' At%.csv').T)




sRGB1_map, sRGB1_color_points, L_scale = uniform_colormap(data = data,
							ax = ax_0_0,
							Lp = Lp,
							angle_0 = angle_0,
							radius = radius, center = center,
							contrast_boost = contrast_boost)
print(L_scale)
L_scale = 1.0

ax_0_0.minorticks_on()
ax_0_0.plot([1,1], [2,2], zorder = -2)


ax_3_0.imshow( deuteranomaly_with_alpha(sRGB1_map) )
ax_3_0.minorticks_on()

if use_fancy_triangle == False:
	sRGB1_map = corner_map(	data = data,
			elements = elements,
			ax = axes[1,0],
			Lp = Lp,
			radius = radius,
			angle_0 = angle_0,
			L_scale = L_scale)
	axes[4,0].imshow( deuteranomaly_with_alpha(sRGB1_map),
		origin = 'lower', extent=[0,1,0,1], interpolation = 'bicubic')


else:

	if True: CIELAB_color_triangle(ax = ax_1_0 ,
				Lp = Lp,
				radius = radius, center = center,
				angle_0 = angle_0,
				#L_scale = L_scale,
				nsegs = nsegs )



	if True: CIELAB_color_triangle(ax = ax_4_0 ,
					Lp = Lp,
					radius = radius, center = center,
					angle_0 = angle_0,
					#L_scale = L_scale,
					nsegs = nsegs,
					use_deuteranomaly = True )


	#inset_fig = figure(figsize = (1.5,1.5), dpi = 150)
	inset_ax_0 = ax_5_0
	#amin, amax = -128, 128
	#bmin, bmax = -128, 128
	CIELAB_triangle_with_LAB_slice(ax = inset_ax_0,
					Lp = Lp,
					radius = radius, center = center,
					angle_0 = angle_0,
					ab_step_size = ab_step_size)
	inset_ax_1 = ax_5_1
	CIELAB_triangle_with_LAB_slice(ax = inset_ax_1,
					Lp = Lp,
					radius = radius, center = center,
					angle_0 = angle_0,
					ab_step_size = ab_step_size,
					use_deuteranomaly = True)
	#inset_fig.tight_layout(pad = 0.1)
	#inset_fig.savefig('jess_inset_for_paper.pdf')
	#inset_fig.savefig('jess_inset_for_paper.svg')
##########  row 1   ##########
sRGB1_map = sRGB1_colormap(   data = data, elements = elements, ax = ax_0_1 , color_points = color_points)
ax_0_1.minorticks_on()

ax_3_1.imshow( deuteranomaly_with_alpha(sRGB1_map) )
ax_3_1.minorticks_on()

if use_fancy_triangle == False:
	sRGB1_map = corner_sRGB1_map( data = data, elements = elements, ax = axes[3,1], color_points = color_points)
	axes[4,1].imshow( deuteranomaly_with_alpha(sRGB1_map),
		origin = 'lower', extent=[0,1,0,1], interpolation = 'bicubic')
else:


	sRGB1_color_triangle(ax = ax_1_1,
					color_points = color_points,
					nsegs = nsegs)

	sRGB1_color_triangle(ax = ax_4_1,
				color_points = color_points,
				nsegs = nsegs,
				use_deuteranomaly = True)
########## scan_profile ####

scans_to_average = 20


sum_intensity = zeros(data[0].shape[1])

# get the overall max signal
normed_data = []
for data_index in range(3):
	normed_data.append(data[data_index]/data[data_index].max())
signal_sum = normed_data[0]+normed_data[1]+normed_data[2]
signal_sum_max = signal_sum.max()
overall_norm = signal_sum/signal_sum_max # must be less than 3

##
for i in range(3):
	mean_data = mean(normed_data[i][-1-scans_to_average:-1],axis =0)
	#line_data = data[i][0]
	sum_intensity += mean_data
	scan_ax.plot(mean_data/signal_sum_max , label = elements[i], color = clip(sRGB1_color_points[i],0,1))
scan_ax.plot(sum_intensity/signal_sum_max , label = 'Total',color = 'k')
scan_ax.set_ylabel("Normalized Intensity")
scan_ax.set_ylim(0,0.65)
scan_ax.set_xlim(0,data[0].shape[1]-1)

from matplotlib.ticker import MultipleLocator
scan_ax.xaxis.set_major_locator(MultipleLocator(base=50))
scan_ax.xaxis.set_minor_locator(MultipleLocator(base=10))


scan_ax.minorticks_on()
if use_separate_line_scan:
	scan_ax.legend(loc= 'upper right')

	scan_fig.tight_layout(pad = 0.1)
	scan_fig.subplots_adjust(right = 0.99, top = 0.99)
	scan_fig.savefig('line_scan_for_paper_style_%i.pdf'%style ,transparent = True, dpi =1200)
	scan_fig.savefig('line_scan_for_paper_style_%i.svg'%style,transparent = True, dpi =1200)
else:
	scan_ax.legend(bbox_to_anchor=(1.3, 1.0))
###################################################
fig.tight_layout(pad = 0.1)
#fig.subplots_adjust(left = 0.11, right = 0.99, hspace = 0.29) # no overlap for tests
fig.subplots_adjust(left = 0.07, right = 0.99, top = 0.99, bottom = 0.06, wspace = 0.16,  hspace = 0.0) # wide enough for actual publication
fig.savefig('data_for_paper_style_%i.pdf'%style,transparent = True, dpi =1200)
fig.savefig('data_for_paper_style_%i.svg'%style,transparent = True, dpi =1200)


show()
