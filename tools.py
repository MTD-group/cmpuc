

default_corner_points = 256
from numpy import array, cos, pi, sin, zeros, ones, clip
from colorspacious import cspace_convert


def bounds_test(sRGB_map):
	nout = [0,0,0,0]
	for ix in range(sRGB_map.shape[0]):
		for iy in range(sRGB_map.shape[1]):
			state = 0
			for color in range(3):
				if sRGB_map[ix,iy,color] < 0.0 or sRGB_map[ix,iy,color] > 1.0:
					state+=1
			nout[state] += 1

	total = sRGB_map.shape[0]*sRGB_map.shape[1]
	print('Out of Bounds Triply: %i     (%.3f%%)'% (nout[3], nout[3]*100.0/total) )
	print('Out of Bounds Doubly: %i     (%.3f%%)'% (nout[2], nout[2]*100.0/total) )
	print('Out of Bounds Singly: %i     (%.3f%%)'% (nout[1], nout[1]*100.0/total) )
	print('Total Out of Bounds : %i/%i (%.3f%%)'% (sum(nout[1:]), total, sum(nout[1:])*100.0/total ) )
	return state



def uniform_colormap(data,
					ax,
					Lp,
					radius,
					angle_0,
					center = (0.0,0.0), pixel_size = 1.0, contrast_boost = 1.0, verbose = True, use_deuteranomaly = False ):

	color_points = []
	color_points.append( array([Lp, radius*cos(pi/180*(angle_0 + 0))  , radius*sin(pi/180*(angle_0 +  0))])  )
	color_points.append( array([Lp, radius*cos(pi/180*(angle_0+ 120)) , radius*sin(pi/180*(angle_0 + 120))])  )
	color_points.append( array([Lp, radius*cos(pi/180*(angle_0+ 240)) , radius*sin(pi/180*(angle_0 + 240))])  )
	cvd_space = {"name": "sRGB1+CVD", "cvd_type":"deuteranomaly", "severity": 50}

	for i in range(3):
		color_point = color_points[i]
		sRGB1_cp = cspace_convert(color_point, "CIELab", "sRGB1")
		if verbose: print('sRGB1 color point %i:' %i, "(%.3f %.3f %.3f)"% tuple(sRGB1_cp))
	#total = data[0] + data[1]+ data[2]

	normed_data = []
	for data_index in range(3):
		normed_data.append(data[data_index]/data[data_index].max() )
		#print(data[data_index].max())

	overall_norm = normed_data[0]+normed_data[1]+normed_data[2]
	overall_norm = overall_norm/overall_norm.max()

	data_shape = data[0].shape
	Lab_map =  zeros(list(data_shape)+[3])
	sRGB1_map = ones(list(data_shape)+[4])

	for i in range(data_shape[0]):
		for j in range(data_shape[1]):
			Lab_map[i,j,0] = Lp * overall_norm[i,j]
			for data_index in  range(3):
				# color points are really vectors, the normed_data is the projection along each one, and overall_norm is the overall scalling
				if True:
					Lab_map[i,j,1] += (color_points[data_index][1] * normed_data[data_index][i,j] + center[0])  * overall_norm[i,j]
					Lab_map[i,j,2] += (color_points[data_index][2] * normed_data[data_index][i,j] + center[1])  * overall_norm[i,j]
				else:
					Lab_map[i,j,1] += color_points[data_index][1] * normed_data[data_index][i,j] * overall_norm[i,j] + center[0]
					Lab_map[i,j,2] += color_points[data_index][2] * normed_data[data_index][i,j] * overall_norm[i,j] + center[1]


	max_L = Lab_map[:,:,0].max()
	if verbose: print ('Max L', Lab_map[:,:,0].max())

	### now we can scale it all to fit tighter
	# step 1, do a log search upwards to find an upper cutoff
	from numpy import any as np_any
	L_scale = 1.0
	Lab_test = L_scale* Lab_map
	sRGB1_test = cspace_convert(Lab_test, "CIELab", "sRGB1")
	while np_any(sRGB1_test>1.0) == False: # loop to find an upper bound
		L_scale = L_scale*1.2 # doubling is somewhat dumb, since we know we are on the right order of magnitude
		Lab_test[:,:,0] = L_scale * Lab_map[:,:,0]
		sRGB1_test = cspace_convert(Lab_test, "CIELab", "sRGB1")


	# this finds the perfect scalling of L for the color points
	cut_off = 1.0/1024 # half a 16 bit color depth unit, very accurate
	L_scale_max = L_scale
	L_scale_min = 0.0
	sRGB1_max = sRGB1_test.max()
	while ( 1.0 - sRGB1_max ) > cut_off or sRGB1_max > 1:
		L_scale = 0.5*(L_scale_max + L_scale_min)
		Lab_test[:,:,0] = L_scale * Lab_map[:,:,0]
		sRGB1_test = cspace_convert(Lab_test, "CIELab", "sRGB1")
		sRGB1_max = sRGB1_test.max()

		if sRGB1_max>1.0:
			L_scale_is_safe = False
			if verbose:print (L_scale_min, L_scale, L_scale_max, L_scale_is_safe)
			L_scale_max = L_scale

		else:
			L_scale_is_safe = True
			if verbose:print (L_scale_min, L_scale, L_scale_max, L_scale_is_safe)
			L_scale_min = L_scale

	if L_scale_is_safe==False:
		L_scale = L_scale_min

	# only use L_scale if it helps!
	if L_scale < 1.0:
		L_scale = 1.0

	Lab_test[:,:,0] = L_scale * Lab_map[:,:,0]
	sRGB1_test = cspace_convert(contrast_boost * Lab_test, "CIELab", "sRGB1")



	sRGB1_map[:,:,0:3] = sRGB1_test[:,:,0:3]
	print ('Max sRGB1', sRGB1_map[:,:,0:3].max())

	if verbose: bounds_test(sRGB1_map)
	if use_deuteranomaly:
		sRGB1_map = clip(cspace_convert(sRGB1_map, cvd_space, "sRGB1"),0,1)
	ax.imshow(clip(sRGB1_map,0,1), origin = 'lower',  extent=[0,sRGB1_map.shape[1]*pixel_size ,0,sRGB1_map.shape[0]*pixel_size])
	ax.set_axis_on()

	sRGB1_color_points = []
	center_array = array([0,center[0], center[1]])
	for i in range(3):
		sRGB1_color_points.append(cspace_convert(color_points[i]+center_array, "CIELab", "sRGB1"))

	if use_deuteranomaly:
		for i in range(3):
			sRGB1_color_points[i] = clip(cspace_convert(sRGB1_color_points[i], cvd_space, "sRGB1"),0,1)

	return  sRGB1_map, sRGB1_color_points, L_scale


def corner_map(	data,
				elements,
				ax,
				Lp,
				radius,
				angle_0,
				L_scale,
				npts = default_corner_points,
				bg_color = array([1, 1, 1, 0]) ):

	color_points = []
	color_points.append( array([Lp,  radius*cos(pi/180*(angle_0+  0))  ,   radius*sin(pi/180*(angle_0+   0))])  )
	color_points.append( array([Lp,  radius*cos(pi/180*(angle_0+ 120)) ,   radius*sin(pi/180*(angle_0+ 120))])  )
	color_points.append( array([Lp,  radius*cos(pi/180*(angle_0+ 240)) ,   radius*sin(pi/180*(angle_0+ 240))])  )

	from numpy import linspace
	values = linspace(0,1,npts)

	Lab_guide = zeros((npts,npts,3))
	fractions = zeros((npts,npts,3))
	sRGB1_map = ones((npts,npts,4))

	for i in range(npts):
		iv = values[i]
		for j in range(npts-i):
			jv = values[j]

			fraction = array( [(1-jv-iv), iv, jv ])
			#print(iv, jv, fraction)
			for data_index in  range(3):
				Lab_guide[j,i] += color_points[data_index] * fraction[data_index]* L_scale

	sRGB1_map[:,:,0:3] = cspace_convert(Lab_guide, "CIELab", "sRGB1")
	for i in range(npts):
		for j in range(npts-i,npts):
			sRGB1_map[j,i] = bg_color

	maxes = []
	for dat in data:
		maxes.append(dat.max())

	ax.imshow(clip(sRGB1_map,0,1), origin = 'lower', extent=[0,1,0,1], interpolation = 'bicubic')
	ax.text(0,0, elements[0]+'\n%0.1f%%'% maxes[0], ha = 'left')
	ax.text(1,0, elements[1]+'\n%0.1f%%'% maxes[1], ha = 'right')
	ax.text(0,1, elements[2]+'\n%0.1f%%'% maxes[2], ha = 'left', va = 'top')

	return clip(sRGB1_map,0,1)




def sRGB1_colormap(data,
					elements,
					ax,
					color_points):

	maxes = []
	for dat in data:
		maxes.append(dat.max())

	data_shape = data[0].shape
	sRGB1_map = zeros(list(data_shape)+[4])
	sRGB1_color_points = array(color_points)

	for i in range(data_shape[0]):
		for j in range(data_shape[1]):
			sRGB1_map[i,j,3] = 1.0
			for data_index in  range(3):
				sRGB1_map[i,j,0:3] += sRGB1_color_points[data_index] * data[data_index][i,j]/maxes[data_index]

	ax.imshow(clip(sRGB1_map,0,1))


	return clip(sRGB1_map,0,1)


def corner_sRGB1_map(
				data,
				elements,
				ax,
				color_points,
				npts = default_corner_points,
				bg_color = array([1, 1, 1, 0]) ):


	from numpy import linspace
	values = linspace(0,1,npts)

	sRGB1_map = ones((npts,npts,4))
	fractions = zeros((npts,npts,3))
	sRGB1_color_points = array(color_points)
	bg_color = array( [1,1,1,0])

	for i in range(npts):
		iv = values[i]
		for j in range(npts-i):
			jv = values[j]

			fraction = array( [(1-jv-iv), iv, jv ])
			#print(iv, jv, fraction)
			sRGB1_map[j,i] = array([0,0,0,1])
			for data_index in  range(3):
				sRGB1_map[j,i,0:3] += sRGB1_color_points[data_index] * fraction[data_index]

		for j in range(npts-i,npts):
			sRGB1_map[j,i] = bg_color

	maxes = []
	for dat in data:
		maxes.append(dat.max())
	ax.imshow(clip(sRGB1_map,0,1), origin = 'lower', extent=[0,1,0,1], interpolation = 'bicubic')
	ax.text(0,0, elements[0]+'\n%0.1f%%'% maxes[0], ha = 'left')
	ax.text(1,0, elements[1]+'\n%0.1f%%'% maxes[1], ha = 'right')
	ax.text(0,1, elements[2]+'\n%0.1f%%'% maxes[2], ha = 'left', va = 'top')

	return clip(sRGB1_map,0,1)


def deuteranomaly_with_alpha(sRGB1_map):
	cvd_space = {"name": "sRGB1+CVD", "cvd_type":"deuteranomaly", "severity": 50}
	sRGB1_map_d = zeros(sRGB1_map.shape)
	sRGB1_map_d[:,:,0:3] = clip(cspace_convert(sRGB1_map[:,:,0:3], cvd_space, "sRGB1"),0,1)
	sRGB1_map_d[:,:,3] = sRGB1_map[:,:,3]
	return sRGB1_map_d

# I'm not sure that 'antialiased' is needed to be off, but it's suggested online
# there needs to be *some* line width or the vectorized representations like pdf have stupid lines
# the joinstyle of bevel or round makes the lines not shoot outwards like miter
triangle_style = {'antialiased': False, 'linewidth' : 0.1, 'joinstyle': 'bevel' }

def sRGB1_color_triangle(ax,  color_points = [[1.0, 0.0, 1.0], [0.0, 1.0, 1.0], [1.0, 1.0, 0.0]] , nsegs = 40, use_deuteranomaly = False, verbose = False):
	from matplotlib.patches import Polygon

	patches = []
	colors = []
	from numpy import sin, tan, array, pi, cos, clip
	color_ps = array([color_points[0], color_points[2], color_points[1]])
	#color_ps = array(color_points)
	delta = (1/nsegs)
	deltay = sin(60*pi/180)/nsegs

	cvd_space = {"name": "sRGB1+CVD", "cvd_type":"deuteranomaly", "severity": 50}
	from colorspacious import cspace_convert

	from numpy.linalg import inv
	map_matrix = array([
						[0,   0,              1],
						[0.5, sin(60*pi/180), 1],
						[1.0, 0,              1]]).T


	imap = inv(map_matrix)

	for i in range(3):
		if verbose: print ('point',i)
		if verbose: print(map_matrix.T[i])
		fracs = imap.dot(map_matrix.T[i])
		if verbose: print ('fracs', fracs)
		color = color_ps[0]*fracs[0]+color_ps[1]*fracs[1]+color_ps[2]*fracs[2]
		if verbose: print('color',color)


	if verbose: print('upwards triangles')
	for iy in range(nsegs):
		for ix in range(nsegs-iy):
			vertices = [
						(delta*(    0.5*iy+ix), deltay *   iy),
						(delta*(0.5+0.5*iy+ix), deltay*(iy+1)),
						(delta*(1.0+0.5*iy+ix), deltay *   iy)]

			xc = delta*(0.5+0.5*iy+ix)
			yc = deltay * iy + 0.5*delta*tan(30*pi/180)
			fracs = imap.dot([xc,yc,1])
			#print (fracs)
			color = color_ps[0]*fracs[0]+color_ps[1]*fracs[1]+color_ps[2]*fracs[2]
			if use_deuteranomaly: color = cspace_convert(color, cvd_space, "sRGB1")
			ax.add_patch(  Polygon(vertices, color =clip(color,0,1), **triangle_style))


	if verbose: print('downwards triangles')
	for iy in range(nsegs-1):
		for ix in range(nsegs-iy-1):
			vertices = [
						(delta*(0.5+0.5*iy+ix), deltay *(iy+1)),
						(delta*(1.0+0.5*iy+ix), deltay * iy   ),
						(delta*(1.5+0.5*iy+ix), deltay *(iy+1))]
			xc = delta*(1.0+0.5*iy+ix)
			yc = deltay *(iy+1) - 0.5*delta*tan(30*pi/180)
			fracs = imap.dot([xc,yc,1])
			#print (fracs)
			color = color_ps[0]*fracs[0]+color_ps[1]*fracs[1]+color_ps[2]*fracs[2]
			if use_deuteranomaly: color = clip(cspace_convert(color, cvd_space, "sRGB1"),0,1)
			ax.add_patch(  Polygon(vertices, color =clip(color,0,1), **triangle_style))

	#ax.set_aspect('equal', 'box')


	ax.set_xlim(0,1)
	ax.set_ylim(0,sin(60*pi/180))
	ax.set_axis_off()
	ax.set_aspect('equal')


def CIELAB_color_triangle(ax,
				Lp,
				radius,
				angle_0,
				center = (0.0, 0.0),
				nsegs = 40,
				L_scale = 1.0,
				labels = [],
				font_options = {},
				label_offset = 0.02,
				use_deuteranomaly = False,verbose = False):
	from matplotlib.patches import Polygon
	from colorspacious import cspace_convert

	cvd_space = {"name": "sRGB1+CVD", "cvd_type":"deuteranomaly", "severity": 50}
	from numpy import sin, tan, array, clip, pi, cos
	color_points = []
	color_points.append( array([Lp,  radius*cos(pi/180*(angle_0+  0))  ,   radius*sin(pi/180*(angle_0+   0))])  )
	color_points.append( array([Lp,  radius*cos(pi/180*(angle_0- 120)) ,   radius*sin(pi/180*(angle_0- 120))])  )
	color_points.append( array([Lp,  radius*cos(pi/180*(angle_0- 240)) ,   radius*sin(pi/180*(angle_0- 240))])  )

	center_array = array([0,center[0], center[1]])

	#patches = []
	#colors = []

	color_ps = array(color_points )
	delta = (1/nsegs)
	deltay = sin(60*pi/180)/nsegs
	from numpy.linalg import inv
	map_matrix = array([
						[0,   0,              1],
						[0.5, sin(60*pi/180), 1],
						[1.0, 0,              1]]).T


	imap = inv(map_matrix)

	for i in range(3):
		if verbose: print ('point',i)
		if verbose: print(map_matrix.T[i])
		fracs = imap.dot(map_matrix.T[i])
		if verbose: print ('fracs', fracs)
		color = color_ps[0]*fracs[0]+color_ps[1]*fracs[1]+color_ps[2]*fracs[2]
		if verbose: print('color',color)


	if verbose: print('upwards triangles')
	for iy in range(nsegs):
		for ix in range(nsegs-iy):
			vertices = [
						(delta*(    0.5*iy+ix), deltay *   iy),
						(delta*(0.5+0.5*iy+ix), deltay*(iy+1)),
						(delta*(1.0+0.5*iy+ix), deltay *   iy)]

			xc = delta*(0.5+0.5*iy+ix)
			yc = deltay * iy + 0.5*delta*tan(30*pi/180)
			fracs = imap.dot([xc,yc,1])
			#print (fracs)
			color = color_ps[0]*fracs[0]+color_ps[1]*fracs[1]+color_ps[2]*fracs[2] + center_array
			sRGB1_color = cspace_convert(color, "CIELab", "sRGB1")
			if use_deuteranomaly: sRGB1_color = clip(cspace_convert(sRGB1_color, cvd_space, "sRGB1"),0,1)
			ax.add_patch(  Polygon(vertices, color =clip(sRGB1_color,0,1), **triangle_style))


	if verbose: print('downwards triangles')
	for iy in range(nsegs-1):
		for ix in range(nsegs-iy-1):
			vertices = [
						(delta*(0.5+0.5*iy+ix), deltay *(iy+1)),
						(delta*(1.0+0.5*iy+ix), deltay * iy   ),
						(delta*(1.5+0.5*iy+ix), deltay *(iy+1))]
			xc = delta*(1.0+0.5*iy+ix)
			yc = deltay *(iy+1) - 0.5*delta*tan(30*pi/180)
			fracs = imap.dot([xc,yc,1])
			#print (fracs)
			color = color_ps[0]*fracs[0]+color_ps[1]*fracs[1]+color_ps[2]*fracs[2] + center_array
			sRGB1_color = cspace_convert(color, "CIELab", "sRGB1")
			if use_deuteranomaly: sRGB1_color = clip(cspace_convert(sRGB1_color, cvd_space, "sRGB1"),0,1)
			ax.add_patch(  Polygon(vertices, color =clip(sRGB1_color,0,1), **triangle_style))


	if len(labels) ==3:

		ax.text(   0,               -label_offset, labels[0], ha = 'center', va = 'top',    **font_options)
		ax.text(   1,               -label_offset, labels[1], ha = 'center', va = 'top',    **font_options)
		ax.text( 0.5, sin(60*pi/180)+label_offset, labels[2], ha = 'center', va = 'bottom', **font_options)

	ax.set_aspect('equal','box')
	#ax.set_xlim(0,1)
	#ax.set_ylim(0,sin(60*pi/180))
	ax.set_axis_off()


def make_sRGB1_image_in_CIELAB_space(
					Lp,
					ab_step_size = 2.0,
					amin = -128, amax = 128,
					bmin = -128, bmax = 128,
					color_for_undisplayable = [1.0, 1.0, 1.0, 0.0], L_scale = 1.0):

	from numpy import  array, ones, zeros, meshgrid, arange, linspace
	from colorspacious import cspace_convert
	### i didn't like this way, i want to ensure square pixels
	#da = (amax-amin)/na_grid
	#db = (bmax-bmin)/nb_grid
	#a_points = linspace(amin+da/2, amax-da/2, na_grid)
	#b_points = linspace(bmin+db/2, bmax-db/2, nb_grid)

	### better with square pixels
	da = db = ab_step_size
	a_points = arange(amin, amax, ab_step_size)
	b_points = arange(bmin, bmax, ab_step_size)

	######### setting up grids
	L_max = Lp*L_scale
	na_grid, nb_grid = len(a_points), len(b_points)
	lab_map = L_max * ones((nb_grid, na_grid, 3))
	a_grid, b_grid = meshgrid(a_points,b_points, indexing = 'xy')
	lab_map[:,:,1] = a_grid
	lab_map[:,:,2] = b_grid

	bg_color = array(color_for_undisplayable)

	sRGB1_map = zeros(( nb_grid, na_grid, 4)) # default is transparent
	sRGB1_map[:,:,0:3] = cspace_convert(lab_map, "CIELab", "sRGB1")
	for a_index in range(na_grid):
		for b_index in range(nb_grid):
			sRGB1_color = sRGB1_map[b_index,a_index,0:3]
			if  all(0<sRGB1_color) and all(sRGB1_color<1):
				sRGB1_map[b_index,a_index,3] = 1
			else:
				sRGB1_map[b_index,a_index] = bg_color

	return sRGB1_map


def CIELAB_triangle_with_LAB_slice(
						ax,
						Lp,
						radius,
						angle_0,
						center = (0.0, 0.0),
						ab_step_size = 2.0,
						amin = -128, amax = 128,
						bmin = -128, bmax = 128, L_scale = 1.0,
						use_deuteranomaly = False, color_for_undisplayable = [1.0, 1.0, 1.0, 0.0],
						):



	from numpy import  cos, sin, array, pi
	a_points = []
	b_points = []
	for i in range(3):
		a_points.append( L_scale * radius*cos(pi/180*(angle_0 - 120*i)) + center[0])
		b_points.append( L_scale * radius*sin(pi/180*(angle_0 - 120*i)) + center[1])
		ax.plot([0, a_points[-1]],  [0, b_points[-1]], color = 'dimgrey')

	ax.plot(a_points + [a_points[0]],  b_points + [b_points[0]], color = 'k')

	ax.set_xlabel('$a^{*}$')
	ax.set_ylabel('$b^{*}$', labelpad = -12)
	ax.minorticks_on()

	## these overide the input for tighter  fit to the triangle
	#buffer = 0.4
	#amin, amax = min(a_points) - buffer*radius, max(a_points) + buffer*radius
	#bmin, bmax = min(b_points) - buffer*radius, max(b_points) + buffer*radius


	LAB_slice_image = make_sRGB1_image_in_CIELAB_space( Lp = Lp,
							ab_step_size = ab_step_size,
							amin = amin, amax = amax,
							bmin = bmin, bmax = bmax, L_scale = L_scale, color_for_undisplayable = color_for_undisplayable  )
	extent = [	amin, amin + ab_step_size*LAB_slice_image.shape[1],
				bmin, bmin + ab_step_size*LAB_slice_image.shape[0]]

	if use_deuteranomaly : LAB_slice_image = deuteranomaly_with_alpha(LAB_slice_image)

	ax.imshow(LAB_slice_image , extent = extent, origin = 'lower')#, interpolation = 'bicubic')
	ax.set_xlim(amin, amax)
	ax.set_ylim(bmin, bmax)

def maximize_triangle_radius( Lp, angle_0, center = (0.0, 0.0), cut_off = 0.01,
				radius_min = 0.0,  radius_max = 256.0, verbose=False):

	from numpy import zeros, cos, sin, array, pi
	from colorspacious import cspace_convert
	from numpy import any as np_any
	angles = angle_0 + array([0,120,240])
	a_unit_vec = cos(pi/180*angles)
	b_unit_vec = sin(pi/180*angles)

	lab_sequence = zeros((3,3))
	lab_sequence[:,0] = Lp
	while (radius_max - radius_min) > cut_off:
		radius = 0.5*(radius_max + radius_min)

		lab_sequence[:,1] = radius*a_unit_vec + center[0]
		lab_sequence[:,2] = radius*b_unit_vec + center[1]
		sRGB1_sequence = cspace_convert(lab_sequence, "CIELab", "sRGB1")

		if np_any(sRGB1_sequence>1.0) or np_any(sRGB1_sequence<0.0):
			radius_is_safe = False
			if verbose:print (radius_min, radius, radius_max, radius_is_safe)
			radius_max = radius

		else:
			radius_is_safe = True
			if verbose:print (radius_min, radius, radius_max, radius_is_safe)
			radius_min = radius

	if radius_is_safe==False:
		radius = radius_min

	return radius


def maximize_triangle_radius_and_angle_0( Lp, angle_0, center = (0.0, 0.0),  angle_step = 1.0, angle_range = 120,
				cut_off = 0.01, radius_min = 0.0,  radius_max = 128.0, verbose=False):

	from numpy import arange
	angles = angle_0 + arange(-angle_range/2.0, angle_range/2.0, angle_step)

	best_angle = angle_0
	best_radius = maximize_triangle_radius( Lp, angle_0 = angle_0, center = center,  cut_off = cut_off ,
				radius_min = 0.0,  radius_max = 128.0, verbose=False)

	for angle in angles:
		radius = maximize_triangle_radius( Lp, angle_0 = angle, center = center, cut_off = cut_off ,
					radius_min = 0.0,  radius_max = 128.0, verbose=False)
		if radius > best_radius:
			best_angle = angle
			best_radius = radius

	return best_angle, best_radius


if __name__ == "__main__":


	from matplotlib.pylab import figure, gca, axis, savefig, show, imsave, imshow
	figure()
	ax=gca()
	sRGB1_color_triangle(ax,
					color_points = [[1.0, 0.0, 1.0], [0.0, 1.0, 1.0],
					[1.0, 1.0, 0.0]], nsegs = 7)


	savefig('sRGB_triangle.pdf',transparent = True)



	####################
	figure()
	ax=gca()

	Lp = 74
	radius = 45
	angle_0 = -60
	CIELAB_color_triangle(ax,
				Lp = Lp,
				radius = radius,
				angle_0 = angle_0,
				L_scale = 1.0,
				nsegs = 7)


	savefig('CIELab_triangle.pdf',transparent = True)

	##################
	figure()
	ax = gca()

	ab_step_size = 1.0
	amin, amax = -128, 128
	bmin, bmax = -128, 128

	Lp = 70
	angle_0 = -180
	L_scale = 1.0

	radius = maximize_triangle_radius( Lp, angle_0, verbose = False)
	print(radius)

	CIELAB_triangle_with_LAB_slice( ax,
				Lp = Lp,
				radius = radius,
				angle_0 = angle_0,
				ab_step_size = ab_step_size,
				amin = amin, amax = amax,
				bmin = bmin, bmax = bmax, L_scale = L_scale)
	ax.set_title("angle_0 = %.3f degrees" % angle_0)


	savefig('CIELAB_triangle_with_LAB_slice.pdf',transparent = True)
	############
	figure()
	ax = gca()

	Lp = 61.5

	angle_0, radius = maximize_triangle_radius_with_angle_0( Lp = Lp, angle_0 = angle_0, angle_step = 0.2)
	print( angle_0, radius)
	ax.set_title("angle_0 = %.3f degrees" % angle_0)
	CIELAB_triangle_with_LAB_slice( ax,
					Lp = Lp,
					radius = radius,
					angle_0 = angle_0,
					ab_step_size = ab_step_size,
					amin = amin, amax = amax,
					bmin = bmin, bmax = bmax, L_scale = L_scale)

	savefig('maximized_CIELAB_triangle_with_LAB_slice.pdf',transparent = True)


	###############

	figure()
	ax = gca()
	ax.grid(b=True, which='major', color='k', linestyle='-')
	ax.grid(b=True, which='minor', color='grey', linestyle='-')

	angle_step = 10 #0.1
	cut_off = 0.2 #0.01
	Lp_step = 5 #0.5

	angle_0 = -180
	angle_0_list = []
	radius_list = []
	from numpy import arange
	Lp_list = arange(100,0,-Lp_step)
	for Lp in  Lp_list:
		angle_0, radius = maximize_triangle_radius_with_angle_0( Lp = Lp, angle_0 = angle_0, angle_step = angle_step, cut_off = cut_off)

		angle_0_list.append(angle_0)
		radius_list.append(radius)

		print(Lp, angle_0, radius)

	ax.plot(Lp_list, radius_list)
	ax.set_xlabel('L*')
	ax.set_ylabel('Max Radius')

	ax.minorticks_on()

	#savefig('maximized_radius_vs_Lp.pdf',transparent = True)
	show()
