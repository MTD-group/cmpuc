

default_corner_points = 256
import numpy as np
from colorspacious import cspace_convert

default_color_space = "CIELab"
cvd_space = {"name": "sRGB1+CVD", "cvd_type":"deuteranomaly", "severity": 50}


# I'm not sure that 'antialiased' is needed to be off, but it's suggested online
# there needs to be *some* line width or the vectorized representations like pdf have stupid lines
# the joinstyle of bevel or round makes the lines not shoot outwards like miter
triangle_style = {'antialiased': False, 'linewidth' : 0.1, 'joinstyle': 'bevel' }




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


class uniform_chemical_map:
	def __init__(self, radius, angle_0, L_plane, center = (0.0,0.0), color_space = default_color_space):
		# center is in the L_plane and also get scaled to 0 for lower signals
		self.radius  = radius
		self.angle_0 = angle_0
		self.L_plane = L_plane
		self.center = np.asarray(center)
		self.color_space = color_space

	def get_color_vectors(self, norm = 1.0):
		color_vectors = []
		for i in range(3):
			color_vectors.append(
				norm * np.array([
					self.L_plane,
					self.radius*np.cos(np.pi/180*(self.angle_0 + i*120) ),
					self.radius*np.sin(np.pi/180*(self.angle_0 + i*120) )
					])  )
		return color_vectors

	def get_color_points(self, norm = 1.0):
		cv = self.get_color_vectors(norm)
		center_array = np.array([0, self.center[0],self.center[1]])
		color_points = []
		for i in range(3):
			color_points.append(cv[i]+center_array)
		return color_points

	def get_sRGB1_color_points(self, use_deuteranomaly = False, norm = 1.0):
		cp = self.get_color_points(norm)
		sRGB1_color_points = []
		for i in range(3):
			sRGB1_color_points.append(cspace_convert(cp[i], self.color_space, "sRGB1"))
		if use_deuteranomaly:
			for i in range(3):
				sRGB1_color_points[i] = np.clip(cspace_convert(sRGB1_color_points[i], cvd_space, "sRGB1"),0,1)
		return sRGB1_color_points


	def __call__(self, data, auto_tune_L_plane_to_data = False,  contrast_boost = 1.0, verbose = True, use_deuteranomaly = False, scale_alpha = False ):


		sRGB1_cp = self.get_sRGB1_color_points()
		if verbose:
			for i in range(3):
				print('sRGB1 color point %i:' %i, "(%.3f %.3f %.3f)"% tuple(sRGB1_cp[i]))


		normed_data = []
		for data_index in range(3):
			normed_data.append(data[data_index]/data[data_index].max() )

		overall_norm = normed_data[0]+normed_data[1]+normed_data[2]
		overall_norm = overall_norm/overall_norm.max()

		data_shape = data[0].shape
		Lab_map =  np.zeros(list(data_shape)+[3])
		sRGB1_map = np.ones(list(data_shape)+[4])

		color_vectors = self.get_color_vectors()
		center_vector = np.array([0, self.center[0],self.center[1]])

		for i in range(data_shape[0]):
			for j in range(data_shape[1]):
				Lab_map[i,j,0] = self.L_plane * overall_norm[i,j]
				for data_index in  range(3):
					# color points are really vectors, the normed_data is the projection along each one, and overall_norm is the overall scalling
					#if True:
					# taken out of the loop since all vectors have the same overall_norm for L
					#Lab_map[i,j,0] += (color_vectors[data_index][0] * normed_data[data_index][i,j] + center_vector[0])  * overall_norm[i,j]
					Lab_map[i,j,1] += (color_vectors[data_index][1] * normed_data[data_index][i,j] + center_vector[1])  * overall_norm[i,j]
					Lab_map[i,j,2] += (color_vectors[data_index][2] * normed_data[data_index][i,j] + center_vector[2])  * overall_norm[i,j]
					#else:
					#	Lab_map[i,j,1] += self.color_vectors[data_index][1] * normed_data[data_index][i,j] * overall_norm[i,j] + self.center[0]
					#	Lab_map[i,j,2] += self.color_vectors[data_index][2] * normed_data[data_index][i,j] * overall_norm[i,j] + self.center[1]


		L_scale = 1.0
		Lab_test = L_scale * Lab_map

		if auto_tune_L_plane_to_data:
			max_L = Lab_map[:,:,0].max()
			if verbose: print ('Max L', Lab_map[:,:,0].max())

			### now we can scale it all to fit tighter
			# step 1, do a log search upwards to find an upper cutoff
			sRGB1_test = cspace_convert(Lab_test, self.color_space, "sRGB1")
			while np.any(sRGB1_test>1.0) == False: # loop to find an upper bound
				L_scale = L_scale*1.2 # doubling is somewhat dumb, since we know we are on the right order of magnitude
				Lab_test[:,:,0] = L_scale * Lab_map[:,:,0]
				sRGB1_test = cspace_convert(Lab_test, self.color_space, "sRGB1")


			# this finds the perfect scalling of L for the color points
			cut_off = 1.0/1024 # half a 16 bit color depth unit, very accurate
			L_scale_max = L_scale
			L_scale_min = 0.0
			sRGB1_max = sRGB1_test.max()
			while ( 1.0 - sRGB1_max ) > cut_off or sRGB1_max > 1:
				L_scale = 0.5*(L_scale_max + L_scale_min)
				Lab_test[:,:,0] = L_scale * Lab_map[:,:,0]
				sRGB1_test = cspace_convert(Lab_test, self.color_space, "sRGB1")
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

			self.L_plane = L_scale * self.L_plane
			Lab_test[:,:,0] = L_scale * Lab_map[:,:,0]


		if contrast_boost==1.0:
			sRGB1_test = cspace_convert( Lab_test, self.color_space, "sRGB1")
		else:
			sRGB1_test = cspace_convert(contrast_boost * Lab_test, self.color_space, "sRGB1")
			self.center  = contrast_boost * self.center
			self.radius  = contrast_boost * self.radius
			self.L_plane = contrast_boost * self.L_plane

		sRGB1_map[:,:,0:3] = sRGB1_test[:,:,0:3]
		print ('Max sRGB1', sRGB1_map[:,:,0:3].max())

		if scale_alpha:
			sRGB1_map[:,:,3] =  overall_norm

		if verbose: bounds_test(sRGB1_map)
		if use_deuteranomaly:
			sRGB1_map = np.clip(cspace_convert(sRGB1_map, cvd_space, "sRGB1"),0,1)
		return  sRGB1_map

	def maximize_triangle_radius(self, cut_off = 0.01, radius_min = 0.0,
							radius_max = 256.0, verbose=False):

		angles = self.angle_0 + np.array([0,120,240])
		a_unit_vec = np.cos(np.pi/180*angles)
		b_unit_vec = np.sin(np.pi/180*angles)

		lab_sequence = np.zeros((3,3))
		lab_sequence[:,0] = self.L_plane
		while (radius_max - radius_min) > cut_off:
			radius = 0.5*(radius_max + radius_min)

			lab_sequence[:,1] = radius*a_unit_vec + self.center[0]
			lab_sequence[:,2] = radius*b_unit_vec + self.center[1]
			sRGB1_sequence = cspace_convert(lab_sequence, self.color_space, "sRGB1")

			if np.any(sRGB1_sequence>1.0) or np.any(sRGB1_sequence<0.0):
				radius_is_safe = False
				if verbose:print (radius_min, radius, radius_max, radius_is_safe)
				radius_max = radius

			else:
				radius_is_safe = True
				if verbose:print (radius_min, radius, radius_max, radius_is_safe)
				radius_min = radius

		if radius_is_safe==False:
			radius = radius_min

		self.radius = radius
		return radius


	def maximize_triangle_radius_and_angle_0(self,  angle_step = 0.2, angle_range = 120,
					cut_off = 0.01, radius_min = 0.0,  radius_max = 128.0, verbose=False):

		angles = self.angle_0 + np.arange(-angle_range/2.0, angle_range/2.0, angle_step)

		best_angle = self.angle_0
		best_radius = self.maximize_triangle_radius( cut_off = cut_off ,
					radius_min = 0.0,  radius_max = 128.0, verbose=False)

		for angle in angles:
			self.angle_0 = angle
			radius = self.maximize_triangle_radius(
						radius_min = 0.0,  radius_max = 128.0, verbose=False)
			if radius > best_radius:
				best_angle = angle
				best_radius = radius

		self.radius = best_radius
		self.angle_0 = best_angle
		return best_angle, best_radius



def deuteranomaly_with_alpha(sRGB1_map):
	#cvd_space = {"name": "sRGB1+CVD", "cvd_type":"deuteranomaly", "severity": 50}
	sRGB1_map_d = np.zeros(sRGB1_map.shape)
	sRGB1_map_d[:,:,0:3] = np.clip(cspace_convert(sRGB1_map[:,:,0:3], cvd_space, "sRGB1"),0,1)
	sRGB1_map_d[:,:,3] = sRGB1_map[:,:,3]
	return sRGB1_map_d



def isoluminant_triangle(ax,
				chemical_map,
				nsegs = 20,
				labels = ['Part 1','Part 2','Part 3' ],
				norm = 1.0,
				font_options = {},
				label_offset = 0.02,
				use_deuteranomaly = False,verbose = False):
	from matplotlib.patches import Polygon


	cv = chemical_map.get_color_vectors()

	center_array = np.array([0,chemical_map.center[0], chemical_map.center[1]])

	#patches = []
	#colors = []

	color_ps = np.array(cv )
	delta = (1/nsegs)
	deltay = np.sin(60*np.pi/180)/nsegs
	from numpy.linalg import inv
	map_matrix = np.array([
						[0,   0,              1],
						[0.5, np.sin(60*np.pi/180), 1],
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
			yc = deltay * iy + 0.5*delta*np.tan(30*np.pi/180)
			fracs = imap.dot([xc,yc,1])
			#print (fracs)
			color = norm * (color_ps[0]*fracs[0]+color_ps[1]*fracs[1]+color_ps[2]*fracs[2] + center_array)
			sRGB1_color = cspace_convert(color, chemical_map.color_space, "sRGB1")
			if use_deuteranomaly: sRGB1_color = np.clip(cspace_convert(sRGB1_color, cvd_space, "sRGB1"),0,1)
			ax.add_patch(  Polygon(vertices, color =np.clip(sRGB1_color,0,1), **triangle_style))


	if verbose: print('downwards triangles')
	for iy in range(nsegs-1):
		for ix in range(nsegs-iy-1):
			vertices = [
						(delta*(0.5+0.5*iy+ix), deltay *(iy+1)),
						(delta*(1.0+0.5*iy+ix), deltay * iy   ),
						(delta*(1.5+0.5*iy+ix), deltay *(iy+1))]
			xc = delta*(1.0+0.5*iy+ix)
			yc = deltay *(iy+1) - 0.5*delta*np.tan(30*np.pi/180)
			fracs = imap.dot([xc,yc,1])
			#print (fracs)
			color = norm * (color_ps[0]*fracs[0]+color_ps[1]*fracs[1]+color_ps[2]*fracs[2] + center_array)
			sRGB1_color = cspace_convert(color, chemical_map.color_space, "sRGB1")
			if use_deuteranomaly: sRGB1_color = np.clip(cspace_convert(sRGB1_color, cvd_space, "sRGB1"),0,1)
			ax.add_patch(  Polygon(vertices, color =np.clip(sRGB1_color,0,1), **triangle_style))


	if len(labels) ==3:

		ax.text(   0,               -label_offset, labels[0], ha = 'center', va = 'top',    **font_options)
		ax.text(   1,               -label_offset, labels[1], ha = 'center', va = 'top',    **font_options)
		ax.text( 0.5, np.sin(60*np.pi/180)+label_offset, labels[2], ha = 'center', va = 'bottom', **font_options)

	ax.set_aspect('equal','box')
	#ax.set_xlim(0,1)
	#ax.set_ylim(0,sin(60*pi/180))
	ax.set_axis_off()


def isoluminant_slice(
					L_plane,
					ab_step_size = 1.0,
					amin = -128, amax = 128,
					bmin = -128, bmax = 128,
					color_for_undisplayable = [1.0, 1.0, 1.0, 0.0]):

	#from numpy import  array, ones, zeros, meshgrid, arange, linspace
	from colorspacious import cspace_convert
	### i didn't like this way, i want to ensure square pixels
	#da = (amax-amin)/na_grid
	#db = (bmax-bmin)/nb_grid
	#a_points = linspace(amin+da/2, amax-da/2, na_grid)
	#b_points = linspace(bmin+db/2, bmax-db/2, nb_grid)

	### better with square pixels
	da = db = ab_step_size
	a_points = np.arange(amin, amax, ab_step_size)
	b_points = np.arange(bmin, bmax, ab_step_size)

	######### setting up grids

	na_grid, nb_grid = len(a_points), len(b_points)
	lab_map = L_plane * np.ones((nb_grid, na_grid, 3))
	a_grid, b_grid = np.meshgrid(a_points,b_points, indexing = 'xy')
	lab_map[:,:,1] = a_grid
	lab_map[:,:,2] = b_grid

	bg_color = np.array(color_for_undisplayable)

	sRGB1_image = np.zeros(( nb_grid, na_grid, 4)) # default is transparent
	sRGB1_image[:,:,0:3] = cspace_convert(lab_map, default_color_space, "sRGB1")
	for a_index in range(na_grid):
		for b_index in range(nb_grid):
			sRGB1_color = sRGB1_image[b_index,a_index,0:3]
			if  all(0<sRGB1_color) and all(sRGB1_color<1):
				sRGB1_image[b_index,a_index,3] = 1
			else:
				sRGB1_image[b_index,a_index] = bg_color

	return sRGB1_image


def triangle_on_isoluminant_slice(
						ax,
						chemical_map,
						ab_step_size = 1.0,
						amin = -128, amax = 128,
						bmin = -128, bmax = 128,
						use_deuteranomaly = False, color_for_undisplayable = [1.0, 1.0, 1.0, 0.0],
						):

	cv = chemical_map.get_color_vectors()
	center_array = np.array([0,chemical_map.center[0], chemical_map.center[1]])

	a_points = []
	b_points = []
	for i in range(3):
		point = cv[i]+center_array
		a_points.append( point[1])
		b_points.append( point[2])
		ax.plot([0, a_points[-1]],  [0, b_points[-1]], color = 'dimgrey')

	ax.plot(a_points + [a_points[0]],  b_points + [b_points[0]], color = 'k')

	ax.set_xlabel('$a^{*}$')
	ax.set_ylabel('$b^{*}$', labelpad = -12)
	ax.minorticks_on()

	## these overide the input for tighter  fit to the triangle
	#buffer = 0.4
	#amin, amax = min(a_points) - buffer*radius, max(a_points) + buffer*radius
	#bmin, bmax = min(b_points) - buffer*radius, max(b_points) + buffer*radius


	LAB_slice_image = isoluminant_slice( L_plane = chemical_map.L_plane,
							ab_step_size = ab_step_size,
							amin = amin, amax = amax,
							bmin = bmin, bmax = bmax, color_for_undisplayable = color_for_undisplayable  )
	extent = [	amin, amin + ab_step_size*LAB_slice_image.shape[1],
				bmin, bmin + ab_step_size*LAB_slice_image.shape[0]]

	if use_deuteranomaly : LAB_slice_image = deuteranomaly_with_alpha(LAB_slice_image)

	ax.imshow(LAB_slice_image , extent = extent, origin = 'lower')#, interpolation = 'bicubic')
	ax.set_xlim(amin, amax)
	ax.set_ylim(bmin, bmax)


def sRGB1_colormap(data,
					elements,
					ax,
					color_points):

	maxes = []
	for dat in data:
		maxes.append(dat.max())

	data_shape = data[0].shape
	sRGB1_map = np.zeros(list(data_shape)+[4])
	sRGB1_color_points = np.array(color_points)

	for i in range(data_shape[0]):
		for j in range(data_shape[1]):
			sRGB1_map[i,j,3] = 1.0
			for data_index in  range(3):
				sRGB1_map[i,j,0:3] += sRGB1_color_points[data_index] * data[data_index][i,j]/maxes[data_index]

	ax.imshow(np.clip(sRGB1_map,0,1))


	return np.clip(sRGB1_map,0,1)


def sRGB1_color_triangle(ax,  color_points = [[1.0, 0.0, 1.0], [0.0, 1.0, 1.0], [1.0, 1.0, 0.0]] , nsegs = 20, use_deuteranomaly = False, verbose = False):
	from matplotlib.patches import Polygon

	patches = []
	colors = []

	color_ps = np.array([color_points[0], color_points[2], color_points[1]])
	#color_ps = np.array(color_points)
	delta = (1/nsegs)
	deltay = np.sin(60*np.pi/180)/nsegs

	cvd_space = {"name": "sRGB1+CVD", "cvd_type":"deuteranomaly", "severity": 50}
	from colorspacious import cspace_convert

	from numpy.linalg import inv
	map_matrix = np.array([
						[0,   0,              1],
						[0.5, np.sin(60*np.pi/180), 1],
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
			yc = deltay * iy + 0.5*delta*np.tan(30*np.pi/180)
			fracs = imap.dot([xc,yc,1])
			#print (fracs)
			color = color_ps[0]*fracs[0]+color_ps[1]*fracs[1]+color_ps[2]*fracs[2]
			if use_deuteranomaly: color = cspace_convert(color, cvd_space, "sRGB1")
			ax.add_patch(  Polygon(vertices, color =np.clip(color,0,1), **triangle_style))


	if verbose: print('downwards triangles')
	for iy in range(nsegs-1):
		for ix in range(nsegs-iy-1):
			vertices = [
						(delta*(0.5+0.5*iy+ix), deltay *(iy+1)),
						(delta*(1.0+0.5*iy+ix), deltay * iy   ),
						(delta*(1.5+0.5*iy+ix), deltay *(iy+1))]
			xc = delta*(1.0+0.5*iy+ix)
			yc = deltay *(iy+1) - 0.5*delta*np.tan(30*np.pi/180)
			fracs = imap.dot([xc,yc,1])
			#print (fracs)
			color = color_ps[0]*fracs[0]+color_ps[1]*fracs[1]+color_ps[2]*fracs[2]
			if use_deuteranomaly: color = np.clip(cspace_convert(color, cvd_space, "sRGB1"),0,1)
			ax.add_patch(  Polygon(vertices, color =np.clip(color,0,1), **triangle_style))

	#ax.set_aspect('equal', 'box')


	ax.set_xlim(0,1)
	ax.set_ylim(0,np.sin(60*np.pi/180))
	ax.set_axis_off()
	ax.set_aspect('equal')




if __name__ == "__main__":


	from matplotlib.pylab import figure, gca, axis, savefig, show, imsave, imshow
	figure()
	ax=gca()
	sRGB1_color_triangle(ax,
					color_points = [[1.0, 0.0, 1.0], [0.0, 1.0, 1.0],
					[1.0, 1.0, 0.0]], nsegs = 7)


	#savefig('sRGB_triangle.pdf',transparent = True)



	####################
	figure()
	ax=gca()

	L_plane = 74
	# First hue angle of triangle vertices
	angle_0 = -60.0
	## a radius for the triangle
	radius = 45
	###### these you'll likely never change
	## Center of Triangle
	center = (0,0)

	my_map =   uniform_chemical_map(
								L_plane = L_plane,
								angle_0 = angle_0,
								radius = radius, center = center)



	#savefig('CIELab_triangle.pdf',transparent = True)

	##################
	figure()
	ax = gca()


	radius = my_map.maximize_triangle_radius()
	print(radius)

	triangle_on_isoluminant_slice(ax, my_map)
	ax.set_title("angle_0 = %.3f degrees" % angle_0)


	#savefig('CIELAB_triangle_with_LAB_slice.pdf',transparent = True)
	############
	figure()
	ax = gca()

	angle_0, radius = my_map.maximize_triangle_radius_and_angle_0()
	print( angle_0, radius)
	ax.set_title("angle_0 = %.3f degrees" % angle_0)
	triangle_on_isoluminant_slice(ax, my_map)

	#savefig('maximized_CIELAB_triangle_with_LAB_slice.pdf',transparent = True)


	###############

	figure()
	ax = gca()
	ax.grid(b=True, which='major', color='k', linestyle='-')
	ax.grid(b=True, which='minor', color='grey', linestyle='-')



	Lp_step = 10 #0.5

	angle_0 = -180
	angle_0_list = []
	radius_list = []
	from numpy import arange
	Lp_list = arange(95,0,-Lp_step)
	for Lp in  Lp_list:
		my_map.L_plane = Lp
		angle_0, radius = my_map.maximize_triangle_radius_and_angle_0()

		angle_0_list.append(angle_0)
		radius_list.append(radius)

		print(Lp, angle_0, radius)

	ax.plot(Lp_list, radius_list)
	ax.set_xlabel('L*')
	ax.set_ylabel('Max Radius')

	ax.minorticks_on()

	#savefig('maximized_radius_vs_Lp.pdf',transparent = True)
	show()
