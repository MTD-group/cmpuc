
import numpy as np
default_corner_points = 128


def corner_sRGB1_map(
				data,
				elements,
				ax,
				color_points,
				npts = default_corner_points,
				bg_color = np.array([1, 1, 1, 0]) ):


	from numpy import linspace
	values = linspace(0,1,npts)

	sRGB1_map = np.ones((npts,npts,4))
	fractions = np.zeros((npts,npts,3))
	sRGB1_color_points = np.array(color_points)
	bg_color = np.array( [1,1,1,0])

	for i in range(npts):
		iv = values[i]
		for j in range(npts-i):
			jv = values[j]

			fraction = np.array( [(1-jv-iv), iv, jv ])
			#print(iv, jv, fraction)
			sRGB1_map[j,i] = np.array([0,0,0,1])
			for data_index in  range(3):
				sRGB1_map[j,i,0:3] += sRGB1_color_points[data_index] * fraction[data_index]

		for j in range(npts-i,npts):
			sRGB1_map[j,i] = bg_color

	maxes = []
	for dat in data:
		maxes.append(dat.max())
	ax.imshow(np.clip(sRGB1_map,0,1), origin = 'lower', extent=[0,1,0,1], interpolation = 'bicubic')
	ax.text(0,0, elements[0]+'\n%0.1f%%'% maxes[0], ha = 'left')
	ax.text(1,0, elements[1]+'\n%0.1f%%'% maxes[1], ha = 'right')
	ax.text(0,1, elements[2]+'\n%0.1f%%'% maxes[2], ha = 'left', va = 'top')

	return np.clip(sRGB1_map,0,1)


def corner_map(	data,
				elements,
				ax,
				Lp,
				radius,
				angle_0,
				L_scale,
				npts = default_corner_points,
				bg_color = np.array([1, 1, 1, 0]) ):

	color_points = []
	color_points.append( np.array([Lp,  radius*np.cos(np.pi/180*(angle_0+  0))  ,   radius*np.sin(np.pi/180*(angle_0+   0))])  )
	color_points.append( np.array([Lp,  radius*np.cos(np.pi/180*(angle_0+ 120)) ,   radius*np.sin(np.pi/180*(angle_0+ 120))])  )
	color_points.append( np.array([Lp,  radius*np.cos(np.pi/180*(angle_0+ 240)) ,   radius*np.sin(np.pi/180*(angle_0+ 240))])  )

	from numpy import linspace
	values = linspace(0,1,npts)

	Lab_guide = np.zeros((npts,npts,3))
	fractions = np.zeros((npts,npts,3))
	sRGB1_map = np.ones((npts,npts,4))

	for i in range(npts):
		iv = values[i]
		for j in range(npts-i):
			jv = values[j]

			fraction = np.array( [(1-jv-iv), iv, jv ])
			#print(iv, jv, fraction)
			for data_index in  range(3):
				Lab_guide[j,i] += color_points[data_index] * fraction[data_index]* L_scale

	sRGB1_map[:,:,0:3] = cspace_convert(Lab_guide, default_color_space, "sRGB1")
	for i in range(npts):
		for j in range(npts-i,npts):
			sRGB1_map[j,i] = bg_color

	maxes = []
	for dat in data:
		maxes.append(dat.max())

	ax.imshow(np.clip(sRGB1_map,0,1), origin = 'lower', extent=[0,1,0,1], interpolation = 'bicubic')
	ax.text(0,0, elements[0]+'\n%0.1f%%'% maxes[0], ha = 'left')
	ax.text(1,0, elements[1]+'\n%0.1f%%'% maxes[1], ha = 'right')
	ax.text(0,1, elements[2]+'\n%0.1f%%'% maxes[2], ha = 'left', va = 'top')

	return np.clip(sRGB1_map,0,1)
