import numpy as np



delta = 0.1# 0.025
imsize = 3
x = y = np.arange(-imsize , imsize , delta)
X, Y = np.meshgrid(x, y)

#negative signals don't make sense

C1 = np.exp(-(X+1)**2 - (Y-1)**2) - 0.5*np.exp(-(X - 0)**2 - (Y - 0)**2)
C1 = C1 - C1.min()

C2 = np.exp(-2*(X-1)**2 - (Y/2)**2) - 1/(2.7)*np.exp(-(X/2 - 1/2)**2 - (Y - 1)**2)
C2 = C2 - C2.min()

C3 = np.exp(-(X-0.1)**2 - (Y+1)**2) - 0.8*np.exp(-(X + 1)**2 - (Y + 1)**2)
C3 = C3 - C3.min()

test_compositions = [C1,C2,C3]

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from matplotlib import cm

    C1 = (C1-C1.min()) / (C1.max() - C1.min())
    C2 = (C2-C2.min()) / (C2.max() - C2.min())
    C3 = (C3-C3.min()) / (C3.max() - C3.min())

    RGB_style = np.array([C1,C2,C3])
    RGB_style = np.moveaxis(RGB_style,0,-1)

    fig, axes = plt.subplots(ncols = 4)
    axes[0].set_title('Red')
    axes[0].imshow(C1, vmax=1.0, vmin=0, cmap = cm.gray,
        origin='lower', extent=[-imsize, imsize, -imsize, imsize])

    axes[1].set_title('Green')
    axes[1].imshow(C2, vmax=1.0, vmin=0, cmap = cm.gray,
        origin='lower', extent=[-imsize, imsize, -imsize, imsize])

    axes[2].set_title('Blue')
    axes[2].imshow(C3, vmax=1.0, vmin=0, cmap = cm.gray,
        origin='lower', extent=[-imsize, imsize, -imsize, imsize])

    axes[3].imshow(RGB_style, vmax=1.0, vmin=0,
        origin='lower', extent=[-imsize, imsize, -imsize, imsize])

    plt.show()
