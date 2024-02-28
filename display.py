import uproot as up
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import argparse
from os.path import join

CellWidthX = 2.5
CellWidthY = 2.5
LayerThick = 4.0
nCellX = 21
nCellY = 21
nLayer = 11


def read_file(fname: str, tree: str, event_index: int, staggered: bool):
    with up.open(fname) as f:
        tree = f[tree]

        Hit_X = tree['Hit_X'].array(library='np')
        Hit_Y = tree['Hit_Y'].array(library='np')
        Hit_Z = tree['Hit_Z'].array(library='np')
        Hit_Energy = tree['Hit_Energy'].array(library='np')

        xtemp = Hit_X[event_index]
        ytemp = Hit_Y[event_index]
        x = np.zeros_like(xtemp)
        y = np.zeros_like(ytemp)
        z = np.round(Hit_Z[event_index] / LayerThick).astype(int)
        energy = Hit_Energy[event_index]

        assert len(x) == len(y)
        assert len(x) == len(z)
        assert len(x) == len(energy)

        if staggered:
            for i in np.arange(len(z)):
                if z[i] % 2 == 0:
                    x[i] = np.round(xtemp[i] / CellWidthX - 0.25).astype(int)
                    y[i] = np.round(ytemp[i] / CellWidthY - 0.25).astype(int)
#                    if abs(xtemp[i] / CellWidthX - 0.25) >= 10 or abs(ytemp[i] / CellWidthY - 0.25) >= 10:
#                        print(xtemp[i] / CellWidthX - 0.25, ytemp[i] / CellWidthY - 0.25)
                else:
                    x[i] = np.round(xtemp[i] / CellWidthX + 0.25).astype(int)
                    y[i] = np.round(ytemp[i] / CellWidthY + 0.25).astype(int)
#                    if abs(xtemp[i] / CellWidthX + 0.25) >= 10 or abs(ytemp[i] / CellWidthY + 0.25) >= 10:
#                        print(xtemp[i] / CellWidthX + 0.25, ytemp[i] / CellWidthY + 0.25)
        else:
            x[i] = np.round(xtemp[i] / CellWidthX).astype(int)
            y[i] = np.round(ytemp[i] / CellWidthY).astype(int)

        znew, ynew, xnew, enew = (np.array(a) for a in zip(*sorted(zip(z, y, x, energy), reverse = True)))

        return xnew, ynew, znew, enew


def plot(fname: str, tree: str, event_index: int, title: str, staggered: bool):
    x, y, z, energy = read_file(fname, tree, event_index, staggered)
    energy_norm = energy / np.max(energy)

    nhits = len(x)

    fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
    cmap = cm.OrRd

    if staggered:
        WidthX = (nCellX + 0.5) * CellWidthX
        WidthY = (nCellY + 0.5) * CellWidthY
    else:
        WidthX = nCellX * CellWidthX
        WidthY = nCellY * CellWidthY

    plt.gca().set_box_aspect((LayerThick * (nLayer + 1) / WidthX, 1, LayerThick * (nLayer + 1) / WidthY))

    for i in np.arange(nhits):
        if staggered:
            if z[i] % 2 == 0:
                xnew = np.arange(x[i] - 1, x[i] + 1) + 0.75
                ynew = np.arange(y[i] - 1, y[i] + 1) + 0.75
            else:
                xnew = np.arange(x[i] - 1, x[i] + 1) + 0.25
                ynew = np.arange(y[i] - 1, y[i] + 1) + 0.25
        else:
            xnew = np.arange(x[i] - 1, x[i] + 1)
            ynew = np.arange(y[i] - 1, y[i] + 1)

        xnew, ynew = np.meshgrid(xnew, ynew)
        znew = z[i] * np.ones(xnew.shape)
        enew = energy_norm[i] * np.ones(xnew.shape)

        ax.plot_surface(xnew, znew, ynew, cmap=cmap, facecolors=cmap(enew), edgecolor='k', alpha=0.8, lw=0.05, rstride=1, cstride=1, antialiased=False)

    ax.set_xlim(-0.5 * nCellX, 0.5 * nCellY)
    ax.set_ylim(0, nLayer + 1)
    ax.set_zlim(-0.5 * nCellY, 0.5 * nCellY)
    ax.set_xticks(np.linspace(-10, 10, 5), CellWidthX * np.linspace(-10, 10, 5))
    ax.set_yticks(np.linspace(0, nLayer + 1, 5))
    ax.set_zticks(np.linspace(-10, 10, 5), CellWidthY * np.linspace(-10, 10, 5))
    ax.set_aspect(aspect='equalxz')
    ax.grid(False)

    fig.suptitle(title, size='xx-large')
    ax.invert_xaxis()
    ax.set_xlabel("X [cm]", size='x-large')
    ax.set_ylabel("Z [layer]", size='x-large')
    ax.set_zlabel("Y [cm]", size='x-large')

    m = plt.cm.ScalarMappable(cmap=cmap)
    m.set_array(energy)
    plt.colorbar(m, pad=0.2, ax=plt.gca()).set_label(label="Hit Energy [MeV]", size='x-large')

    ax.view_init(elev=20, azim=-35, roll=0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--file",      type=str, default='',   required=True,  help="Input ROOT file")
    parser.add_argument("-t", "--tree",      type=str, default='dp', help="Input tree name (default: dp)")
    parser.add_argument("-g", "--staggered", type=int, default=1,    choices=[0, 1], help="Staggered structure")
    parser.add_argument("-i", "--title",     type=str, default='',   help="Title of display figure")
    parser.add_argument("-e", "--event",     type=int, default=0,    help="The event to be displayed")
    parser.add_argument("-d", "--dir",       type=str, default=None, help="Directory to save the plot")
    parser.add_argument("-o", "--output",    type=str, default=None, help="Output file name")
    parser.add_argument("-s", "--show",      type=int, default=1,    choices=[0, 1], help="Instantly display or not")
    args = parser.parse_args()

    filename = args.file
    tree = args.tree
    staggered = args.staggered
    title = args.title
    event_index = args.event
    save_dir = args.dir
    output = args.output
    show = args.show

    plot(filename, tree, event_index, title, staggered)

    if save_dir and output:
        plt.savefig(join(save_dir, output))
        print("Figure", join(save_dir, output), "successfully created!")

    if show:
        plt.show()
