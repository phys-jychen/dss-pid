import numpy as np
import matplotlib.pyplot as plt
import uproot as up
import argparse
from os.path import join

CellWidthX = 2.5
CellWidthY = 2.5
LayerThick = 4
nCellX = 21
nCellY = 21
nLayer = 11


def read_file(fname: str, tree: str, event_index: int):
    with up.open(fname) as f:
        tree = f[tree]

        Hit_X = tree['Hit_X'].array(library='np')[event_index]
        Hit_Y = tree['Hit_Y'].array(library='np')[event_index]
        Hit_Z = tree['Hit_Z'].array(library='np')[event_index]
        Hit_Energy = tree['Hit_Energy'].array(library='np')[event_index]

        assert len(Hit_X) == len(Hit_Y)
        assert len(Hit_X) == len(Hit_Z)
        assert len(Hit_X) == len(Hit_Energy)

        return Hit_X, Hit_Y, Hit_Z, Hit_Energy


def plot(fname: str, tree: str, event_index: int, title: str, staggered: bool):
    Hit_X, Hit_Y, Hit_Z, Hit_Energy = read_file(fname, tree, event_index)

    d_xy = dict()
    x_xy = np.array([])
    y_xy = np.array([])
    energy_xy = np.array([])

    # d_xz = dict()
    # x_xz = np.array([])
    # y_xz = np.array([])
    # energy_xz = np.array([])

    # d_yz = dict()
    # x_yz = np.array([])
    # y_yz = np.array([])
    # energy_yz = np.array([])

    # d_x = dict()
    # x = np.array([])
    # energy_x = np.array([])

    # d_y = dict()
    # y = np.array([])
    # energy_y = np.array([])

    for i in np.arange(len(Hit_X)):
        if (Hit_X[i], Hit_Y[i]) not in d_xy:
            d_xy[Hit_X[i], Hit_Y[i]] = Hit_Energy[i]
        else:
            d_xy[Hit_X[i], Hit_Y[i]] += Hit_Energy[i]

        # if (Hit_X[i], Hit_Z[i]) not in d_xz:
        #     d_xz[Hit_X[i], Hit_Z[i]] = Hit_Energy[i]
        # else:
        #     d_xz[Hit_X[i], Hit_Z[i]] += Hit_Energy[i]

        # if (Hit_Y[i], Hit_Z[i]) not in d_yz:
        #     d_yz[Hit_Y[i], Hit_Z[i]] = Hit_Energy[i]
        # else:
        #     d_yz[Hit_Y[i], Hit_Z[i]] += Hit_Energy[i]

        # if Hit_X[i] not in d_x:
        #     d_x[Hit_X[i]] = Hit_Energy[i]
        # else:
        #     d_x[Hit_X[i]] += Hit_Energy[i]

        # if Hit_Y[i] not in d_y:
        #     d_y[Hit_Y[i]] = Hit_Energy[i]
        # else:
        #     d_y[Hit_Y[i]] += Hit_Energy[i]

    for i in d_xy:
        x_xy = np.append(x_xy, i[0])
        y_xy = np.append(y_xy, i[1])
        energy_xy = np.append(energy_xy, d_xy[i])

    # for i in d_xz:
    #     x_xz = np.append(x_xz, i[0])
    #     y_xz = np.append(y_xz, i[1])
    #     energy_xz = np.append(energy_xz, d_xz[i])

    # for i in d_yz:
    #     x_yz = np.append(x_yz, i[0])
    #     y_yz = np.append(y_yz, i[1])
    #     energy_yz = np.append(energy_yz, d_yz[i])

    # for i in d_x:
    #     x = np.append(x, i)
    #     energy_x = np.append(energy_x, d_x[i])

    # for i in d_y:
    #     y = np.append(y, i)
    #     energy_y = np.append(energy_y, d_y[i])

    WidthX = (nCellX + 0.5 * staggered) * CellWidthX
    WidthY = (nCellY + 0.5 * staggered) * CellWidthY
    WidthZ = nLayer * LayerThick

    ratioX = WidthX / WidthY
    ratioY = 1
    ratioZ = 1

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.gca().set_box_aspect((ratioX, ratioY, ratioZ))

    binx = (staggered + 1) * nCellX
    biny = (staggered + 1) * nCellY

    hist, xedges, yedges = np.histogram2d(x_xy, y_xy, bins=[binx, biny], range=[[-0.5 * WidthX, 0.5 * WidthX], [-0.5 * WidthY, 0.5 * WidthY]],  weights=energy_xy)
    xpos, ypos = np.meshgrid(xedges[:-1], yedges[:-1], indexing='ij')
    xpos = xpos.ravel()
    ypos = ypos.ravel()
    zpos = 0

    dx = CellWidthX / (staggered + 1)
    dy = CellWidthY / (staggered + 1)
    dz = hist.ravel()

    cmap = plt.cm.RdPu
    max_height = np.max(dz)
    min_height = np.min(dz)
    rgba = [cmap((k - min_height) / (max_height - min_height)) for k in dz]

    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=rgba, edgecolor='none')

    clb = plt.cm.ScalarMappable(cmap=cmap)
    clb.set_array(dz)
    fig.colorbar(clb, pad=0.2, ax=plt.gca())

    ax.xaxis._axinfo['grid'].update(linestyle=':', linewidth=0.5)
    ax.yaxis._axinfo['grid'].update(linestyle=':', linewidth=0.5)
    ax.zaxis._axinfo['grid'].update(linestyle=':', linewidth=0.5)
    ax.set_xlim(-0.5 * WidthX, 0.5 * WidthX)
    ax.set_ylim(-0.5 * WidthY, 0.5 * WidthY)

    ax.set_xlabel("X [cm]", size='x-large')
    ax.set_ylabel("Y [cm]", size='x-large')
    ax.set_zlabel("Energy [MeV]", size='x-large')
    fig.suptitle(title, size='xx-large')

    ax.view_init(elev=15, azim=-60, roll=0)


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
        plt.savefig(join(save_dir, output), bbox_inches='tight')
        print("Figure ", join(save_dir, output), " successfully created!")

    if show:
        plt.show()
