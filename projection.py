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

WidthZ = nLayer * LayerThick

RdPu = plt.cm.RdPu
OrRd = plt.cm.OrRd


def read_file(fname: str, tree: str, event_index: int):
    with up.open(fname) as f:
        tree = f[tree]

        EventID = tree['EventNumber'].array(library='np')

        try:
            event = np.where(EventID == event_index)[0][0]
        except:
            print('Event ID does not exist in the ROOT file!')
            raise

        print(f'Entry ID: {event}')
        print(f'Event ID: {event_index}')

        Hit_X = tree['Hit_X'].array(library='np')[event]
        Hit_Y = tree['Hit_Y'].array(library='np')[event]
        Hit_Z = tree['Hit_Z'].array(library='np')[event]
        Hit_Energy = tree['Hit_Energy'].array(library='np')[event]

        assert len(Hit_X) == len(Hit_Y) == len(Hit_Z) == len(Hit_Energy)

        return Hit_X, Hit_Y, Hit_Z, Hit_Energy


def plot_axis(fname: str, tree: str, run_index: int, event_index: int, title: str, staggered: bool):
    Hit_X, Hit_Y, Hit_Z, Hit_Energy = read_file(fname, tree, event_index)

    d_x = dict()
    d_y = dict()
    d_z = dict()

    for i in np.arange(len(Hit_X)):
        if Hit_X[i] not in d_x:
            d_x[Hit_X[i]] = Hit_Energy[i]
        else:
            d_x[Hit_X[i]] += Hit_Energy[i]

        if Hit_Y[i] not in d_y:
            d_y[Hit_Y[i]] = Hit_Energy[i]
        else:
            d_y[Hit_Y[i]] += Hit_Energy[i]

        if Hit_Z[i] not in d_z:
            d_z[Hit_Z[i]] = Hit_Energy[i]
        else:
            d_z[Hit_Z[i]] += Hit_Energy[i]

    WidthX = (nCellX + 0.5 * staggered) * CellWidthX
    WidthY = (nCellY + 0.5 * staggered) * CellWidthY
    binx = (staggered + 1) * nCellX
    biny = (staggered + 1) * nCellY

    x = np.array(list(d_x.keys()))
    energy_x = np.array(list(d_x.values()))

    y = np.array(list(d_y.keys()))
    energy_y = np.array(list(d_y.values()))

    z = np.array(list(d_z.keys()))
    energy_z = np.array(list(d_z.values()))

    fig = plt.figure()
    spec = fig.add_gridspec(2, 2)

    ax_z = fig.add_subplot(spec[:, 0])
    ax_x = fig.add_subplot(spec[0, 1])
    ax_y = fig.add_subplot(spec[1, 1])

    # Projection on z axis
    n_z, bins_z, patches_z = ax_z.hist(z, bins=nLayer, range=[0, WidthZ], weights=energy_z)
    height = n_z / n_z.max()
    for i in np.arange(len(patches_z)):
        patches_z[i].set_facecolor(RdPu((4 * height[i] + 1) * 0.2))
    ax_z.grid(ls=':')
    ax_z.set_axisbelow(True)
    ax_z.set_xlim(0, WidthZ)
    ax_z.set_xlabel("Z [cm]", size='large')
    ax_z.set_ylabel("Energy [MeV]", size='large')

    # Projection on x axis
    n_x, bins_x, patches_x = ax_x.hist(x, bins=binx, range=[-0.5 * WidthX, 0.5 * WidthX], weights=energy_x)
    height = n_x / n_x.max()
    for i in np.arange(len(patches_x)):
        patches_x[i].set_facecolor(OrRd((4 * height[i] + 1) * 0.2))
    ax_x.grid(ls=':')
    ax_x.set_axisbelow(True)
    ax_x.set_xlim(-0.5 * WidthX, 0.5 * WidthX)
    ax_x.set_xticks(np.linspace(-20, 20, 5))
    ax_x.set_xlabel("X [cm]", size='large')
    ax_x.set_ylabel("Energy [MeV]", size='large')

    # Projection on y axis
    n_y, bins_y, patches_y = ax_y.hist(y, bins=biny, range=[-0.5 * WidthY, 0.5 * WidthY], weights=energy_y)
    height = n_y / n_y.max()
    for i in np.arange(len(patches_y)):
        patches_y[i].set_facecolor(OrRd((4 * height[i] + 1) * 0.2))
    ax_y.grid(ls=':')
    ax_y.set_axisbelow(True)
    ax_y.set_xlim(-0.5 * WidthY, 0.5 * WidthY)
    ax_y.set_xticks(np.linspace(-20, 20, 5))
    ax_y.set_xlabel("Y [cm]", size='large')
    ax_y.set_ylabel("Energy [MeV]", size='large')

    fig.suptitle(title, size='xx-large')

    fig.text(0.05, 0.9, f'Run ID: {run_index}\nEvent ID: {event_index}')
    fig.text(0.75, 0.9, r'$E_\mathrm{total} =$' + f'{np.sum(Hit_Energy):.3f} MeV\n' + r'$E_\mathrm{max} =$' + f'{np.max(Hit_Energy):.3f} MeV')

    plt.tight_layout()


def plot_plane(fname: str, tree: str, event_index: int, title: str, staggered: bool):
    Hit_X, Hit_Y, Hit_Z, Hit_Energy = read_file(fname, tree, event_index)

    d_xy = dict()
    d_xz = dict()
    d_yz = dict()

    for i in np.arange(len(Hit_X)):
        if (Hit_X[i], Hit_Y[i]) not in d_xy:
            d_xy[Hit_X[i], Hit_Y[i]] = Hit_Energy[i]
        else:
            d_xy[Hit_X[i], Hit_Y[i]] += Hit_Energy[i]

        if (Hit_X[i], Hit_Z[i]) not in d_xz:
            d_xz[Hit_X[i], Hit_Z[i]] = Hit_Energy[i]
        else:
            d_xz[Hit_X[i], Hit_Z[i]] += Hit_Energy[i]

        if (Hit_Y[i], Hit_Z[i]) not in d_yz:
            d_yz[Hit_Y[i], Hit_Z[i]] = Hit_Energy[i]
        else:
            d_yz[Hit_Y[i], Hit_Z[i]] += Hit_Energy[i]

    x_xy = np.array(list(d_xy.keys()))[:, 0]
    y_xy = np.array(list(d_xy.keys()))[:, 1]
    energy_xy = np.array(list(d_xy.values()))

    x_xz = np.array(list(d_xz.keys()))[:, 0]
    z_xz = np.array(list(d_xz.keys()))[:, 1]
    energy_xz = np.array(list(d_xz.values()))

    y_yz = np.array(list(d_yz.keys()))[:, 0]
    z_yz = np.array(list(d_yz.keys()))[:, 1]
    energy_yz = np.array(list(d_yz.values()))

    WidthX = (nCellX + 0.5 * staggered) * CellWidthX
    WidthY = (nCellY + 0.5 * staggered) * CellWidthY

    ratioX_xy = WidthX / WidthY
    ratioX_xz = WidthX / WidthZ
    ratioX_yz = WidthY / WidthZ
    ratioY = 1
    ratioZ = 1

    fig = plt.figure()
    spec = fig.add_gridspec(2, 2)

    ax_xy = fig.add_subplot(spec[:, 0], projection='3d')
    plt.gca().set_box_aspect((ratioX_xy, ratioY, ratioZ))
    ax_xz = fig.add_subplot(spec[0, 1], projection='3d')
    plt.gca().set_box_aspect((ratioX_xz, ratioY, 0.85 * ratioZ))
    ax_yz = fig.add_subplot(spec[1, 1], projection='3d')
    plt.gca().set_box_aspect((ratioX_yz, ratioY, 0.85 * ratioZ))

    binx = (staggered + 1) * nCellX
    biny = (staggered + 1) * nCellY

    # Projection on xOy plane
    hist_xy, xedges_xy, yedges_xy = np.histogram2d(x_xy, y_xy, bins=[binx, biny], range=[[-0.5 * WidthX, 0.5 * WidthX], [-0.5 * WidthY, 0.5 * WidthY]],  weights=energy_xy)
    xpos_xy, ypos_xy = np.meshgrid(xedges_xy[:-1], yedges_xy[:-1], indexing='ij')
    xpos_xy = xpos_xy.ravel()
    ypos_xy = ypos_xy.ravel()
    zpos_xy = 0

    dx_xy = CellWidthX / (staggered + 1)
    dy_xy = CellWidthY / (staggered + 1)
    dz_xy = hist_xy.ravel()

    max_height_xy = np.max(dz_xy)
    min_height_xy = np.min(dz_xy)
    rgba = [RdPu((k - min_height_xy) / (max_height_xy - min_height_xy)) for k in dz_xy]

    ax_xy.bar3d(xpos_xy, ypos_xy, zpos_xy, dx_xy, dy_xy, dz_xy, color=rgba, edgecolor='none')

    clb = plt.cm.ScalarMappable(cmap=RdPu)
    clb.set_array(dz_xy)
    fig.colorbar(clb, pad=0.2, ax=ax_xy, orientation='horizontal')

    ax_xy.xaxis._axinfo['grid'].update(linestyle=':', linewidth=0.5)
    ax_xy.yaxis._axinfo['grid'].update(linestyle=':', linewidth=0.5)
    ax_xy.zaxis._axinfo['grid'].update(linestyle=':', linewidth=0.5)
    ax_xy.set_xlim(-0.5 * WidthX, 0.5 * WidthX)
    ax_xy.set_ylim(-0.5 * WidthY, 0.5 * WidthY)
    ax_xy.set_xlabel("X [cm]", size='x-large')
    ax_xy.set_ylabel("Y [cm]", size='x-large')
    ax_xy.set_zlabel("Energy [MeV]", size='x-large')
    ax_xy.view_init(elev=15, azim=-60, roll=0)

    # Projection on xOz plane
    hist_xz, xedges_xz, yedges_xz = np.histogram2d(x_xz, z_xz, bins=[binx, nLayer], range=[[-0.5 * WidthX, 0.5 * WidthX], [0, WidthZ]],  weights=energy_xz)
    xpos_xz, ypos_xz = np.meshgrid(xedges_xz[:-1], yedges_xz[:-1], indexing='ij')
    xpos_xz = xpos_xz.ravel()
    ypos_xz = ypos_xz.ravel()
    zpos_xz = 0

    dx_xz = CellWidthX / (staggered + 1)
    dy_xz = LayerThick
    dz_xz = hist_xz.ravel()

    max_height_xz = np.max(dz_xz)
    min_height_xz = np.min(dz_xz)
    rgba = [OrRd((k - min_height_xz) / (max_height_xz - min_height_xz)) for k in dz_xz]

    ax_xz.bar3d(xpos_xz, ypos_xz, zpos_xz, dx_xz, dy_xz, dz_xz, color=rgba, edgecolor='none')

    clb = plt.cm.ScalarMappable(cmap=OrRd)
    clb.set_array(dz_xz)

    ax_xz.xaxis._axinfo['grid'].update(linestyle=':', linewidth=0.5)
    ax_xz.yaxis._axinfo['grid'].update(linestyle=':', linewidth=0.5)
    ax_xz.zaxis._axinfo['grid'].update(linestyle=':', linewidth=0.5)
    ax_xz.set_xlim(-0.5 * WidthX, 0.5 * WidthX)
    ax_xz.set_ylim(0, WidthZ)
    ax_xz.set_xlabel("X [cm]", size='large')
    ax_xz.set_ylabel("Z [cm]", size='large')
    ax_xz.set_zlabel("Energy [MeV]", size='large')
    ax_xz.view_init(elev=10, azim=-40, roll=0)

    # Projection on yOz plane
    hist_yz, xedges_yz, yedges_yz = np.histogram2d(y_yz, z_yz, bins=[biny, nLayer], range=[[-0.5 * WidthY, 0.5 * WidthY], [0, WidthZ]],  weights=energy_yz)
    xpos_yz, ypos_yz = np.meshgrid(xedges_yz[:-1], yedges_yz[:-1], indexing='ij')
    xpos_yz = xpos_yz.ravel()
    ypos_yz = ypos_yz.ravel()
    zpos_yz = 0

    dx_yz = CellWidthY / (staggered + 1)
    dy_yz = LayerThick
    dz_yz = hist_yz.ravel()

    max_height_yz = np.max(dz_yz)
    min_height_yz = np.min(dz_yz)
    rgba = [OrRd((k - min_height_yz) / (max_height_yz - min_height_yz)) for k in dz_yz]

    ax_yz.bar3d(xpos_yz, ypos_yz, zpos_yz, dx_yz, dy_yz, dz_yz, color=rgba, edgecolor='none')

    clb = plt.cm.ScalarMappable(cmap=OrRd)
    clb.set_array(dz_yz)

    ax_yz.xaxis._axinfo['grid'].update(linestyle=':', linewidth=0.5)
    ax_yz.yaxis._axinfo['grid'].update(linestyle=':', linewidth=0.5)
    ax_yz.zaxis._axinfo['grid'].update(linestyle=':', linewidth=0.5)
    ax_yz.set_xlim(-0.5 * WidthY, 0.5 * WidthY)
    ax_yz.set_ylim(0, WidthZ)
    ax_yz.set_xlabel("Y [cm]", size='large')
    ax_yz.set_ylabel("Z [cm]", size='large')
    ax_yz.set_zlabel("Energy [MeV]", size='large')
    ax_yz.view_init(elev=10, azim=-40, roll=0)

    fig.suptitle(title, size='xx-large')

    fig.text(0.05, 0.85, f'Run ID: {run_index}\nEvent ID: {event_index}')
    fig.text(0.75, 0.85, r'$E_\mathrm{total} =$' + f'{np.sum(Hit_Energy):.3f} MeV\n' + r'$E_\mathrm{max} =$' + f'{np.max(Hit_Energy):.3f} MeV')

    plt.tight_layout()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--file",      type=str, default='',   required=True,  help="Input ROOT file")
    parser.add_argument("-t", "--tree",      type=str, default='dp', help="Input tree name (default: dp)")
    parser.add_argument("-g", "--staggered", type=int, default=1,    choices=[0, 1], help="Staggered structure")
    parser.add_argument("-i", "--title",     type=str, default='',   help="Title of display figure")
    parser.add_argument("-r", "--run",       type=int, default=0,    help="Run ID of the event to be displayed")
    parser.add_argument("-e", "--event",     type=int, default=0,    help="Event ID")
    parser.add_argument("-d", "--dir",       type=str, default=None, help="Directory to save the plot")
    parser.add_argument("-o", "--output",    type=str, default=None, help="Output file name")
    parser.add_argument("-s", "--show",      type=int, default=1,    choices=[0, 1], help="Instantly display or not")
    args = parser.parse_args()

    filename = args.file
    tree = args.tree
    staggered = args.staggered
    title = args.title
    run_index = args.run
    event_index = args.event
    save_dir = args.dir
    output = args.output
    show = args.show

    plot_axis(filename, tree, run_index, event_index, title, staggered)

    if save_dir and output:
        plt.savefig(join(save_dir, "Axis" + output))
        print("Figure ", join(save_dir, "Axis" + output), " successfully created!")

    plot_plane(filename, tree, event_index, title, staggered)

    if save_dir and output:
        plt.savefig(join(save_dir, "Plane" + output))
        print("Figure ", join(save_dir, "Plane" + output), " successfully created!")

    if show:
        plt.show()
