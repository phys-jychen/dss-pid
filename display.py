import uproot as up
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import argparse
from os.path import join

CellWidthX = 2.5
CellWidthY = 2.5
LayerThick = 4
nCellX = 21
nCellY = 21
nLayer = 11


def read_file(fname: str, tree: str, event_index: int, staggered: bool):
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

        xtemp = tree['Hit_X'].array(library='np')[event]
        ytemp = tree['Hit_Y'].array(library='np')[event]
        x = np.zeros_like(xtemp)
        y = np.zeros_like(ytemp)
        z = np.round(tree['Hit_Z'].array(library='np')[event] / LayerThick).astype(int)
        energy = tree['Hit_Energy'].array(library='np')[event]

        assert len(x) == len(y) == len(z) == len(energy)

        if staggered:
            for i in np.arange(len(z)):
                if z[i] % 2 == 0:
                    x[i] = np.round(xtemp[i] / CellWidthX - 0.25).astype(int)
                    y[i] = np.round(ytemp[i] / CellWidthY - 0.25).astype(int)
                else:
                    x[i] = np.round(xtemp[i] / CellWidthX + 0.25).astype(int)
                    y[i] = np.round(ytemp[i] / CellWidthY + 0.25).astype(int)
        else:
            x[i] = np.round(xtemp[i] / CellWidthX).astype(int)
            y[i] = np.round(ytemp[i] / CellWidthY).astype(int)

        znew, ynew, xnew, enew = (np.array(a) for a in zip(*sorted(zip(z, y, x, energy), reverse=True)))

        return xnew, ynew, znew, enew


def plot(fname: str, tree: str, run_index: int, event_index: int, title: str, staggered: bool):
    x, y, z, energy = read_file(fname, tree, event_index, staggered)
    energy_norm = energy / np.max(energy)

    nhits = len(x)

    WidthX = (nCellX + 0.5 * staggered) * CellWidthX
    WidthY = (nCellY + 0.5 * staggered) * CellWidthY

    ratioX = WidthX / (LayerThick * (nLayer + 1))
    ratioY = 1
    ratioZ = WidthY / (LayerThick * (nLayer + 1))

    fig = plt.figure()
    spec = fig.add_gridspec(2, 2)

    ax_xz = fig.add_subplot(spec[0, 1], projection='3d')
    plt.gca().set_box_aspect((ratioX, ratioY, ratioZ))
    ax_xy = fig.add_subplot(spec[1, 1], projection='3d')
    plt.gca().set_box_aspect((ratioX, ratioY, ratioZ))
    ax = fig.add_subplot(spec[:, 0], projection='3d')
    plt.gca().set_box_aspect((ratioX, ratioY, ratioZ))

    cmap = cm.OrRd

    for i in np.arange(nhits):
        if z[i] % 2 == 0:
            xnew = np.arange(x[i] - 1, x[i] + 1) + 0.75 * staggered
            ynew = np.arange(y[i] - 1, y[i] + 1) + 0.75 * staggered
        else:
            xnew = np.arange(x[i] - 1, x[i] + 1) + 0.25 * staggered
            ynew = np.arange(y[i] - 1, y[i] + 1) + 0.25 * staggered

        xnew, ynew = np.meshgrid(xnew, ynew)
        znew = z[i] * np.ones(xnew.shape)
        enew = energy_norm[i] * np.ones(xnew.shape)

        for axis in (ax, ax_xz, ax_xy):
            axis.plot_surface(xnew, znew, ynew, cmap=cmap, facecolors=cmap(enew), edgecolor='k', alpha=0.8, lw=0.05, rstride=1, cstride=1, antialiased=False)
            ax.computed_zorder = False

    for axis in (ax, ax_xz, ax_xy):
        axis.set_xlim(-0.5 * nCellX, 0.5 * nCellY)
        axis.set_ylim(0, nLayer)
        axis.set_zlim(-0.5 * nCellY, 0.5 * nCellY)
        axis.set_aspect(aspect='equalxz')
        axis.grid(False)

        axis.invert_xaxis()

        if axis == ax:
            axis.set_xticks(np.linspace(-10, 10, 5), CellWidthX * np.linspace(-10, 10, 5))
            axis.set_yticks(np.linspace(0, nLayer + 1, 5))
            axis.set_zticks(np.linspace(-10, 10, 5), CellWidthY * np.linspace(-10, 10, 5))
            axis.set_xlabel("X [cm]", size='x-large')
            axis.set_ylabel("Z [layer]", size='x-large')
            axis.set_zlabel("Y [cm]", size='x-large')

            m = plt.cm.ScalarMappable(cmap=cmap)
            m.set_array(energy)
            plt.colorbar(m, pad=0.2, ax=plt.gca(), orientation='horizontal').set_label(label="Hit Energy [MeV]", size='large')
        else:
            axis.set_xticks([])
            axis.set_yticks([])
            axis.set_zticks([])

    ax_xz.set_xlabel("X", size='large', labelpad=-10)
    ax_xz.set_ylabel("Z", size='large', labelpad=-10)
    ax_xy.set_xlabel("X", size='large', labelpad=-10)
    ax_xy.set_zlabel("Y", size='large', labelpad=-10)

    fig.suptitle(title, size='xx-large')

    ax.text2D(0.05, 0.95, f'Run ID: {run_index}\nEvent ID: {event_index}', transform=ax.transAxes)
    ax.text2D(0.8, 0.95, r'$E_\mathrm{total} =$' + f'{np.sum(energy):.3f} MeV\n' + r'$E_\mathrm{max} =$' + f'{np.max(energy):.3f} MeV', transform=ax.transAxes)

    ax.view_init(elev=20, azim=-35, roll=0)
    ax_xz.view_init(elev=90, azim=0, roll=0)
    ax_xy.view_init(elev=0, azim=-90, roll=0)


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

    plot(filename, tree, run_index, event_index, title, staggered)

    if save_dir and output:
        plt.savefig(join(save_dir, output), bbox_inches='tight')
        print("Figure ", join(save_dir, output), " successfully created!")

    if show:
        plt.show()
