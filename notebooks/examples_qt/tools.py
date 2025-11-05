import numpy as np
import matplotlib.pyplot as plt
import qutip as qt


# calculated sorted eigenvalues
def eigensorter(H, threshold=1e-12):
    # calculate eigenvalues and eigenvectorss
    evals, evecs = np.linalg.eig(H)
    # get the list of sorted indices from the eigenvalues
    ids = np.argsort(evals)
    # sort the eigenvalues
    evals = evals[ids]
    # sort the eigenvectors
    evecs = evecs[:, ids]
    # transpose the vectors
    evecs = evecs.T
    # test the result
    error = np.linalg.norm(
        [evals[k] - vec.conj().T @ H @ vec for k, vec in enumerate(evecs)]
    )
    if error > threshold:
        # inform if error is too large
        print("Error size: ", error)
        return error
    else:
        # otherwise return results
        return evals, evecs


# import tqdm to monitor progress if the package is installed
def progressbar(iterator):
    # check if the package exists
    try:
        # import progress bar method
        from tqdm import tqdm

        # return progress bar iterator
        return tqdm(iterator)
    except ModuleNotFoundError:
        # otherwise leave iterator unchanged
        return iterator


# plot sphere
def plt_sphere(ax, radius=1, color="black", alpha=0.05):
    # draw sphere
    u, v = np.mgrid[0 : 2 * np.pi : 50j, 0 : np.pi : 50j]
    x = radius * np.cos(u) * np.sin(v)
    y = radius * np.sin(u) * np.sin(v)
    z = radius * np.cos(v)
    ax.plot_surface(x, y, z, color=color, alpha=alpha)
    u = np.linspace(0, 2 * np.pi, 100)
    x = radius * np.cos(u) * np.sin(-np.pi / 2)
    y = radius * np.sin(u) * np.sin(-np.pi / 2)
    z = radius * np.cos(-np.pi / 2 * np.ones(len(u)))
    ax.plot(x, y, z, lw=1, color=color, alpha=0.2)
    ax.plot(z, y, x, lw=1, color=color, alpha=0.2)
    R, N, N = np.linspace(-radius, radius, 100), np.zeros(100), np.zeros(100)
    ax.plot(R, N, N, lw=1, c="r", alpha=0.2)
    ax.plot(N, R, N, lw=1, c="g", alpha=0.2)
    ax.plot(N, N, R, lw=1, c="b", alpha=0.2)
    return


def plot_bloch_trajectories(
    states_list,
    labels=None,
    colors=None,
    linestyles=None,
    elev=15,
    azim=-45,
    offset=2.25,
    xy_offset_factor=0.6,
    proj_alpha=0.7,
    sphere_alpha=0.1,
    skip_density=200,
    fig_size=(10, 9),
):
    """Plot multiple Bloch trajectories on a single Bloch sphere with projections onto
    the three orthogonal planes.

    Parameters
    - states_list: list of list-of-states (each inner list is a time-ordered sequence of Qobj states)
    - labels: list of labels for each trajectory
    - colors: list of color strings
    - linestyles: list of linestyles
    - elev, azim: view angles
    - offset: distance of the side planes from the origin
    - xy_offset_factor: multiplier for the x-y plane offset (so it can be closer)
    - proj_alpha: alpha for projection lines/scatter
    - sphere_alpha: alpha for the Bloch sphere surface
    - skip_density: approximate number of points to display (controls subsampling)
    - fig_size: figure size tuple

    Returns (fig, ax)
    """
    sx = qt.sigmax()
    sy = qt.sigmay()
    sz = qt.sigmaz()

    # default labels/colors/linestyles
    n = len(states_list)
    if labels is None:
        labels = [f"Traj {i+1}" for i in range(n)]
    if colors is None:
        default_colors = [
            "tab:blue",
            "tab:green",
            "tab:orange",
            "tab:red",
            "tab:purple",
        ]
        colors = [default_colors[i % len(default_colors)] for i in range(n)]
    if linestyles is None:
        linestyles = ["-"] * n

    # compute Bloch vectors for each states list
    all_bvs = []
    for states in states_list:
        bv = np.real(
            np.array(
                [[qt.expect(sx, s), qt.expect(sy, s), qt.expect(sz, s)] for s in states]
            )
        )
        all_bvs.append(bv)

    # plotting
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    fig = plt.figure(figsize=fig_size)
    ax = fig.add_subplot(111, projection="3d")

    # draw Bloch sphere with QuTiP helper
    b = qt.Bloch(fig=fig, axes=ax)
    b.render()

    # set sphere alpha (try to only affect polygon/collection artists)
    for coll in list(ax.collections):
        try:
            coll.set_alpha(sphere_alpha)
        except Exception:
            pass

    def bv_to_xyz(bv):
        bv = np.asarray(bv)
        x = -bv[:, 1]
        y = bv[:, 0]
        z = bv[:, 2]
        return x, y, z

    # set view
    ax.view_init(elev=elev, azim=azim)
    az = np.deg2rad(azim)
    el = np.deg2rad(elev)
    cam_dir = np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])

    # compute plane offsets (x-y plane may be closer)
    z_plane = -np.sign(cam_dir[2]) * offset * xy_offset_factor
    x_plane = -np.sign(cam_dir[0]) * offset
    y_plane = -np.sign(cam_dir[1]) * offset

    # projection sampling
    skip = max(1, int(len(all_bvs[0]) / skip_density)) if len(all_bvs[0]) > 0 else 1

    main_handles = []
    for bv, color, ls, label in zip(all_bvs, colors, linestyles, labels):
        x, y, z = bv_to_xyz(bv)
        # main trajectory
        line = ax.plot(x, y, z, color=color, linestyle=ls, linewidth=2, label=label)[0]
        main_handles.append(line)

        # projections (x-y on z_plane)
        ax.plot(
            x,
            y,
            np.full_like(x, z_plane),
            color=color,
            linestyle=ls,
            alpha=proj_alpha,
            linewidth=1,
        )
        ax.scatter(
            x[::skip],
            y[::skip],
            np.full_like(x[::skip], z_plane),
            color=color,
            alpha=proj_alpha,
            s=6,
        )

        # x-z on fixed y plane
        ax.plot(
            x,
            np.full_like(x, y_plane),
            z,
            color=color,
            linestyle=ls,
            alpha=proj_alpha,
            linewidth=1,
        )
        ax.scatter(
            x[::skip],
            np.full_like(x[::skip], y_plane),
            z[::skip],
            color=color,
            alpha=proj_alpha,
            s=6,
        )

        # y-z on fixed x plane
        ax.plot(
            np.full_like(y, x_plane),
            y,
            z,
            color=color,
            linestyle=ls,
            alpha=proj_alpha,
            linewidth=1,
        )
        ax.scatter(
            np.full_like(y[::skip], x_plane),
            y[::skip],
            z[::skip],
            color=color,
            alpha=proj_alpha,
            s=6,
        )

        # start markers
        if len(x) > 0:
            ax.scatter(x[0], y[0], z[0], color=color, s=36)
            ax.scatter(x[0], y[0], z_plane, color=color, s=12, alpha=0.9)

    # draw projected circle borders (unit circles in plotting coords)
    T = np.linspace(0, 2 * np.pi, 300)
    cx = np.cos(T)
    cy = np.sin(T)
    ax.plot(
        cx,
        cy,
        np.full_like(cx, z_plane),
        color="gray",
        linewidth=1.25,
        alpha=proj_alpha,
    )
    ax.plot(
        cx,
        np.full_like(cx, y_plane),
        np.sin(T),
        color="gray",
        linewidth=1.25,
        alpha=proj_alpha,
    )
    ax.plot(
        np.full_like(cx, x_plane),
        cx,
        np.sin(T),
        color="gray",
        linewidth=1.25,
        alpha=proj_alpha,
    )

    # appearance
    ax.set_xlim(-1.3, 1.5)
    ax.set_ylim(-1.3, 1.5)
    ax.set_zlim(-1.3, 1.3)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    # legend only for main trajectories
    ax.legend(handles=main_handles, loc="upper left", bbox_to_anchor=(0.01, 0.98))
    plt.tight_layout()

    return fig, ax


def is_psd_dm(rho: qt.Qobj, tol=1e-12):
    """
    Return True iff rho is (numerically) a valid density matrix and PSD.
    tol is the allowed negative slack due to roundoff.
    """
    # Basic sanity checks (you can relax these if needed)
    if rho.type == "ket" or rho.type == "bra":
        rho = qt.ket2dm(rho)
    # Use Hermitian eigen-solver (eigh) via QuTiP
    evals = rho.eigenenergies()  # sorted real eigenvalues for Hermitian input
    # Trace ~ 1 (optional; skip if you only care about PSD)
    tr_ok = abs((rho.tr()).real - 1.0) <= tol
    # Positive semidefinite: all eigenvalues >= -tol
    psd_ok = np.min(evals) >= -tol
    # Within unit trace and non-negative (optional extra: >=0 and <=1)
    return tr_ok and psd_ok
