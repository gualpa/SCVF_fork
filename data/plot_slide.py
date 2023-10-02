#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection


def read_data_file(file_name, column_names):
    """
    Reads a data file with the following column format:
    - Column 1: Float
    - Column 2: Float
    - Column 3: Float
    - Column 4: Float
    - Column 5: Float
    - Column 6: Float
    - Column 7: Float

    Args:
        file_name (list): name of the file to read
        column_names (list): list of names of the columns of 
                             the dataframe to be returne.

    Returns:
        Dataframe: dataframe with the read data of the file

    Example:
        >>> columns = ['x', 'y', 'z']
        >>> file_name = 'file.dat'
        >>> read_data_file(file_name, columns)
                x      y      z
            0   1.1    2.3    8.9
            1   4.1    5.5    3.7
    """
    data = []
    with open(file_name, 'r') as file:
        for line in file:
            parts = line.strip().split()
            data_line = [float(parts[0])] + [float(x) for x in parts[1:]]
            data.append(data_line)
    df = pd.DataFrame(data, columns=column_names)
    return df


def read_pops(file_name):
    with open(file_name, 'r') as f:
        npop = int(f.readline())
        popcorn = []
        for _ in range(npop):
            head_pop = list(map(float, f.readline().split()))
            id = int(head_pop[0])
            nmem = int(head_pop[1])
            vol = head_pop[2]
            npart = int(head_pop[3])
            x = []
            y = []
            z = []
            r = []
            if nmem > 0:
                for _ in range(nmem):
                    popmem = list(map(float, f.readline().split()))
                    x.append(popmem[0])
                    y.append(popmem[1])
                    z.append(popmem[2])
                    r.append(popmem[3])
            if npart > 0:
                for _ in range(npart):
                    f.readline()
            pop = pd.DataFrame({'x': x, 'y': y, 'z': z, 'r': r})
            df = {'id': id, 'nmem': nmem, 'vol': vol, 'npart': npart, 'pop': pop}
            popcorn.append(df)
    return popcorn

def read_data_slide(hname, sphname, fname):
    #hname = "./halos_cut_scorpio.dat"
    #sphname = "./sphvoids_tol_0.00.dat"
    #fname = "./clean_popvds.dat"
    
    #sphd = pd.read_table(sphname, names=["id","rad","x","y","z","vx","vy","vz","delta","dummy1","dummy2"])
    sphd = read_data_file(sphname, ["id","rad","x","y","z","vx","vy","vz","delta","dummy1","dummy2"])

    #pts = pd.read_table(hname, names=["n","x","y","z","vx","vy","vz"])
    pts = read_data_file(hname, ["n","x","y","z","vx","vy","vz"])


    pops = read_pops(fname)
    return {'sphd': sphd, 'pts': pts, 'pops': pops}

def fcircle(xc, yc, rc, cl, bc):
    tt = np.linspace(0, 2*np.pi, 100)
    xt = rc*np.cos(tt) + xc
    yt = rc*np.sin(tt) + yc
    return Circle((xc, yc), rc, edgecolor=bc, facecolor=cl, linewidth=2)

def plot_slide(zslide, data, title, output_file):
    sphd = data['sphd'] # create_zero_hname_dataframe()#  circulos
    pts = data['pts'] # create_zero_sphname_dataframe() # estrellas
    pops = data['pops'] #create_zero_pops_dataframe() #data['pops'] # circulos de colores




    sphd_sorted = sphd.sort_values('z').reset_index()
    sph_rad = sphd_sorted['rad']
    sph_x = sphd_sorted['x']
    sph_y = sphd_sorted['y']
    sph_z = sphd_sorted['z']

    color_map = plt.cm.get_cmap('tab20')
    coco = color_map(np.linspace(0, 1, 30))
    
    x = []
    y = []
    z = []
    rad = []
    idv = []
    nmem = []
    
    for i in range(len(pops)):
        pp = pops[i]['pop']
        nn = len(pp)
        for j in range(nn):
            x.append(pp['x'][j])
            y.append(pp['y'][j])
            z.append(pp['z'][j])
            rad.append(pp['r'][j])
            idv.append(i)
            nmem.append(nn)
    
    mR = 20
    
    cut = ((np.array(z) > (zslide - mR)) & (np.array(z) < (zslide + mR)))
    x = np.array(x)[cut]
    y = np.array(y)[cut]
    z = np.array(z)[cut]
    rad = np.array(rad)[cut]
    idv = np.array(idv)[cut]
    nmem = np.array(nmem)[cut]
    
    fig, ax = plt.subplots()
    ax.set_xlim([0, 500])
    ax.set_ylim([0, 500])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title(f"Slide at Z={zslide} {title}")
    
    patches = []
    colors = []
    
    for id in np.unique(idv):
        im = np.where(idv == id)[0]
        #colo = coco[np.random.randint(0, 30)]
        colo = coco[np.random.randint(0, 30)]

        if len(im) >= 1 and nmem[im[0]] >= 1:
            for ii in im:
                rho = z[ii] - zslide
                if abs(rho) < rad[ii]:
                    radio = np.sqrt(rad[ii]**2 - rho**2)
                    circle = fcircle(x[ii], y[ii], radio, colo, colo)
                    patches.append(circle)
            colors.extend([colo] * len(im))
    
    collection = PatchCollection(patches, edgecolors=colors, facecolors=colors)
    ax.add_collection(collection)
    
    #print("len(sph_rad)",len(sph_rad))
    for ii in range(len(sph_rad)):


        rho = sph_z[ii] - zslide
        #print("rho ",rho," sph_z[ii]  ",sph_z[ii] ," zslide ", zslide)
        if abs(rho) < sph_rad[ii]:
            radio = np.sqrt(sph_rad[ii]**2 - rho**2)
            circle = fcircle(sph_x[ii], sph_y[ii], radio, 'none', 'black')
            ax.add_artist(circle)
    
    i = np.where(np.abs(pts['z'] - zslide) < 20.0)[0]
    ax.plot(pts['x'][i], pts['y'][i], '.', color='black', markersize=0.5)
    
    plt.savefig(output_file)
    plt.close()


def create_zero_hname_dataframe():
    """
    Creates a DataFrame with the specified format, setting all values to zero.

    Returns:
        pandas.DataFrame: DataFrame with zero values.
    """
    data = {
        'id': [0],
        'rad': [0.0],
        'x': [0.0],
        'y': [0.0],
        'z': [0.0],
        'vx': [0.0],
        'vy': [0.0],
        'vz': [0.0],
        'delta': [0.0],
        'dummy1': [0.0],
        'dummy2': [0.0]
    }
    df = pd.DataFrame(data)
    return df

def create_zero_sphname_dataframe():
    """
    Creates a DataFrame with the specified format (n, x, y, z, vx, vy, vz),
    setting all values to zero.

    Returns:
        pandas.DataFrame: DataFrame with zero values.
    """
    data = {
        'n': [0.0],
        'x': [0.0],
        'y': [0.0],
        'z': [0.0],
        'vx': [0.0],
        'vy': [0.0],
        'vz': [0.0]
    }
    df = pd.DataFrame(data)
    return df

def create_zero_pops_dataframe():
    """
    Creates a DataFrame with the specified format:
    [{'id': int, 'nmem': int, 'vol': float, 'npart': int, 'pop': DataFrame}]

    Returns:
        pandas.DataFrame: DataFrame with zero values.
    """
    data = {
        'id': [1],
        'nmem': [0],
        'vol': [0.0],
        'npart': [0],
        'pop': [pd.DataFrame({'x': [0.0], 'y': [0.0], 'z': [0.0], 'r': [0.0]})]
    }
    df = pd.DataFrame(data)
    return df

def create_zero_pops_dataframe():
    """
    Creates a list of dictionaries with the specified format:
    [{'id': int, 'nmem': int, 'vol': float, 'npart': int, 'pop': DataFrame}]

    Returns:
        list: List of dictionaries with the specified values.
    """
    data = [
        {
            'id': 8,
            'nmem': 4,
            'vol': 48137.82,
            'npart': 7,
            'pop': pd.DataFrame({
                'x': [0, 0, 0, 0],
                'y': [0, 0, 0, 0],
                'z': [0, 0, 0, 0],
                'r': [0, 0, 0, 0]
            })
        }
    ]
    return data


def plot_slide_3d(zslide, data, title, output_file):
    sphd = data['sphd']
    pts = data['pts']
    pops = data['pops']

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f"{title}")

    #for i in range(len(pts)):
    #    ax.scatter(pts['x'][i], pts['y'][i], pts['z'][i], c='black', s=0.1)
        #ax.scatter(pts['x'][i], pts['y'][i], pts['z'][i], c='black', s=1)

    sphd_sorted = sphd.sort_values('z').reset_index()
    sph_rad = sphd_sorted['rad']
    sph_x = sphd_sorted['x']
    sph_y = sphd_sorted['y']
    sph_z = sphd_sorted['z']

    for i in range(len(sph_rad)):
        rho = sph_z[i] - zslide
        if abs(rho) < sph_rad[i]:
            radio = np.sqrt(sph_rad[i] ** 2 - rho ** 2)
            u = np.linspace(0, 2 * np.pi, 100)
            v = np.linspace(0, np.pi, 100)
            x = sph_x[i] + radio * np.outer(np.cos(u), np.sin(v))
            y = sph_y[i] + radio * np.outer(np.sin(u), np.sin(v))
            z = sph_z[i] + radio * np.outer(np.ones(np.size(u)), np.cos(v))
            ax.plot_surface(x, y, z, color='none', edgecolor='black')


    for i in range(len(pops)):
        pp = pops[i]['pop']
        nn = len(pp)
        for j in range(nn):
            ax.scatter(pp['x'][j], pp['y'][j], pp['z'][j], c='red', s=2)
            #ax.scatter(pp['x'][j], pp['y'][j], pp['z'][j], c='red', s=5)
    plt.savefig(output_file)
    plt.close()



def plot_slide_3d_2(zslide, data, title, output_file):
    sphd = data['sphd']
    pts = data['pts']
    pops = data['pops']

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f"{title}")

    ax.scatter(pts['x'], pts['y'], pts['z'], c='black', s=0.1)

    sphd_sorted = sphd.sort_values('z').reset_index()
    sph_rad = sphd_sorted['rad']
    sph_x = sphd_sorted['x']
    sph_y = sphd_sorted['y']
    sph_z = sphd_sorted['z']

    rho_values = sph_z - zslide

    for i in range(len(sph_rad)):
        rho = rho_values[i]
        if abs(rho) < sph_rad[i]:
            radio = np.sqrt(sph_rad[i] ** 2 - rho ** 2)
            u = np.linspace(0, 2 * np.pi, 100)
            v = np.linspace(0, np.pi, 100)
            cos_u = np.cos(u)
            sin_v = np.sin(v)
            x = sph_x[i] + radio * np.outer(cos_u, sin_v)
            #y = sph_y[i] + radio * np.outer(sin_u, sin_v)
            y = sph_y[i] + radio * np.outer(np.sin(u), np.sin(v))

            z = sph_z[i] + radio * np.outer(np.ones(np.size(u)), np.cos(v))
            ax.plot_surface(x, y, z, color='none', edgecolor='black')

    for i in range(len(pops)):
        pp = pops[i]['pop']
        ax.scatter(pp['x'], pp['y'], pp['z'], c='red', s=2)

    plt.savefig(output_file)
    plt.close()

hname = "./halos_cut_scorpio.dat"
sphname = "./sphvoids_tol_0.00.dat"
fname = "./clean_popvds.dat"

data = read_data_slide(hname, sphname, fname)

#plot_slide(50.0, data, "MWE, z=50", "slide1.pdf")
#plot_slide(250.0, data, "MWE, z=250", "slide2.pdf")
#plot_slide(400.0, data, "MWE, z=400", "slide3.pdf")
#plot_slide_3d(250.0, data, "MWE, z=250", "slide2_3D.pdf")
plot_slide_3d(250.0, data, "MWE, z=250", "slide2_3D_4.pdf")

