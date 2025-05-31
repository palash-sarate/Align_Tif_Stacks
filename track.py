#%%
import trackpy as tp
import pandas as pd
import matplotlib.pyplot as plt
import freud
import numpy as np


def find_particle_locations(image_path, diam=5, max_iterations=10, minmass=800, separation=5, save_path="particles.csv",overlay=False):
    """
    Find particle locations in an image using trackpy.

    Parameters:
    image : ndarray
        The image in which to find particles.
    diameter : int
        The expected diameter of the particles.
    minmass : float, optional
        The minimum integrated brightness for a particle to be considered valid.

    Returns:
    DataFrame
        A DataFrame containing the coordinates and properties of the located particles.
    """
    diameter = [diam, diam]
    image = plt.imread(image_path)
    # preprocess the image (threshold, denoise, etc.) if necessary
    # threshold = 30  # adjust based on your images
    # image[image < threshold] = 0
    
    # print(image.shape)
    # Locate particles (adjust diameter based on your particles' size)
    particles = tp.locate(
        image,
        diameter=diameter,
        max_iterations=max_iterations,
        minmass=minmass,
        separation=separation,
    )
    particles_csv = particles.to_csv(index=False)
    # save the csv file
    with open(save_path, "w") as f:
        f.write(particles_csv)
    
    if not overlay:
        return True
    
    # overlay the particles on the image
    fig, ax = plt.subplots()
    tp.annotate(particles, image, ax=ax, plot_style={'markersize': 0.3})
    ax.set_xlim(300, 800)
    ax.set_ylim(200, 800)
    plt.show()

def calculate_rdf(particles_path, r_max=100, dr=1):
    """
    Calculate the radial distribution function (RDF) for a set of particles.

    Parameters:
    particles : DataFrame
        A DataFrame containing the coordinates of the particles.
    r_max : int, optional
        The maximum radius to consider.
    dr : float, optional
        The bin size for the radii.

    Returns:
    Series
        A Series containing the RDF.
    """
    # load particles from csv
    particles = pd.read_csv(particles_path)
    # Calculate the RDF
    rdf = tp.pair_correlation_2d(particles, cutoff=r_max, dr=dr)
    return rdf

def compute_ql(particles_path= "F:\\shake_table_data\\N4\\4hz_hopperflow\\60deg\\10cm\\particle_locations_1.csv",
               bead_dia = 7.1789, l=6, nhood = 2, box_size=800):
    # TODO : Add padding to the box for border corrections as the si 6 should have a mean of 1
    """
    Compute the local order parameter QL for a set of particles.

    Parameters:
    particles : DataFrame
        A DataFrame containing the coordinates of the particles.

    Returns:
    Series
        A Series containing the QL values.
    """
 
    r_max = bead_dia * nhood
    # load particles from csv
    particles = pd.read_csv(particles_path)
    points = particles[['x', 'y']].values
    # Add a z-coordinate of 0 to make the points 3D
    points = np.concatenate((points, np.zeros((points.shape[0], 1))), axis=1)
    # shift the points to the center of the box
    points[:, 0] -= box_size/2
    points[:, 1] -= box_size/2
    # Create a 2D freud box object
    box = freud.box.Box.square(L=box_size)
    system = (box, points)
    
    # ax = plt.figure().add_subplot()
    # system1 = freud.AABBQuery(box, points)
    # system1.plot(ax)
    # plt.show()
    
    # Create a freud Ql object
    Ql = freud.order.Steinhardt(l)

    # Compute the RDF
    ql = Ql.compute(system, neighbors={'r_max': r_max})
    psi = ql.particle_order
    return psi

# compute bond orientation order parameter
def compute_bond_orientation_order(particles_path= "F:\\shake_table_data\\N4\\4hz_hopperflow\\60deg\\10cm\\particle_locations_1.csv",
                                   l=6, num_neighbors=6):
    """
    Compute the bond orientation order parameter for a set of particles.

    Parameters:
    particles : DataFrame
        A DataFrame containing the coordinates of the particles.

    Returns:
    Series
        A Series containing the bond orientation order parameter values.
    """
    # load particles from csv
    particles = pd.read_csv(particles_path)
    points = particles[['x', 'y']].values
    # Add a z-coordinate of 0 to make the points 3D
    points = np.concatenate((points, np.zeros((points.shape[0], 1))), axis=1)
    
    order_parameter = np.zeros(len(points))
    
    for i in range(len(points)):
        # find the closest n neighbors of the particle
        distances = np.linalg.norm(points - points[i], axis=1)
        # keep only first n neighbors
        neighbors = np.argsort(distances)[1:num_neighbors+1]
        # compute the bond orientation order parameter for the particle
        ql = 0
        for j in neighbors:
            # angle between the line connecting the particle to its neighbors and the x-axis
            theta = np.arctan2(points[j][1] - points[i][1], points[j][0] - points[i][0])
            # compute the bond orientation order parameter for the particle
            ql += np.exp(1j * l * theta)
        ql /= num_neighbors
        # store the bond orientation order parameter for the particle
        order_parameter[i] = np.abs(ql)
        
    return order_parameter

#%% test find_particle_locations function
# minmass=800
# image_path = "F:\\shake_table_data\\N4\\4hz_hopperflow\\60deg\\10cm\\1\\1_00002.TIF"
# find_particle_locations(image_path, diam=5, max_iterations=10,
#                         minmass=minmass, separation=5, 
#                         save_path="particles.csv",overlay=True)
# image_path = "F:\\shake_table_data\\N4\\4hz_hopperflow\\60deg\\10cm\\1\\1_00968.TIF"
# find_particle_locations(image_path, diam=5, max_iterations=10,
#                         minmass=minmass, separation=5, 
#                         save_path="particles.csv",overlay=True)
#%% test compute_ql function
# psi = compute_ql()
# # plot histogram of psi
# plt.hist(psi, bins=100, density=True)
# plt.xlabel('QL')
# plt.ylabel('Density')
# plt.title('Distribution of QL')
# plt.xlim(0, 1)
# plt.show()
#%% test compute_bond_orientation_order function
# bond_orders = compute_bond_orientation_order(num_neighbors=12)
# # get the bond orientation order parameter
# bond_order = np.mean(bond_orders)
# # get the standard deviation of the bond orientation order parameter
# bond_order_std = np.std(bond_orders)
# print(f"Bond orientation order parameter: {bond_order:.3f} +/- {bond_order_std:.3f}")

# # plot histogram of bond orientation order parameter
# plt.hist(bond_orders, bins=100, density=True)
# plt.xlabel('Bond orientation order parameter')

# plt.ylabel('Density')
# plt.title('Distribution of Bond orientation order parameter')
# # plt.xlim(0, 1)
# plt.show()

#%% import matplotlib.pyplot as plt
# nhood_values = [4, 6, 8, 10 , 12, 14, 16]
# colors = ['red', 'blue', 'green', 'purple', 'orange', 'pink', 'brown']

# for nhood, color in zip(nhood_values, colors):
#     psi = compute_ql(nhood=nhood)
#     mean = np.mean(psi)
#     std_dev = np.std(psi)
#     print(f"nhood={nhood}: mean={mean:.3f}, std_dev={std_dev:.3f}")
#     plt.hist(psi, density=True, alpha=0.5, color=color, label=f'nhood={nhood}')

# plt.xlim(0, 1)
# plt.legend()
# plt.xlabel('QL')
# plt.ylabel('Density')
# plt.title('Distribution of QL for different nhood values')
# plt.show()
#%%
# sc = freud.data.UnitCell.sc()
# sc_system = sc.generate_system(5)
# test the function
# image = 'E:\\shake_table_data\\N4\\4hz_hopperflow\\60deg\\10cm\\2\\1_00001.TIF'
# particle_radius = 15
# particles_path = "E:\\shake_table_data\\N4\\4hz_hopperflow\\60deg\\10cm\\particle_locations_1.csv"
# rdf = calculate_rdf(particles_path, r_max=8*particle_radius, dr=particle_radius/2)
# # # Display the located particles
# plt.figure()
# # tp.annotate(particles, image, plot_style={'markersize': 1.3})
# rs = rdf[0]
# g = rdf[1]
# plt.plot(rs[1:]/particle_radius, g)
# plt.show()

# %%
