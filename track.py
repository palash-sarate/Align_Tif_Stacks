#%%
import trackpy as tp
import pandas as pd
import matplotlib.pyplot as plt
import freud
import numpy as np


def find_particle_locations(image_path, diam=21, max_iterations=10, minmass=1, separation=15, save_path="particles.csv"):
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
        
    return True
    # return particles

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
               bead_dia = 7.1789, l=6, nhood = 6):
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
    # Create a 2D freud box object
    box = freud.box.Box.square(L=800)
    system = (box, points)
    # Create a freud Ql object
    Ql = freud.order.Steinhardt(l)

    # Compute the RDF
    ql = Ql.compute(system, neighbors={'r_max': r_max})
    psi = ql.particle_order
    return psi
    
#%%
# import matplotlib.pyplot as plt
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
