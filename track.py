import trackpy as tp
import pandas as pd
import matplotlib.pyplot as plt

def find_particle_locations(image, diameter=[21, 21], max_iterations=10, minmass=1, separation=15, save_path="particles.csv"):
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


# test the function
# image = plt.imread("E:\\shake_table_data\\N48\\20hz_hopperflow\\60deg\\10cm\\1\\1_00002.tif")
# particles = find_particle_locations(image, diameter=[5, 5], max_iterations=10, minmass=1, separation=5)

# # Display the located particles
# plt.figure()
# tp.annotate(particles, image, plot_style={'markersize': 1.3})
# plt.show()