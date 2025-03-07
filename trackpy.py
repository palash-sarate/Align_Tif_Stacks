import trackpy as tp
import pandas as pd
import matplotlib.pyplot as plt

def find_particle_locations(image, diameter, minmass=100):
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
    # Locate particles in the image
    particles = tp.locate(image, diameter, minmass=minmass)
    
    return particles

# Example usage:
# image = ...  # Load your image here
# diameter = 11  # Example particle diameter
# particles = find_particle_locations(image, diameter)
# print(particles)
# tp.annotate(particles, image)
# plt.imshow(image)
# plt.show()
