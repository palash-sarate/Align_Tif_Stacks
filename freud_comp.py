#%%
import freud

box, points = freud.data.make_random_system(10, 100, seed=0)
ql = freud.order.Steinhardt(l=6)
ql.compute((box, points), {'r_max':3})

# plot the ql values
import matplotlib.pyplot as plt
plt.plot(ql.particle_order)
plt.show()

#%%
import freud.box as Box
# lengths = []
# angles = []
# box = Box.from_box_lengths_and_angles(lengths, angles)
# particles_path = "E:\\shake_table_data\\N4\\4hz_hopperflow\\60deg\\10cm\\particle_locations_1.csv"
# particles = pd.read_csv(particles_path)
# points = particles[['x', 'y']].values
# ql = freud.order.Steinhardt(l=6)
# ql.compute((box, points), {'r_max':3})