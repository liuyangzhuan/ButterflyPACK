import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import label, generate_binary_structure

def generate_shapes(shape_background, shape, offset, num_shapes, slowness_min, slowness_max, slowness_back):
    volume = np.zeros(shape, dtype=bool)

    def add_cube(center, size):
        x_start, x_end = max(0, center[0] - size // 2), min(shape[0], center[0] + size // 2)
        y_start, y_end = max(0, center[1] - size // 2), min(shape[1], center[1] + size // 2)
        z_start, z_end = max(0, center[2] - size // 2), min(shape[2], center[2] + size // 2)
        volume[x_start:x_end, y_start:y_end, z_start:z_end] = True

    def add_ellipsoid(center, axes):
        x, y, z = np.ogrid[:shape[0], :shape[1], :shape[2]]
        axes = np.where(axes < 2, 2, axes)
        term = ((x - center[0])**2 / axes[0]**2) + ((y - center[1])**2 / axes[1]**2) + ((z - center[2])**2 / axes[2]**2)
        volume[term <= 1] = True

    for _ in range(num_shapes):
        center = np.random.randint(0, min(shape), 3)
        if np.random.rand() > 0.66:
            add_cube(center, np.random.randint(shape[0]*0.1, shape[0]*0.2))
        elif np.random.rand() > 0.33:
            add_ellipsoid(center, np.random.randint(shape[0]*0.05, shape[0]*0.15, 3))
        else:
            axes = np.random.randint(shape[0]*0.05, shape[0]*0.15, 3)
            axes[1] = axes[0] // 2
            add_ellipsoid(center, axes)

    labeled_volume, num_features = label(volume, structure=generate_binary_structure(3, 1))
    slowness_map=np.full_like(volume, slowness_back, dtype=np.float64)
    for label_num in range(1, num_features + 1):
        slowness_map[labeled_volume == label_num] = np.random.uniform(slowness_min, slowness_max)

    volume_background = np.zeros(shape_background, dtype=bool)
    slowness_map_tot = np.full_like(volume_background, slowness_back, dtype=np.float64)
    slowness_map_tot[offset[0]:offset[0]+shape[0],offset[1]:offset[1]+shape[1],offset[2]:offset[2]+shape[2],] = slowness_map
    s1=str(slowness_min)
    if '.' in s1:
        s1 = s1.rstrip('0').rstrip('.')
    s2=str(slowness_max)
    if '.' in s2:
        s2 = s2.rstrip('0').rstrip('.')
    slowness_map_tot.tofile('slowness_map_shapeoutter'+str(shape_background[0])+'x'+str(shape_background[1])+'x'+str(shape_background[2])+'_shape'+str(shape[0])+'x'+str(shape[1])+'x'+str(shape[2])+'_off'+str(offset[0])+'x'+str(offset[1])+'x'+str(offset[2])+'_range'+s1+'_'+s2+'_nshape'+str(num_shapes)+'.bin')

    return slowness_map_tot

def main():
    parser = argparse.ArgumentParser(description="Generate a 3D array with randomized, piece-wise constant slowness values with random geometric shapes.")
    parser.add_argument('--shape', nargs=3, type=int, default=[100, 100, 100], help="Shape of the 3D volume as three integers (depth, height, width). Default is 100x100x100.")
    parser.add_argument('--shape_background', nargs=3, type=int, default=[120, 120, 120], help="Shape of the background 3D volume as three integers (depth, height, width). Default is 120x120x120.")
    parser.add_argument('--offset', nargs=3, type=int, default=[5, 5, 5], help="Offsets between the lowerleft corners of the volume and the background.")
    parser.add_argument('--num_shapes', type=int, default=50, help="Number of random shapes. Default is 50.")
    parser.add_argument('--plot', type=int, default=1, help="Whether to visualize.")
    parser.add_argument('--slowness_min', type=float, default=1.0, help="Minimum slowness value. Default is 1.0.")
    parser.add_argument('--slowness_max', type=float, default=3.0, help="Maximum slowness value. Default is 3.0.")
    parser.add_argument('--slowness_back', type=float, default=2.0, help="Background slowness value. Default is 2.0.")

    args = parser.parse_args()

    slowness_map_tot = generate_shapes(tuple(args.shape_background),tuple(args.shape), tuple(args.offset), args.num_shapes, args.slowness_min, args.slowness_max, args.slowness_back)


    # Visualization: Show a slice of the slowness map
    if(args.plot==1):
        plt.figure(figsize=(10, 5))
        plt.imshow(slowness_map_tot[args.shape_background[0] // 2], cmap='Spectral')
        plt.colorbar(label='Slowness (s/m)')
        plt.title('Central Slice of Slowness Map')
        plt.show()

if __name__ == "__main__":
    main()