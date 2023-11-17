import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
from deap import base, creator, tools, algorithms
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
from scipy.special import comb
import mrcfile
from pickpoints import generate_data_points
from cupyx.scipy.signal import correlate2d
import random
import multiprocessing
from multiprocessing import Pool
import json
from scipy.ndimage import zoom


# def bezier_curve(control_points, t):
#     B = np.outer((1 - t) ** 3, control_points[0]) + \
#         np.outer(3 * (1 - t) ** 2 * t, control_points[1]) + \
#         np.outer(3 * (1 - t) * t ** 2, control_points[2]) + \
#         np.outer(t ** 3, control_points[3])
#     return B.squeeze()
def bezier_curve(control_points, t):
    n = len(control_points) - 1
    B = np.zeros_like(control_points[0], dtype=float)
    for i, point in enumerate(control_points):
        B += comb(n, i) * (1 - t) ** (n - i) * t ** i * point
    return B
# def bezier_curve_derivative(control_points, t):
#     control_points = np.array(control_points)
#     B_prime = 3 * (1 - t) ** 2 * (control_points[1] - control_points[0]) + \
#               6 * (1 - t) * t * (control_points[2] - control_points[1]) + \
#               3 * t ** 2 * (control_points[3] - control_points[2])
#     return B_prime
def bezier_curve_derivative(control_points, t):
    n = len(control_points) - 1
    B_prime = np.zeros(2)
    for i in range(n):
        coef = n * comb(n-1, i) * t**i * (1-t)**(n-1-i)
        B_prime += coef * (control_points[i+1] - control_points[i])
    return B_prime

# def bezier_curvature(control_points, t):
#     dB0 = -3 * (1 - t) ** 2
#     dB1 = 3 * (1 - t) ** 2 - 6 * t * (1 - t)
#     dB2 = 6 * t * (1 - t) - 3 * t ** 2
#     dB3 = 3 * t ** 2

#     ddB0 = 6 * (1 - t)
#     ddB1 = 6 - 18 * t
#     ddB2 = 18 * t - 6
#     ddB3 = 6 * t

#     p = control_points
#     dx = sum([p[i, 0] * [dB0, dB1, dB2, dB3][i] for i in range(4)])
#     dy = sum([p[i, 1] * [dB0, dB1, dB2, dB3][i] for i in range(4)])
#     ddx = sum([p[i, 0] * [ddB0, ddB1, ddB2, ddB3][i] for i in range(4)])
#     ddy = sum([p[i, 1] * [ddB0, ddB1, ddB2, ddB3][i] for i in range(4)])

#     curvature = abs(dx * ddy - dy * ddx) / (dx * dx + dy * dy) ** 1.5
#     return curvature
def bezier_curvature(control_points, t, threshold=1e-6, high_curvature_value=1e6):
    n = len(control_points) - 1
    
    dB = np.zeros(2)
    ddB = np.zeros(2)
    
    for i in range(n):
        coef = n * comb(n-1, i) * t**i * (1-t)**(n-1-i)
        dB += coef * (control_points[i+1] - control_points[i])
        
    for i in range(n-1):
        coef = n * (n-1) * comb(n-2, i) * t**i * (1-t)**(n-2-i)
        ddB += coef * (control_points[i+2] - 2*control_points[i+1] + control_points[i])

    dx, dy = dB
    ddx, ddy = ddB
    
    magnitude_squared = dx * dx + dy * dy
    
    # 规避除数接近零的问题
    if magnitude_squared < threshold:
        return high_curvature_value

    curvature = abs(dx * ddy - dy * ddx) / magnitude_squared ** 1.5
    return curvature
def gaussian_pdf(x, mean, std):
    coefficient = 1.0 / (std * np.sqrt(2 * np.pi))
    exponential = np.exp(- (x - mean) ** 2 / (2 * std ** 2))
    return coefficient * exponential

def gaussian2(x, membrane_dist, std1, std2):
    mean1 = - membrane_dist / 2
    mean2 = membrane_dist / 2
    gaussian1 = gaussian_pdf(x, mean1, std1)
    gaussian2 = gaussian_pdf(x, mean2, std2)
    g = gaussian1 + gaussian2
    g = g / np.max(g)
    return g


def points_along_normal(control_points, t_values):
    derivatives = cp.asarray([bezier_curve_derivative(control_points, t) for t in t_values])
    # Normalize the derivatives to get the unit tangent vectors
    tangents = derivatives / cp.linalg.norm(derivatives, axis=-1)[:, np.newaxis]
    # Compute the normals from the tangents
    normals = cp.zeros_like(tangents)
    normals[:, 0] = -tangents[:, 1]
    normals[:, 1] = tangents[:, 0]
    return normals

def coarsefit_evaluate_individual(ind, obj, data_points):
    return (obj.loss_function(np.array(ind), data_points),)


def evaluate_individual(ind, obj, base_image):
    # assert len(ind) == 8, f"Unexpected length of individual: {len(ind)}"
    return (obj.cross_correlation_fitness(np.array(ind).reshape(obj.degree+1, 2), base_image, obj.penalty_threshold),)
def init_cuda():
    x = cp.array([1, 2, 3])
    del x

def cupy_cdist(XA, XB):
    """
    Compute the pairwise distances between two sets of points using CuPy.
    
    Args:
        XA (cupy.ndarray): An array of shape (M, K).
        XB (cupy.ndarray): An array of shape (N, K).
        
    Returns:
        cupy.ndarray: An array of shape (M, N) representing the pairwise distances.
    """
    # Using broadcasting and the fact that (a-b)^2 = a^2 - 2*a*b + b^2
    sq1 = cp.sum(XA**2, axis=1).reshape(-1, 1)
    sq2 = cp.sum(XB**2, axis=1).reshape(1, -1)
    inner_product = cp.dot(XA, XB.T)
    
    return cp.sqrt(sq1 - 2*inner_product + sq2)

class Coarsefit:
    def __init__(self, image, num_points, degree, iteration, num_cpu=20):
        self.image = image
        self.num_points = num_points
        self.img_sz = image.shape[0]
        self.degree = degree
        self.num_cpu = num_cpu
        self.iteration = iteration
    def __call__(self):
        data_points = generate_data_points(self.image, self.num_points)
        control_points = self.coarse_fitting_ga(data_points)
        return control_points
    def loss_function(self, control_points_flat, data_points):
        control_points = control_points_flat.reshape(self.degree+1, 2)
        t_values = np.linspace(0, 1, len(data_points))
        curve_points = np.array([bezier_curve(control_points, t) for t in t_values])
        return np.sum((curve_points - data_points) ** 2)
    def coarse_fitting_ga(self, data_points):
        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
        creator.create("Individual", np.ndarray, fitness=creator.FitnessMin)
        toolbox = base.Toolbox()
        toolbox.register("attr_float", np.random.uniform, -self.img_sz // 2, self.img_sz // 2)
        toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_float, n=(self.degree+1)*2)  # 四个控制点，每个控制点有两个坐标，所以n=8
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)
        # toolbox.register("evaluate", lambda ind: (self.loss_function(np.array(ind), data_points),))
        toolbox.register("evaluate", coarsefit_evaluate_individual, obj=self, data_points=data_points)
        toolbox.register("mate", tools.cxBlend, alpha=0.5)
        toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=1, indpb=0.2)
        toolbox.register("select", tools.selTournament, tournsize=5)
        pop = toolbox.population(n=100)
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("avg", np.mean)
        stats.register("std", np.std)
        stats.register("min", np.min)
        stats.register("max", np.max)
        pool = Pool(processes=self.num_cpu)
        toolbox.register("map", pool.map)
        algorithms.eaSimple(pop, toolbox, cxpb=0.6, mutpb=0.3, ngen=self.iteration, 
                                stats=stats, halloffame=None, verbose=True)
        pool.close()
        pool.join()
        best_ind = tools.selBest(pop, 1)[0]
        return np.array(best_ind).reshape(self.degree+1, 2)

class GA_Refine:
    def __init__(self, image, pixel_size, penalty_threshold, dithering_range, iterations, num_cpu=20):
        self.image = image
        self.img_sz = image.shape[0]
        physical_membrane_dist = 35.0
        self.pixel_size = pixel_size
        self.penalty_threshold = penalty_threshold
        self.mem_dist = int(physical_membrane_dist / self.pixel_size)
        self.dithering_range = dithering_range
        self.iterations = iterations
        self.num_cpu = num_cpu

    def __call__(self, initial_control_points, image):
        self.degree = len(initial_control_points) - 1
        self.initial_control_points = initial_control_points
        refined_control_points = self.ga_refine_controlpoints(image, initial_control_points)
        return refined_control_points

    # def average_1d(self, image, fitted_points, normals, extra_mem_dist):
    #     average_1d_lst = []
    #     for membrane_dist in range(-extra_mem_dist, extra_mem_dist+1):
    #         normals_points = fitted_points + membrane_dist * normals
    #         # Ensure the points are within the image boundaries
    #         mask = (normals_points[:, 0] > 0) & (normals_points[:, 0] < image.shape[1]) & \
    #             (normals_points[:, 1] > 0) & (normals_points[:, 1] < image.shape[0])
    #         normals_points = normals_points[mask]
    #         # Get interpolated gray values for the normal points
    #         interpolated_values = bilinear_interpolation(image, normals_points[:, 0], normals_points[:, 1])
    #         # convert nan to 0
    #         interpolated_values = np.nan_to_num(interpolated_values)
    #         average_1d_lst.append(np.mean(interpolated_values))
    #     return average_1d_lst
    def generate_2d_average(self, image, fitted_points, average_1d_lst, membrane_distance):
        fitted_points = cp.array(fitted_points)
        # Create an empty image of the same size as the original
        new_image = cp.zeros_like(image)
        # Create a meshgrid for the image
        y, x = cp.mgrid[:image.shape[0], :image.shape[1]]
        coords = cp.stack((x, y), axis=-1).reshape(-1, 2)
        # Calculate distances from each pixel to the fitted_curve_points
        distances = cupy_cdist(coords, fitted_points)
        min_distances = cp.min(distances, axis=1).reshape(image.shape)
        edge_sigma = 5
        mask_within_distance = cp.abs(min_distances) <= membrane_distance
        gray_value = cp.exp(-(cp.abs(min_distances) - membrane_distance)**2 / (2 * edge_sigma**2))
        mask_small_gray_value = gray_value < 0.001
        mask_outside_distance = ~mask_within_distance
        membrane_mask = cp.zeros_like(min_distances)
        membrane_mask[mask_within_distance] = 1
        membrane_mask[mask_outside_distance & ~mask_small_gray_value] = gray_value[mask_outside_distance & ~mask_small_gray_value]

        # Use the distances to interpolate values from average_1d_lst
        f = interp1d(np.arange(-membrane_distance, membrane_distance+1), average_1d_lst, kind='linear', bounds_error=False, fill_value=0)
        new_image = cp.asarray(f(min_distances.get())) * membrane_mask

        return new_image, membrane_mask
    
    def cross_correlation_fitness(self, control_points, base_image, penalty_threshold):
        base_image = cp.asarray(base_image)
        fitted_curve_points, filtered_t_values = generate_curve_within_boundaries(control_points, base_image.shape, 0.01)
        curvatures = [bezier_curvature(control_points, t) for t in filtered_t_values]
        curvatures_abs = [abs(curvature) for curvature in curvatures]
        control_points_out_of_bounds = [max(0, point[0] - self.img_sz, point[1] - self.img_sz, -point[0], -point[1]) for point in control_points]
        control_points_penalty = sum([10*(2**out_of_bound) for out_of_bound in control_points_out_of_bounds])
        if any(curvature > penalty_threshold for curvature in curvatures_abs):
            cur_penalty = 1e4*max(curvatures_abs)
        else:
            cur_penalty = 0
        std = 3
        new_image, membrane_mask = self.generate_2d_average(base_image, fitted_curve_points, gaussian2(np.arange(-self.mem_dist*2, self.mem_dist*2+1), self.mem_dist, std, std), self.mem_dist*2)
        cc_value = correlate2d(new_image, base_image*membrane_mask, mode='valid')
        return cc_value.get() - cur_penalty - control_points_penalty
    def attr_around_initial(self, index):
        return self.initial_control_points.flatten()[index] + np.random.uniform(-self.dithering_range, self.dithering_range)
    def custom_mutGaussian(self, individual, mu, sigma, indpb):
        for i in range(len(individual)):
            if random.random() < indpb:
                individual[i] += random.gauss(mu, sigma)
        return individual,
    def ga_refine_controlpoints(self, base_image, initial_control_points):
        initial_control_points = np.array(initial_control_points)
        creator.create("FitnessMax", base.Fitness, weights=(1.0,))
        creator.create("Individual", list, fitness=creator.FitnessMax)

        toolbox = base.Toolbox()
        # toolbox.register("attr_float", self.attr_around_initial, index=np.arange(8))
        toolbox.register("attr_float", self.attr_around_initial, index=np.arange((self.degree+1)*2))
        toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.attr_float)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)
        toolbox.register("evaluate", evaluate_individual, obj = self, base_image = base_image)
        toolbox.register("mate", tools.cxUniform, indpb=0.5)
        toolbox.register("mutate", self.custom_mutGaussian, mu=0, sigma=1, indpb=0.2)
        toolbox.register("select", tools.selTournament, tournsize=5)
        pop = toolbox.population(n=30)
        pop.append(creator.Individual(initial_control_points.flatten()))

        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("avg", np.mean)
        stats.register("std", np.std)
        stats.register("min", np.min)
        stats.register("max", np.max)
        pool = Pool(processes=self.num_cpu, initializer=init_cuda)
        toolbox.register("map", pool.map)
        algorithms.eaSimple(pop, toolbox, cxpb=0.6, mutpb=0.4, ngen=self.iterations, 
                                stats=stats, halloffame=None, verbose=True)
        pool.close()
        pool.join()
        best_individual = tools.selBest(pop, 1)[0]
        best_control_points = np.array(best_individual).reshape(self.degree+1, 2)
        # print("Best Control Points:")
        # print(best_control_points)
        return best_control_points

def generate_curve_within_boundaries(control_points, image_shape, step):
    t_values = []
    t = 0

    # Find the first t value within the image boundaries
    while t < 2:  # Limit to avoid infinite loop
        point = bezier_curve(control_points, t)
        if 0 <= point[0] < image_shape[0] and 0 <= point[1] < image_shape[1]:
            t_values.append(t)
            break
        t += step

    # Generate curve points moving forward
    while True:
        t += step
        point = bezier_curve(control_points, t)
        if 0 <= point[0] < image_shape[0] and 0 <= point[1] < image_shape[1]:
            t_values.append(t)
        else:
            break

    # Reset t to the initial point and generate curve points moving backward
    t = t_values[0] - step
    while t > -2:  # Limit to avoid infinite loop
        point = bezier_curve(control_points, t)
        if 0 <= point[0] < image_shape[0] and 0 <= point[1] < image_shape[1]:
            t_values.insert(0, t)  # Insert at the beginning
            t -= step
        else:
            break
    fitted_curve_points = np.array([bezier_curve(control_points, t_val) for t_val in t_values])
    return np.array(fitted_curve_points), np.array(t_values)

if __name__ == '__main__':
    multiprocessing.set_start_method('spawn', force=True)
    with mrcfile.open('/data3/kzhang/cryosparc/CS-vsv/J354/templates_selected.mrc') as f:
        image = f.data[2]
    image = zoom(image, 2)
    coarsefit = Coarsefit(image, 600, 3, 300, 20)
    initial_control_points = coarsefit()
    ga_refine = GA_Refine(image, 1.068, 0.05, 50, 700, 18)
    refined_control_points = ga_refine(initial_control_points, image)
    refined_control_points = np.array(refined_control_points)
    fitted_curve_points, t_values = generate_curve_within_boundaries(refined_control_points, image.shape, 0.01)
    # save the control points in JSON format
    with open('control_points.json', 'w') as f:
        json.dump(refined_control_points.tolist(), f)
    plt.imshow(image, cmap='gray')
    plt.plot(fitted_curve_points[:, 0], fitted_curve_points[:, 1], 'r-')
    plt.plot(refined_control_points[:, 0], refined_control_points[:, 1], 'g.')
    plt.show()
