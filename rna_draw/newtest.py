import numpy as np

def normalize_angle(start_angle, end_angle):
    if start_angle < 0:
        start_angle += 2 * np.pi
    if end_angle < 0:
        end_angle += 2 * np.pi
    if start_angle > end_angle:
        end_angle += 2 * np.pi
        
    start_angle = start_angle % (2 * np.pi)
    end_angle = end_angle % (2 * np.pi)

    if start_angle >= end_angle:
        temp = start_angle
        start_angle = end_angle
        end_angle = temp

    return start_angle, end_angle

if __name__ == "__main__":
    # Run tests with randomized angles
    num_tests = 10
    for _ in range(num_tests):
        start_angle = np.random.uniform(-np.pi, np.pi)
        end_angle = np.random.uniform(-np.pi, np.pi)
        
        print(f"Initial start angle: {np.rad2deg(start_angle)}, end angle: {np.rad2deg(end_angle)}")
        norm_start, norm_end = normalize_angle(start_angle, end_angle)
        print(f"Normalized start angle: {np.rad2deg(norm_start)}, end angle: {np.rad2deg(norm_end)}\n")
