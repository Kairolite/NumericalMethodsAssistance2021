import numpy as np


def dataset():
    x = [0, 1, 2, 3, 4, 5]
    y = [2.1, 7.7, 13.6, 27.2, 40.9, 61.1]

    # degree of polynomial regression,
    # m-order polynomial where n >= m+1
    m = 3

    return x, y, m


def equation_matrix():
    # Generate the matrices: A = mat_x; b = mat_y;
    x, y, m = dataset()

    mat_y = [[sum([(xi ** i) * yi for xi, yi in zip(x, y)])] for i in range(m + 1)]

    sum_power_x = [sum([xi ** i for xi in x]) for i in range(2 * m + 1)]

    mat_x = []
    for row in range(m + 1):
        rowlist = []
        for col in range(m + 1):
            rowlist.append(sum_power_x[col + row])
        mat_x.append(rowlist)

    return mat_x, mat_y


def solve_coeff():
    # Solve the coefficient: a = x; Solve Ax = b
    x, y = equation_matrix()

    x = np.array(x)
    y = np.array(y)

    a = np.linalg.solve(x, y)

    return a


def regression():
    # Generate and print the regression equation
    a = solve_coeff()

    lq_fit = f"y = {a[0]} + {a[1]}*x"

    if len(a > 2):
        for i in range(len(a) - 2):
            term = f" + {a[i + 2]}*x^{i + 2}"
            lq_fit = lq_fit + term

    print(f"Least-Squares Fit\n"
          f"{lq_fit}\n")

    return lq_fit


if __name__ == '__main__':
    regression()
