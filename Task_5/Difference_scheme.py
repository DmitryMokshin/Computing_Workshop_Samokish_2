import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

d = 0.25
T = 0.5


def phi(x):
    return 1.0 - np.power(x, 2)


def a(x):
    return x * (1.0 - x)


def alpha(t):
    return (1.0 + t) * (1.0 - t / d) if 0 <= t < d else 0.0


def explicit_difference_scheme(N):
    h = 1.0 / N
    M = int(N ** 2 * (2.0 * T))
    tau = T / M
    nu = tau / np.power(h, 2)

    u_present = np.array([phi(i * h) for i in range(N + 1)])
    u_future = np.zeros(N + 1)
    file_out = open('./Data/result_explicit_N=' + str(N) + '_M=' + str(M) + '.csv', 'w')
    file_out.write(f' ')
    for i in range(N + 1):
        file_out.write(f'| {i * h:.4f}')
    file_out.write('\n')
    file_out.write(f'{0.0}')
    for u_present_i in u_present:
        file_out.write(f'| {u_present_i}')
    file_out.write(f'\n')
    for k in range(1, M + 1):
        u_future[1: N] = u_present[1: N] + np.array(
            [nu * (u_present[i + 1] - 2.0 * u_present[i] + u_present[i - 1]) - a(i * h) * u_present[i] * tau for i in
             range(1, N)])
        u_future[0] = (4.0 * u_future[1] - u_future[2] - 2.0 * h * alpha(k * tau)) / 3.0
        u_future[N] = (4.0 * u_future[N - 1] - u_future[N - 2]) / 3.0
        u_present = u_future
        if k % (M // 10) == 0.0:
            file_out.write(f'{k * tau}')
            for u_present_i in u_present:
                file_out.write(f'| {u_present_i}')
            file_out.write(f'\n')
    return


def tridiagonal_matrix_algorithm(matrix_A, vector_f):
    N = len(vector_f)
    vector_alpha = np.zeros(N)
    vector_beta = np.zeros(N)
    vector_result = np.zeros(N)
    vector_alpha[0] = matrix_A[2, 0] / matrix_A[1, 0]
    vector_beta[0] = -vector_f[0] / matrix_A[1, 0]
    for i in range(1, N):
        m = matrix_A[1, i] - matrix_A[0, i] * vector_alpha[i - 1]
        vector_alpha[i] = matrix_A[2, i] / m
        vector_beta[i] = (matrix_A[0, i] * vector_beta[i - 1] - vector_f[i]) / m
    vector_result[-1] = vector_beta[-1]
    for i in range(1, N):
        vector_result[N - i - 1] = vector_alpha[N - 1 - i] * vector_result[N - i] + vector_beta[N - 1 - i]
    return vector_result


def implicit_difference_scheme(N, sigma):
    h = 1.0 / N
    M = int(N ** 2 * (2.0 * T))
    tau = T / M

    u_present = np.array([phi(i * h) for i in range(N + 1)])

    file_out = open('./Data/result_implicit_N=' + str(N) + '_M=' + str(M) + '.csv', 'w')
    file_out.write(f' ')
    for i in range(N + 1):
        file_out.write(f'| {i * h:.4f}')
    file_out.write('\n')

    file_out.write(f'{0.0}')
    for u_present_i in u_present:
        file_out.write(f'| {u_present_i}')
    file_out.write(f'\n')

    for k in range(1, M + 1):
        D = np.zeros(N + 1)
        A_matrix = np.zeros((3, N + 1))
        A_matrix[0, :] = -sigma / np.power(h, 2)
        A_matrix[1, :] = -np.array([1.0 / tau + 2.0 * sigma / np.power(h, 2) + a(i * h) * sigma for i in range(N + 1)])
        A_matrix[2, :] = -sigma / np.power(h, 2)
        A_matrix[:, 0] = np.array([0.0, 1.0 / h, 1.0 / h])
        A_matrix[:, -1] = np.array([-1.0 / h, -1.0 / h, 0.0])
        D[0] = alpha(k * tau)
        D[1:N] = np.array([u_present[i] / tau + (1.0 - sigma) * (
                (u_present[i + 1] - 2.0 * u_present[i] + u_present[i - 1]) / np.power(h, 2) - a(i * h) * u_present[
            i]) for i in range(1, N)])
        D[N] = 0.0

        u_present = tridiagonal_matrix_algorithm(A_matrix, D)
        if k % (M // 10) == 0.0:
            file_out.write(f'{k * tau}')
            for u_present_i in u_present:
                file_out.write(f'| {u_present_i}')
            file_out.write(f'\n')
    return


def result_graphics(N, method):
    global file
    h = 1.0 / N
    x = [i * h for i in range(N + 1)]
    M = int(N ** 2 * (2.0 * T))

    file_name = 'result_' + method + '_N=' + str(N) + '_M=' + str(M) + '.csv'
    import_directory = './Data/'
    dataframe = pd.read_csv(import_directory + file_name, sep='|', index_col=' ')
    for i in range(11):
        plt.plot(x, dataframe.iloc[i], label='t = ' + str(dataframe.index[i]), linestyle='--')
    plt.legend()
    plt.xlabel('X', fontsize=17)
    plt.ylabel('Y', fontsize=17)
    plt.title('N = ' + str(N) + ', M = ' + str(M))
    plt.grid()
    plt.savefig('./Pictures/' + method + 'graphics_n=' + str(N) + '.png', format='png', dpi=800)
    plt.close()
    columns = dataframe.columns
    file = open("./Result_table/" + 'result_' + method + '_tex_n_' + str(N) + '.txt', 'w')
    for i in range(N // 5 - 1):
        # file = open("./Result_table/" + 'result_' + method + '_tex_n_' + str(i + 1) + '=' + str(N) + '.txt', 'w')
        file.write(dataframe[columns[5 * i:5 * (i + 1)]].to_latex(index=True, longtable=True, float_format='%.9f'))
        file.write('\n')
        # file.close()
    # file = open("./Result_table/" + 'result_' + method + '_tex_n_' + str(N // 5) + '=' + str(N) + '.txt', 'w')
    file.write(dataframe[columns[N - 5:N + 1]].to_latex(index=True, longtable=True, float_format='%.9f'))
    file.close()
    return


n = [10, 20, 30, 40]

[explicit_difference_scheme(i) for i in n]
[result_graphics(i, 'explicit') for i in n]

[implicit_difference_scheme(i, 0.5) for i in n]
[result_graphics(i, 'implicit') for i in n]
