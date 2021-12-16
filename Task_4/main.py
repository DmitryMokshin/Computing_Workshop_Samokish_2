import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.special import legendre

import_directory = './Data/'


def result_graphics_table(n):
    file_name = ''
    file_name_coef = ''
    file_name_error = ''
    file_name_error_method = ''
    if len(str(n)) == 2:
        file_name = 'result' + str(n) + '.csv'
        file_name_coef = 'result_coef' + str(n) + '.csv'
        file_name_error = 'result_error' + str(n) + '.csv'
        file_name_error_method = 'result_error_method' + str(n) + '.csv'
    if len(str(n)) == 1:
        file_name = 'result0' + str(n) + '.csv'
        file_name_coef = 'result_coef0' + str(n) + '.csv'
        file_name_error = 'result_error0' + str(n) + '.csv'
        file_name_error_method = 'result_error_method0' + str(n) + '.csv'
    dataframe = pd.read_csv(import_directory + file_name, sep='|')
    dataframe_coef = pd.read_csv(import_directory + file_name_coef, sep='|')
    dataframe_error = pd.read_csv(import_directory + file_name_error, sep='|')
    dataframe_error_method = pd.read_csv(import_directory + file_name_error_method, sep='|')

    x = np.arange(0, 1.0, 0.01)

    for column in dataframe.columns:
        if column != 'x':
            plt.plot(dataframe['x'], dataframe[column], label=column.replace('y,', ''), linestyle='-.')
    # plt.plot(x, np.cos(x), label='Solution', linestyle='-')
    plt.legend()
    plt.xlabel('X', fontsize=17)
    plt.ylabel('Y', fontsize=17)
    plt.xlim(0, 1.0)
    plt.grid()
    plt.savefig('./Pictures/' + 'Graphicsn=' + str(n) + '.png', format='png', dpi=800)
    plt.close()

    for column in dataframe_error.columns:
        if column != 'x':
            plt.plot(dataframe_error['x'], dataframe_error[column], label=column.replace(r'\hat{K}u-f,', ''),
                     linestyle='-.')

    plt.legend()
    plt.xlabel('X', fontsize=17)
    plt.ylabel(r'$\hat{K}u-f$', fontsize=17)
    plt.xlim(0, 1.0)
    plt.grid()
    plt.savefig('./Pictures/' + 'Graphicserrorn=' + str(n) + '.png', format='png', dpi=800)
    plt.close()

    for column in dataframe_error_method.columns:
        if column != 'x':
            plt.plot(dataframe_error_method['x'], dataframe_error_method[column],
                     label=column.replace(r'\hat{A}u-f,', ''), linestyle='-.')

    plt.legend()
    plt.xlabel('X', fontsize=17)
    plt.ylabel(r'$\hat{A}u-f_1$', fontsize=17)
    plt.xlim(0, 1.0)
    plt.grid()
    plt.savefig('./Pictures/' + 'Graphicserrornmethod=' + str(n) + '.png', format='png', dpi=800)
    plt.close()

    dataframe_error['x'] = dataframe_error['x'].astype(str)
    dataframe_error_method['x'] = dataframe_error_method['x'].astype(str)

    file = open("./Result_tab/" + "result_texn" + str(n) + ".txt", 'w')
    file.write(dataframe.to_latex(index=False, longtable=True, float_format='%.10f'))
    file.close()

    file = open("./Result_tab/" + "resultcoef_texn" + str(n) + ".txt", 'w')
    file.write(dataframe_coef.to_latex(index=False, longtable=True, float_format='%.10f'))
    file.close()

    file = open("./Result_tab/" + "resulterror_texn" + str(n) + ".txt", 'w')
    file.write(dataframe_error.to_latex(index=False, longtable=True, float_format='%.6e'))
    file.close()

    file = open("./Result_tab/" + "resulterrormethod_texn" + str(n) + ".txt", 'w')
    file.write(dataframe_error_method.to_latex(index=False, longtable=True, float_format='%.6e'))
    file.close()

    return


def result_graphics_legendre(n):
    file_name = ''
    if len(str(n)) == 2:
        file_name = 'resultpolynomlegendre' + str(n) + '.csv'
    if len(str(n)) == 1:
        file_name = 'resultpolynomlegendre0' + str(n) + '.csv'
    dataframe = pd.read_csv(import_directory + file_name, sep='|')

    file = open("./Result_tab/" + "Legendreresult_texn" + str(n) + ".txt", 'w')
    file.write(dataframe.to_latex(index=False, longtable=True, float_format='%.10f'))
    file.close()

    for column in dataframe.columns:
        if column != 'x':
            N = int(column.replace('P, N = ', ''))
            # plt.plot(dataframe['x'], dataframe[column] - legendre(N)(2.0 * dataframe['x'] - 1.0),
            #          label=column.replace('P, ', ''), linestyle='-.')
            plt.plot(dataframe['x'], legendre(N)(2.0 * dataframe['x'] - 1.0), label=N, linestyle='-')

    plt.legend()
    plt.xlabel('X', fontsize=17)
    plt.ylabel('Y', fontsize=17)
    plt.xlim(0, 1.0)
    plt.grid()
    plt.savefig('./Pictures/' + 'LegendreGraphics_n=' + str(n) + '.png', format='png', dpi=800)
    plt.close()
    return


n = [3, 5, 7, 9, 10, 12, 15, 18, 20]
# n = [15 ]

for i in range(len(n)):
    result_graphics_table(n[i])
result_graphics_legendre(20)
result_graphics_legendre(5)
