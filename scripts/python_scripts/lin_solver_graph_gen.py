__author__ = 'Kazem'
import csv
import os
import glob
from os.path import isfile, join
import sys
import getopt
import collections
import matplotlib.pyplot as plt
import numpy as np
import fnmatch
import pandas as pd
from enum import Enum
from matplotlib.font_manager import FontProperties
from matplotlib import colors as mcolors
from os.path import basename
from operator import itemgetter, attrgetter, methodcaller
from osqp_utils import compute_failure_rates,compute_performance_profiles,plot_performance_profiles,\
    exclude_items,compute_speedup,exclude_items_general,get_name, compute_speedup_general

INV='N/A'

def dir_contents(path, type):
    matches = []
    for root, dirnames, filenames in os.walk(path):
        for filename in fnmatch.filter(filenames, '*.'+ type):
            matches.append(os.path.join(root, filename))
    return matches

class solver_time:
    def __init__(self):
        self.fact_time = np.array([])
        self.analysis_time = np.array([])
        self.ordering_time = np.array([])
        self.solve_time = np.array([])
        self.refine_time = np.array([])

class MatrixRecord:
    def __init__(self, matrix_name='', nnz=0, dim1=0, dim2=0,
                 nnz_l=0, ordering_method='',
                 tree_height=0, norm2=0.0, residual=0.0, bwd_error=0.0, fwd_error=0.0):
        self.matrix_name = matrix_name
        self.nnz = nnz
        self.dim1 = dim1
        self.dim2 = dim2
        self.nnz_l = nnz_l
        self.ordering_method = ordering_method
        self.parallel_time = solver_time()
        self.serial_time = solver_time()
        self.tree_height = tree_height
        self.norm2 = norm2
        self.residual = residual
        self.bwd_error = bwd_error
        self.fwd_error = fwd_error



class RunRecord:
    def __init__(self, tool_name='', target_name='', num_thread=0):
        self.tool_name = tool_name
        self.num_thread = num_thread
        self.target_name = target_name
        self.mat_rec_list = []

OPTIMAL = 'optimal'
MAX_VAL = 1600
class Visualizer:

    def __init__(self, path):
        self.path = path
        self.run_rec_list = []
        self.headerList = []
        self.empty_dic = []
        self.font0 = FontProperties()
        self.qp2size = {}
        self.in_log = {}
        self.sympiler_header = ['Matrix Name', 'Row', 'Col', 'Nonzero', 'Factor Nonzero',
                                'Ordering Time (sec)', 'Analysis Time (sec)',
                                'Factorization Time (sec)', 'Solve Time (sec)',
                                'Iterative Refinement Time (sec)',
                                'Error Norm2-Before', 'Rel Error Norm2-Before',
                                'Infinity Norm-Before', 'BWD Error-Before', 'Residual-Before',
                                'Number of Threads',
                                'Target Name', 'Total Time (sec)', 'Tool Name',
                                'Iterative Refinement Method', 'Chunk Size',
                                'Error Norm2-After', 'Rel Error Norm2-After',
                                'Infinity Norm-After', 'BWD Error-After', 'Residual-After',
                                'Number of W-partition', 'Initial Cut', 'BLAS Threads',
                                'Final Cuts', 'Height', 'Number of Supernodes', 'Number of Pivots',
                                'Precision-Before', 'Precision-After', 'Number of Refinement Iterations',
                                'Factorization (GFLOPS)', 'QP Time', 'Error Inf'
                                ]

    class LogHeader(Enum):
        MATN = 0
        ROW = 1
        COL = 2
        NNZ = 3
        NNZL = 4
        ORDT = 5
        ANAT = 6
        FACT = 7
        SOLT = 8
        ITERT = 9
        NRM1 = 10
        NRM2 = 11
        NRMI = 12
        BWDE = 13
        RESD = 14
        THRN = 15
        TRGT = 16
        TOTT = 17
        TOOL = 18
        IREF = 19
        CHNK = 20
        NR1A = 21
        NR2A = 22
        NRIA = 23
        BWEA = 24
        REDA = 25
        WPAR = 26
        ICUT = 27
        BLAS = 28
        FCUT = 29
        HGHT = 30
        SUPN = 31
        NPIV = 32
        PREC = 33
        PREA = 34
        NREF = 35
        FLPS = 36
        QPTT = 37
        ERRI = 38

    def read_csv(self, report_path):
        #report = RunRecord()
        report = []
        with open(report_path) as csvFile:
            reader = csv.DictReader(csvFile)
            # cur_row = MatrixRecord()
            self.headerList = reader.fieldnames
            for row in reader:
                report.append(row)
        if len(report) > 0:
            self.empty_dic = report[0].copy()
            for key, value in self.empty_dic.items():
                self.empty_dic[key] = ''
        return report

    def read_sympiler_log(self, report_path, idx_eq, idx_sort, f_sort):
        sym_log = self.read_csv(report_path)
        red_log = self.list_reduction(sym_log, idx_eq, idx_sort, f_sort)
        # extracting matrix name from the path
        for i in range(len(red_log)):
            full_name = red_log[i][self.sympiler_header[self.LogHeader.MATN.value]]
            mat_full_name = full_name.split('/')
            mat_name = mat_full_name[len(mat_full_name) - 1].split('.')[0]
            red_log[i][self.sympiler_header[self.LogHeader.MATN.value]] = mat_name
        red_sort_log = sorted(red_log, key=itemgetter(self.sympiler_header[idx_eq]))
        self.in_log = red_sort_log
        return self.in_log

    def convert_log(self, log_in):
        list_cnvrt = []
        prob_name = self.sympiler_header[self.LogHeader.MATN.value]
        fact_time = self.sympiler_header[self.LogHeader.FACT.value]
        row_str = self.sympiler_header[self.LogHeader.ROW.value]
        nnz_str = self.sympiler_header[self.LogHeader.NNZ.value]

        for it in log_in:
            it[prob_name] = get_name(it[prob_name])
            it[fact_time] = float(it[fact_time])
            it[row_str] = int(it[row_str])
            it[nnz_str] = int(it[nnz_str])
        list_cnvrt = log_in
        return list_cnvrt


    def log_export(self, log_list, problems_type):
        num_converged = 0
        avg_solve = 0
        perf_profile = []
        tool_name = log_list[0][self.sympiler_header[self.LogHeader.TOOLN.value]]
        # file_name = os.path.join(".", "results", problems_type,tool_name,
        #                             "statistics.txt")
        file_name = "results_" + problems_type +"_" + tool_name + "_statistics.txt"
        f = open(file_name, "w")
        f.write(self.sympiler_header[self.LogHeader.PROBN.value])
        f.write(",")
        f.write(self.sympiler_header[self.LogHeader.TOOLN.value])
        f.write(",run_time")
        #f.write(self.qp_header[self.qp_log_header.STIME.value])
        f.write(",status")
        #f.write(self.qp_header[self.qp_log_header.STATS.value])
        f.write("\n")
        for el in log_list:
            f.write(el[self.sympiler_header[self.LogHeader.PROBN.value]])
            f.write(",")
            f.write(el[self.sympiler_header[self.LogHeader.TOOLN.value]])
            f.write(",")
            f.write(el[self.sympiler_header[self.LogHeader.STIME.value]])
            f.write(",")
            f.write(el[self.sympiler_header[self.LogHeader.STATS.value]])
            f.write("\n")
        f.close()
        return num_converged, avg_solve, perf_profile

    def compute_scattered(self, *list_args, idx_x, idx_y, idx_y2, y_label='',
                       file_name='', y_scale="linear"):
        if len(list_args) <= 0:
            print("Wrong lists to visualize.")
            return -1
        n_groups = len(list_args[0])
        tool_name = []

        for lst in list_args:
            if len(lst) != n_groups:
                print("Wrong lists to visualize.")
                return -1
            tool_name.append(lst[0][self.sympiler_header[self.LogHeader.TOOLN.value]])
        num_tool = len(list_args)
        all_data = np.zeros((num_tool,n_groups))
        all_data_rel = np.zeros((num_tool, n_groups))
        mat_list = []
        lst_no = 0
        for lst in list_args:
            for d in range(len(lst)):
                if (lst[d][self.sympiler_header[idx_y2.value]]) == OPTIMAL:
                    all_data[lst_no, d] = lst[d][self.sympiler_header[idx_y.value]]
                else:
                    all_data[lst_no, d] = MAX_VAL
                if lst_no == 0:
                    mat_list.append(lst[d][self.sympiler_header[idx_x.value]])
            lst_no += 1

        #all_data_rel = f_val(all_data)
        # for dim in range(num_tool):
        #     all_data_rel[dim, :] = all_data[dim, :]/all_data[0, :]
        # for v in range(n_groups):
        #     print(mat_list[v],',',all_data[0,v],',',all_data[1,v])

        ## Preparing the graph info
        # Other colors: https://matplotlib.org/gallery/color/named_colors.html
        colors = ('red', 'skyblue', 'indigo', 'olive', 'slateblue', 'magenta',
                  'slategray', 'limegreen', 'maroon', 'teal', 'khaki', 'purple',
                  'r', 'c', 'm', 'y', 'k')
        markers = ('s', 'o', '+', "."	, ",", "v", "^", "<", ">",
                   "1", "2", "3", "4", "8", "s", "p", "P")

        x_label = np.arange(1, len(mat_list) + 1)
        font_axis = self.font0.copy()
        fig = plt.figure(figsize=(20, 5))
        ax1 = fig.add_subplot(111)
        ax1.set_yscale('log')
        for i in range(num_tool):
            ax1.scatter(x_label, all_data[i, :], s=10, c=colors[i],
                        marker=markers[i], label=tool_name[i])
            #ax1.scatter(mat_list[:], all_data[1, :], s=10, c='r', marker="o", label=tool_name[1])
        font_axis.set_weight('black')  # ['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']
        ax1.set_ylabel(y_label, fontsize=12, fontproperties=font_axis)
        #ax1.set_xticks(x_label,minor=True)
        #ax1.set_xticklabels(x_label, fontsize=12)
        plt.xticks(x_label,rotation=55)
        plt.legend(loc='upper right');
        plt.show()

    def compute_bar(self, *list_args, idx_x, idx_y, idx_y2, y_label='',
                       file_name='', y_scale="linear"):
        if len(list_args) <= 0:
            print("Wrong lists to visualize.")
            return -1
        n_groups = len(list_args[0])
        tool_name = []
        for lst in list_args:
            if len(lst) != n_groups:
                print("Wrong lists to visualize.")
                return -1
            tool_name.append(lst[0][self.sympiler_header[self.LogHeader.TOOLN.value]])
        num_tool = len(list_args)
        all_data = np.zeros((num_tool, n_groups))
        all_data_rel = np.zeros((num_tool, n_groups))
        mat_list = []
        lst_no = 0
        for lst in list_args:
            for d in range(len(lst)):
                if (lst[d][self.sympiler_header[idx_y2.value]]) == OPTIMAL:
                    all_data[lst_no, d] = lst[d][self.sympiler_header[idx_y.value]]
                else:
                    all_data[lst_no, d] = MAX_VAL
                if lst_no == 0:
                    mat_list.append(lst[d][self.sympiler_header[idx_x.value]])
            lst_no += 1

        #all_data_rel = f_val(all_data)
        # for dim in range(num_tool):
        #     all_data_rel[dim, :] = all_data[dim, :]/all_data[0, :]
        # for v in range(n_groups):
        #     print(mat_list[v],',',all_data[0,v],',',all_data[1,v])

        ## Preparing the graph info
        # Other colors: https://matplotlib.org/gallery/color/named_colors.html
        colors = ('red', 'skyblue', 'indigo', 'olive', 'slateblue', 'magenta',
                  'slategray', 'limegreen', 'maroon', 'teal', 'khaki', 'purple',
                  'r', 'c', 'm', 'y', 'k')
        markers = ('s', 'o', '+', "."	,"," ,"v"	,"^"	,"<"	,">"	,
                   "1","2"	,"3"	,"4"	,"8"	,"s"	,"p"	,"P"	)
        width = 0.7
        ind = (num_tool + 1) * width * np.arange(len(mat_list))
        x_label = np.arange(1,len(mat_list)+1)
        font_axis = self.font0.copy()
        metric_label = 'solve time(sec)'
        fig = plt.figure(figsize=(30,5))
        ax1 = fig.add_subplot(111)
        ax1.set_yscale('log')
        for i in range(num_tool):
            #ax1.scatter(mat_list[:], all_data[i, :], s=10, c=colors[i],
            #            marker=markers[i], label=tool_name[i])
            ax1.bar(ind + i*width, all_data[i, :], width,
                   label=tool_name[i] + ' ' + metric_label, color=colors[i],
                   align='edge')
            #ax1.scatter(mat_list[:], all_data[1, :], s=10, c='r', marker="o", label=tool_name[1])

        #ax1.set_ylabel(mat_list[:], fontsize=12, fontproperties=font_axis)
        font_axis.set_weight('black')  # ['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']
        ax1.set_ylabel(y_label, fontsize=12, fontproperties=font_axis)
        ax1.set_xticks(ind)
        ax1.set_xticklabels(x_label, fontsize=12)
        plt.xticks(rotation=55)
        plt.legend(loc='upper right')
        plt.show()

    def plot_bar_chart(self, list_logs, idx_x, idx_y, idx_y2, y_label='',
                       file_name='', y_scale="linear"):
        if len(list_logs) <= 0:
            print("Wrong lists to visualize.")
            return -1
        n_groups = len(list_logs[0])
        tool_name = []
        for lst in list_logs:
            if len(lst) != n_groups:
                print("Wrong lists to visualize.")
                return -1
            tool_name.append(lst[0][self.qp_header[self.qp_log_header.TOOLN.value]])
        num_tool = len(list_logs)
        all_data = np.zeros((num_tool,n_groups))
        all_data_rel = np.zeros((num_tool, n_groups))
        mat_list = []
        lst_no = 0
        for lst in list_logs:
            for d in range(len(lst)):
                if (lst[d][self.qp_header[idx_y2.value]]) == OPTIMAL:
                    all_data[lst_no, d] = lst[d][self.qp_header[idx_y.value]]
                else:
                    all_data[lst_no, d] = MAX_VAL
                if lst_no == 0:
                    mat_list.append(lst[d][self.qp_header[idx_x.value]])
            lst_no += 1
        #all_data_rel = f_val(all_data)
        # for dim in range(num_tool):
        #     all_data_rel[dim, :] = all_data[dim, :]/all_data[0, :]
        # for v in range(n_groups):
        #     print(mat_list[v],',',all_data[0,v],',',all_data[1,v])

        ## Preparing the graph info
        # Other colors: https://matplotlib.org/gallery/color/named_colors.html
        colors = ('red', 'skyblue', 'indigo', 'olive', 'slateblue', 'magenta',
                  'slategray', 'limegreen', 'maroon', 'teal', 'khaki', 'purple',
                  'r', 'c', 'm', 'y', 'k')
        markers = ('s', 'o', '+', "."	,"," ,"v"	,"^"	,"<"	,">"	,
                   "1","2"	,"3"	,"4"	,"8"	,"s"	,"p"	,"P"	)
        width = 0.8
        ind = (num_tool + 1) * width * np.arange(len(mat_list))
        x_label = np.arange(1,len(mat_list)+1)
        font_axis = self.font0.copy()
        metric_label = 'solve time(sec)'
        fig = plt.figure(figsize=(100,5))
        ax1 = fig.add_subplot(111)
        ax1.set_yscale('log')
        for i in range(num_tool):
            #ax1.scatter(mat_list[:], all_data[i, :], s=10, c=colors[i],
            #            marker=markers[i], label=tool_name[i])
            ax1.bar(ind + i*width, all_data[i, :], width,
                   label=tool_name[i] + ' ' + metric_label, color=colors[i],
                   align='edge')
            #ax1.scatter(mat_list[:], all_data[1, :], s=10, c='r', marker="o", label=tool_name[1])

        #ax1.set_ylabel(mat_list[:], fontsize=12, fontproperties=font_axis)
        font_axis.set_weight('black')  # ['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']
        ax1.set_ylabel(y_label, fontsize=12, fontproperties=font_axis)
        ax1.set_xticks(ind)
        ax1.set_xticklabels(x_label, fontsize=5)
        #ax1.tick_params(axis='x', which='minor', direction='out', bottom=True)
        plt.xticks(rotation=55)
        plt.legend(loc='upper right')
        results_file = file_name
        print("Saving plots to %s" % results_file)
        plt.savefig(results_file, dpi=100)
        plt.show()

    def visualize_list(self, *list_args, idx_x, idx_y, y_label='',
                       file_name='', f_val, y_scale="linear"):
        """

        :param list_args: The list of reports as a set of dics
        :param idx_x: The index ID for x axis
        :param idx_y: Array of y-axis. It is a bar chart if
        len(idx_y) == 1 otherwise, that is a stacked bar chart.
        :param y_label: The label for y-axis, if it is not given,
        the label from the header array will be used.
        :param file_name: The file name to store the graph
        :param f_val: Doing an operation on the data, e.g., normalizing
        the data on one of the dimensions.
        :param y_scale: The scale of Y-axis e.g., log-scale, linear, ...
        :return: returns the raw plotted data.
        """
        list_args = list_args[0]
        if len(list_args) <= 0:
            print("Wrong lists to visualize.")
            return -1
        n_groups = len(list_args[0])
        tool_name = []
        for lst in list_args:
            if len(lst) != n_groups:
                print("Wrong lists to visualize.")
                return -1
            tool_name.append(lst[0][self.sympiler_header[self.LogHeader.TOOL.value]])
        num_tool = len(list_args)
        num_stacks = len(idx_y)
        all_data = np.zeros((num_tool,n_groups,num_stacks))
        all_data_rel = np.zeros((num_tool, n_groups))
        mat_list = []
        lst_no = 0
        for lst in list_args:
            for d in range(len(lst)):
                st_no = 0
                for st in idx_y:
                    all_data[lst_no, d, st_no] = lst[d][self.sympiler_header[st]]
                    st_no += 1
                if lst_no == 0:
                    mat_list.append(lst[d][self.sympiler_header[idx_x]])
            lst_no += 1

        all_data_rel = f_val(all_data)

        i = 0
        for tn in tool_name:
            if "MKL" in tn:
                break
            i = i + 1
        sup = compute_speedup_general(all_data, i)
        for i in range(num_tool):
            avg = np.average(sup[i, :])
            med = np.median(sup[i, :])
            max = np.max(sup[i, :])
            min = np.min(sup[i, :])
            faster = np.sum( sup[i, :] >= 1, axis=0)
            print(tool_name[i], " max: ", max, ": min", min, " : avg: ", avg, " median:", med, " faster: ", faster, " total: ", n_groups)

        all_data = sup
        # for dim in range(num_tool):
        #     all_data_rel[dim, :] = all_data[dim, :]/all_data[0, :]
        # for v in range(n_groups):
        #     print(mat_list[v],',',all_data[0,v],',',all_data[1,v])

        colors = ('red', 'skyblue', 'indigo', 'olive', 'slateblue', 'magenta',
                  'slategray', 'limegreen', 'maroon', 'teal', 'khaki', 'purple',
                  'r', 'c', 'm', 'y', 'k')
        markers = ('s', 'o', '+', ".", ",", "v", "^", "<", ">",
                   "1", "2", "3", "4", "8", "s", "p", "P")
        width = 0.3
        ind = (num_tool + 1) * width * np.arange(len(mat_list))
        x_label = np.arange(1, len(mat_list) + 1)
        font_axis = self.font0.copy()
        metric_label = 'solve time(sec)'
        fig = plt.figure(figsize=(100, 5))
        ax1 = fig.add_subplot(111)
        ax1.set_yscale('log')
        for i in range(num_tool):
            # ax1.scatter(mat_list[:], all_data[i, :], s=10, c=colors[i],
            #            marker=markers[i], label=tool_name[i])
#            ttmmpp = np.ravel(all_data[i, :, :].T)
            ttmmpp = all_data[i, :]
            ax1.bar(ind + i * width, ttmmpp, width,
                    label=tool_name[i] + ' ' + metric_label, color=colors[i],
                    align='edge')
            # ax1.scatter(mat_list[:], all_data[1, :], s=10, c='r', marker="o", label=tool_name[1])

        # ax1.set_ylabel(mat_list[:], fontsize=12, fontproperties=font_axis)
        font_axis.set_weight('black')  # ['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']
        ax1.set_ylabel(y_label, fontsize=12, fontproperties=font_axis)
        ax1.set_xticks(ind)
        ax1.set_xticklabels(x_label, fontsize=5)
        # ax1.tick_params(axis='x', which='minor', direction='out', bottom=True)
        plt.xticks(rotation=55)
        plt.legend(loc='upper right')
        results_file = file_name
        print("Saving plots to %s" % results_file)
        plt.savefig(results_file, dpi=100)
        plt.show()

    def write_to_csv(self, report, csv_file):
        csv_columns = []
        if len(report) > 0:
            one_row = report[0]
        for key, value in one_row.items():
            csv_columns.append(key)
        def remove_key(d, key):
            r = dict(d)
            del r[key]
            return r
        try:
            with open(csv_file, 'w') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
                writer.writeheader()
                for data in report:
                    #data_p = remove_key(data, None)
                    writer.writerow(data)
        except IOError:
            print("I/O error")


def median_idx(tmp_sort):
    med_idx = np.argsort(tmp_sort)[len(tmp_sort) / 2]
    return med_idx


def min_idx(tmp_sort):
    min_idx = np.argsort(tmp_sort)[0]
    return min_idx


def scaling_array(in_list):
    scaled = np.zeros(in_list.shape)
    for l in range(len(in_list)):
        scaled[l, :] = in_list[l, :]/in_list[0, :]
    return scaled


def nop(in_list):
    return in_list




class reporting_type(Enum):
    LinSolve1 = 0
    LinSolve2 = 1
    SCATQP = 2
    BARQP = 3
    PerfProf = 4
    AllFile = 5


def main1(argv):
    output_path = ''
    input_path1 = ''
    input_path2 = ''
    input_path3 = ''
    input_path4 = ''
    report_type = 1
    del_items = ["bunny2_bunny2", "wolf_falling6", "knight_falling2",
                 "InjBnd_Mapping0", "InjBnd_Mapping1",
                 "InjBnd_Mapping2", "InjBnd_Mapping3",
                 # "InjBnd_MappingImp0","InjBnd_MappingImp1","InjBnd_MappingImp2",
                 "lamb_falling1", "horse_horse5", "horse_horse1",
                 "horse_falling1", "garg3_floor1", "bunny_falling1", "bunny2_falling1", "bunny_bunny2",
                 "brick0_brick", "screwdriver_falling"]
    for i in range(1, 130):
        item = "Beam_floor" + str(i)
        del_items.append(item)
    for i in range(290, 338):
        item = "Beam_floor" + str(i)
        del_items.append(item)
    for i in range(400, 524):
        item = "Cube_floor" + str(i)
        # del_items.append(item)
    for i in range(4, 12):
        item = "lion" + str(i) + "_s"
        del_items.append(item)
    for i in range(5, 800):
        item = "Bar2k_floor" + str(i)
        del_items.append(item)
    for i in range(10, 19):
        item = "nail" + str(i) + "t"
        del_items.append(item)
    try:
        opts , args = getopt.getopt(argv,"hi:f:s:w:d:o:t:",["help","inputPath=","outputPath=","fileName=","inputPath4=",
                                                            "directory=", "outputPath=", "reportType="])
    except getopt.GetoptError:
        print('visualize.py -i <sympiler report path> -s <second sym report> -f <mkl report path> -o <output Path for graphs> ')
        sys.exit(-1)
    for opt, arg in opts:
        if opt == '-h':
            print('visualize.py -i <sympiler report path> -f <mkl report path> -o <output Path for graphs> ')
        elif opt in ("-i", "--inputPath"):
            input_path1 = arg
        elif opt in ("-f", "--fileName"):
            input_path2 = arg.strip()
        elif opt in ("-s", "--sym2"):
            input_path3 = arg.strip()
        elif opt in ("-w", "--inputPath4"):
            input_path4 = arg.strip()
        elif opt in ("-d", "--directory"):
            input_dir = arg.strip()
        elif opt in ("-o", "--outputPath"):
            output_path = arg
        elif opt in ("-t", "--reportType"):
            report_type = int(arg)

    graphGenerator = Visualizer(output_path)
    eq_idx = graphGenerator.LogHeader.NNZ.value
    key_item = graphGenerator.sympiler_header[Visualizer.LogHeader.MATN.value]
    print('Reading CSV files from the input directory: ', input_path1)
    #report_type = 2
    def condition_diag(dic_in):
        str_dim = graphGenerator.sympiler_header[Visualizer.LogHeader.ROW.value]
        str_size = graphGenerator.sympiler_header[Visualizer.LogHeader.NNZ.value]
        p_dim = int(dic_in[str_dim])
        p_size = int(dic_in[str_size])
        if p_size > p_dim or p_size < 50:
            return True
        else:
            return False
    if report_type == reporting_type.LinSolve2.value:
        # reading QP time with no initial time
        eps = 'eps_abs'
        th = '1e-06'
        splitted = []
        tmp_qp_time = graphGenerator.read_csv(input_path1)
        df = pd.DataFrame(tmp_qp_time)
        gb = df.groupby([eps])
        gbg = [gb.get_group(x) for x in gb.groups]
        sw = False
        for g in gbg:
            for gg in g[eps]:
                if gg != th:
                    sw = False
                else:
                    sw = True
                break
            if sw:
                tmp = []
                for idx in g.index:
                    tmp.append(tmp_qp_time[idx])
                tmp_qp_time = tmp
                break
            else:
                continue
        qp_time = []
        qp_init_time = []
        lst_len = len(tmp_qp_time)
        frac_time = np.zeros(lst_len)
        for it in range(lst_len):
            tmp_qp_time[it]['Problem Name'] = get_name(tmp_qp_time[it]['Problem Name'])
            tmp_qp_time[it]['Time (s)'] = float(tmp_qp_time[it]['Time (s)'])
            tmp_qp_time[it]['Init Time (s)'] = float(tmp_qp_time[it]['Init Time (s)'])
            tmp_qp_time[it]['Iter Time (s)'] = tmp_qp_time[it]['Time (s)'] - tmp_qp_time[it]['Init Time (s)']
            frac_time[it] = tmp_qp_time[it]['Init Time (s)'] / tmp_qp_time[it]['Time (s)']

        #print(th, ": Max: ", np.max(frac_time), " Min: ", np.min(frac_time), " Avg: ", np.average(frac_time), " Med: ", np.median(frac_time))
        files = [f for f in glob.glob(input_dir + "*.csv", recursive=False)]
        if len(files) <= 0:
            print("No csv file found in the path")
            sys.exit(-1)
        all_logs = []
        all_logs_red = []
        tool_names = []
        for f in files:
            tmp = graphGenerator.read_csv(f)
            tmp = graphGenerator.convert_log(tmp)
            all_logs.append(tmp)



        for l in all_logs:
            log_tmp = exclude_items_general(l, condition_diag)
            for del_item in del_items:
                log_tmp = exclude_items(log_tmp, del_item, key_item)
            log_tmp = sorted(log_tmp, key=itemgetter(graphGenerator.sympiler_header[eq_idx]))
            all_logs_red.append(log_tmp)

        for l in range(len(all_logs_red)):
            for k in range(len(all_logs_red[l])):
                t_name = all_logs_red[l][k]['Matrix Name']
                found_item = list(filter(lambda person: person['Problem Name'] == t_name, tmp_qp_time))
                if len(found_item) == 0:
                    print("Errrr")
                time_iter = found_item[0]['Iter Time (s)']
                all_logs_red[l][k]['QP Time'] = time_iter + all_logs_red[l][k]['Factorization Time (sec)']


        # Writing the CSV file if necessary
        # graphGenerator.write_to_csv(mkl_log,'out_csv_mkl.csv')
        # graphGenerator.write_to_csv(sym_log4, 'out_csv_sym.csv')
        print('Drawing a graph')
        # idx_y_val = np.array([
        #     Visualizer.LogHeader.ANAT.value,
        #     Visualizer.LogHeader.FACT.value,
        #     Visualizer.LogHeader.SOLT.value,
        #     Visualizer.LogHeader.ITERT.value
        # ])


        idx_y_val = np.array([Visualizer.LogHeader.FACT.value])
        # idx_y_val = np.array([Visualizer.LogHeader.QPTT.value])
        graphGenerator.visualize_list(all_logs_red,
                                      idx_x=Visualizer.LogHeader.MATN.value,
                                      idx_y=idx_y_val,
                                      #y_label="MKL/ParSy Fact time",
                                      file_name="accuracy1",
                                      #y_label="MKL Pardiso/ ParSy LDL Time",
                                      f_val=nop
                                      #, y_scale="log"
                                      )
    sys.exit(0)

def main(argv):
    output_path = ''
    input_path1 = ''
    input_path2 = ''
    input_path3 = ''
    input_path4 = ''
    report_type = 1
    try:
        opts , args = getopt.getopt(argv,"hi:f:s:w:d:o:t:",["help","inputPath=","outputPath=","fileName=","inputPath4=",
                                                            "directory=", "outputPath=", "reportType="])
    except getopt.GetoptError:
        print('visualize.py -i <sympiler report path> -s <second sym report> -f <mkl report path> -o <output Path for graphs> ')
        sys.exit(-1)
    for opt, arg in opts:
        if opt == '-h':
            print('visualize.py -i <sympiler report path> -f <mkl report path> -o <output Path for graphs> ')
        elif opt in ("-i", "--inputPath"):
            input_path1 = arg
        elif opt in ("-f", "--fileName"):
            input_path2 = arg.strip()
        elif opt in ("-s", "--sym2"):
            input_path3 = arg.strip()
        elif opt in ("-w", "--inputPath4"):
            input_path4 = arg.strip()
        elif opt in ("-d", "--directory"):
            input_dir = arg.strip()
        elif opt in ("-o", "--outputPath"):
            output_path = arg
        elif opt in ("-t", "--reportType"):
            report_type = int(arg)

    files = [f for f in glob.glob(input_dir + "*.csv", recursive=False)]
    for f in files:
        if 'ql' in f or 'quadprog' in f:
            continue
        logs_tmp = graphGenerator.read_qp_log_split(f, sorting_idx,
                                                    Visualizer.qp_log_header.STIME,
                                                    tol, tol_str, problem_group)

if __name__ == "__main__":
   main(sys.argv[1:])