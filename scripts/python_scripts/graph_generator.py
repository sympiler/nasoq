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
    exclude_items,compute_speedup,exclude_items_general

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
                                'Factorization (GFLOPS)'
                                ]
        self.qp_header = ['Tool Name',
                           'Problem Name', 'Time (s)', 'Residual Lagrangian inf',
                          'Constraint Satisfaction Inf',
                          'Status', 'Problem Size', 'Hessian NNZ', 'Eq Constraint NNZ',
                          'Ineq Constraint NNZ', 'Outer GMRES Iter', 'GMRES Tol', 'Diagonal Pert', 'eps_abs',
                          'Inner GMRES Iter', 'Hessian dim','Non-negativity Inf','Complementarity Inf',
                          'Problem Type', 'Hessian dim', '# of Ineq Const', '# of Eq Constraints'
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

    class qp_log_header(Enum):
        TOOLN = 0
        PROBN = 1
        STIME = 2
        LRNRM = 3
        CSNRM = 4
        STATS = 5
        PROBS = 6
        HSSNZ = 7
        CNSNZ = 8
        INQNZ = 9
        INIR = 10
        TOLIR = 11
        PERTB = 12
        EPSAB = 13
        OUTIR = 14
        HSSDM = 15
        NGINF = 16
        CMINF = 17
        CLSTY = 18
        PROBD = 19
        CNSID = 20
        CNSED = 21
        INVAL = 100



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

    def fact_time_sympiler(self, sym_report_path, cmp_idx=LogHeader.FACT.value):
        mat_list = []
        sym_report = self.read_csv(sym_report_path)
        i1 = 0
        sw = True
        nxt_mat = ''
        while i1 < len(sym_report):
            sr = sym_report[i1]
            cur_mat = sr[self.sympiler_header[self.LogHeader.MATN.value]]  # matrix name
            i2 = 0
            if i1 + 1 < len(sym_report):
                nxt_mat = sym_report[i1 + 1][self.sympiler_header[self.LogHeader.MATN.value]]
            else:
                sw = False
            min_time = sys.maxsize
            min_idx = 0
            while cur_mat == nxt_mat and sw:
                cur_time = float(nxt_mat[self.sympiler_header[cmp_idx]])
                if min_time > cur_time:
                    min_time = cur_time
                    min_idx = i2
                i2 += 1
            mat_list.append(sym_report[i1 + min_idx])
            i1 += (i2 + 1)
        return mat_list

    def fact_time_mkl(self, mkl_report_path):
        mat_list = []
        new_dic = []
        file_token = 'Matrix Name***'
        time_afs = 'Total time spent'
        time_solve = 'Time spent in direct solver at solve step (solve)'
        extra_info_s = '+++extrainfo:'
        sup_node = 'number of supernodes:'
        iter_refine = 'Time spent in additional calculations'
        file_found = False
        analysis_time = False
        fact_time = False
        solve_time = False
        iter_ref_Time = False
        super_node = False
        extra_info = False
        #log = os.popen(mkl_report_path)
        f = open(mkl_report_path, "r")
        log = f.read()
        totoal_time = 0
        for line in log.splitlines(True):
            if not file_found:
                idx_f = line.find(file_token)
                if idx_f >= 0:
                    file_found = True
                    totoal_time = 0
                    new_dic = self.empty_dic.copy()
                    new_dic[self.sympiler_header[self.LogHeader.TOOL.value]] = 'MKL Pardiso'
                    mat_full_name = line.split()[2].split('/')
                    mat_name = mat_full_name[len(mat_full_name)-1].split('.')[0]
                    new_dic[self.sympiler_header[self.LogHeader.MATN.value]] = mat_name
                    continue
            if file_found and not analysis_time:
                idx_f = line.find(time_afs)
                if idx_f >= 0:
                    analysis_time = True
                    a_time = line.split()[4]
                    totoal_time += float(a_time)
                    new_dic[self.sympiler_header[self.LogHeader.ANAT.value]] = a_time
                    continue
            if analysis_time and not fact_time:
                idx_f = line.find(time_afs)
                if idx_f >= 0:
                    fact_time = True
                    a_time = line.split()[4]
                    totoal_time += float(a_time)
                    new_dic[self.sympiler_header[self.LogHeader.FACT.value]] = a_time
                    continue
            if fact_time and not solve_time:
                idx_f = line.find(time_solve)
                if idx_f >= 0:
                    solve_time = True
                    a_time = line.split()[10]
                    totoal_time += float(a_time)
                    new_dic[self.sympiler_header[self.LogHeader.SOLT.value]] = a_time
                    continue
            if solve_time and not iter_ref_Time:
                idx_f = line.find(iter_refine)
                if idx_f >= 0:
                    iter_ref_Time = True
                    a_time = line.split()[6]
                    totoal_time += float(a_time)
                    new_dic[self.sympiler_header[self.LogHeader.ITERT.value]] = a_time
                    continue
            if iter_ref_Time and not super_node:
                idx_f = line.find(sup_node)
                if idx_f >= 0:
                    super_node = True
                    tmp_buf = line.split()[3]
                    new_dic[self.sympiler_header[self.LogHeader.SUPN.value]] = tmp_buf
                    continue
            if super_node and not extra_info:
                idx_f = line.find(extra_info_s)
                if idx_f >= 0:
                    file_found = False
                    analysis_time = False
                    fact_time = False
                    solve_time = False
                    super_node = False
                    iter_ref_Time = False
                    super_node = False
                    tmp_buf = line.split(':')[1].split(',')
                    totoal_time += float(tmp_buf[0])
                    new_dic[self.sympiler_header[self.LogHeader.TOTT.value]] = str(totoal_time)
                    new_dic[self.sympiler_header[self.LogHeader.ORDT.value]] = tmp_buf[0]
                    new_dic[self.sympiler_header[self.LogHeader.NNZL.value]] = tmp_buf[1]
                    new_dic[self.sympiler_header[self.LogHeader.FLPS.value]] = tmp_buf[2]
                    new_dic[self.sympiler_header[self.LogHeader.ROW.value]] = tmp_buf[3]
                    new_dic[self.sympiler_header[self.LogHeader.NREF.value]] = tmp_buf[4]
                    new_dic[self.sympiler_header[self.LogHeader.NNZ.value]] = tmp_buf[5]
                    new_dic[self.sympiler_header[self.LogHeader.THRN.value]] = tmp_buf[6]
                    new_dic[self.sympiler_header[self.LogHeader.NRM1.value]] = tmp_buf[7]
                    new_dic[self.sympiler_header[self.LogHeader.NRM2.value]] = tmp_buf[8]
                    new_dic[self.sympiler_header[self.LogHeader.NRMI.value]] = tmp_buf[9]
                    new_dic[self.sympiler_header[self.LogHeader.NR1A.value]] = tmp_buf[7]
                    new_dic[self.sympiler_header[self.LogHeader.NR2A.value]] = tmp_buf[8]
                    new_dic[self.sympiler_header[self.LogHeader.NRIA.value]] = tmp_buf[9]
                    new_dic[self.sympiler_header[self.LogHeader.TRGT.value]] = tmp_buf[10]
                    #updating analysis time
                    anal_time = float(new_dic[self.sympiler_header[self.LogHeader.ANAT.value]])
                    ord_time = float(new_dic[self.sympiler_header[self.LogHeader.ORDT.value]])
                    total_anal_time = anal_time + ord_time
                    new_dic[self.sympiler_header[self.LogHeader.ANAT.value]] = str(total_anal_time)
                    mat_list.append(new_dic)
        mat_list_sorte = sorted(mat_list,
                                key=itemgetter(self.sympiler_header[Visualizer.LogHeader.MATN.value]))

        return mat_list_sorte

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

        return red_sort_log


    def read_qp_log(self, report_path, idx_eq, idx_min=qp_log_header.INVAL, tol=1e-12):
        sym_log = self.read_csv(report_path)
        red_log = self.qp_lib_list_reduction(sym_log,idx_eq,idx_min,tol)
        #eee_tmp = self.qp_lib_list_reduction_2(sym_log,idx_eq,idx_min,tol)
        #return eee_tmp
        # extracting matrix name from the path
        for i in range(len(red_log)):
            full_name = red_log[i][self.qp_header[self.qp_log_header.PROBN.value]]
            mat_full_name = full_name.split('/')
            mat_name = mat_full_name[len(mat_full_name) - 1].split('.')[0]
            #mat_name = mat_name.split('_')
            us_pos = mat_name.rfind('_')
            if us_pos>0:
                mat_name = mat_name[0:us_pos]
            red_log[i][self.qp_header[self.qp_log_header.PROBN.value]] = mat_name
            psize = 0
            if not red_log[i][self.qp_header[self.qp_log_header.HSSNZ.value]] == INV:
                psize = int(red_log[i][self.qp_header[self.qp_log_header.HSSNZ.value]])
                #FIXME: remove this after the logs are fixed
                if "NASOQ" in red_log[i][self.qp_header[self.qp_log_header.TOOLN.value]]:
                    dim = int(red_log[i][self.qp_header[self.qp_log_header.HSSDM.value]])
                    psize = 2*psize - dim
            if not red_log[i][self.qp_header[self.qp_log_header.CNSNZ.value]] == INV:
                psize += int(red_log[i][self.qp_header[self.qp_log_header.CNSNZ.value]])
            if not red_log[i][self.qp_header[self.qp_log_header.INQNZ.value]] == INV:
                psize += int(red_log[i][self.qp_header[self.qp_log_header.INQNZ.value]])
            red_log[i][self.qp_header[self.qp_log_header.PROBS.value]] = psize

        red_sort_log = sorted(red_log, key=itemgetter(self.qp_header[idx_eq.value]))
        return red_sort_log

    def qp_lib_list_reduction(self,list_log, idx_eq,idx_min, tol):
        min_time = 1e8
        red_sort_log = sorted(list_log, key=itemgetter(self.qp_header[idx_eq.value]))
        reduced_list = []
        if len(list_log) <= 1:
            return list_log
        tmp_list = []
        tmp_sort = np.array([])
        i = 0
        r_i = 0 # num of reduced elements
        while i < len(red_sort_log):
            num_rep = 0
            cur_prob = red_sort_log[i][self.qp_header[idx_eq.value]]
            sat_norm = float(red_sort_log[i][self.qp_header[self.qp_log_header.CSNRM.value]])
            lag_norm = float(red_sort_log[i][self.qp_header[self.qp_log_header.LRNRM.value]])
            comp_norm = float(red_sort_log[i][self.qp_header[self.qp_log_header.CMINF.value]])
            nng_norm = float(red_sort_log[i][self.qp_header[self.qp_log_header.NGINF.value]])
            min_time = float(red_sort_log[i][self.qp_header[idx_min.value]])
            status = red_sort_log[i][self.qp_header[self.qp_log_header.STATS.value]]
            if sat_norm < tol and lag_norm < tol and min_time < MAX_VAL and nng_norm < tol  \
                    and (status == 'solved' or status == '1' or status == 'The optimality conditions are satisfied'):
                red_sort_log[i][self.qp_header[self.qp_log_header.STATS.value]] = OPTIMAL
            else:
                min_time = 1e8
                red_sort_log[i][self.qp_header[self.qp_log_header.STATS.value]] = 'Not converged'
            tmp_list.append(red_sort_log[i])
            if i+1 < len(red_sort_log):
                nxt_prob = red_sort_log[i+1][self.qp_header[idx_eq.value]]
            else:
                break;
            i = i+1
            while cur_prob == nxt_prob:
                nxt_time = float(red_sort_log[i][self.qp_header[idx_min.value]])
                nxt_sat_norm = float(red_sort_log[i][self.qp_header[self.qp_log_header.CSNRM.value]])
                nxt_lag_norm = float(red_sort_log[i][self.qp_header[self.qp_log_header.LRNRM.value]])
                if sat_norm < tol and lag_norm < tol and min_time < MAX_VAL and nng_norm < tol  \
                    and (status == 'solved' or status == '1'):
                    min_time = nxt_time
                    red_sort_log[i][self.qp_header[self.qp_log_header.STATS.value]] = OPTIMAL
                    tmp_list[r_i] = red_sort_log[i]
                i = i+1
                if i < len(red_sort_log):
                    nxt_prob = red_sort_log[i][self.qp_header[idx_eq.value]]
                else:
                    break
            r_i = r_i + 1
        return tmp_list


    def compute_qp_size(self, red_log):
        psize = 0
        if not red_log[self.qp_header[self.qp_log_header.HSSNZ.value]] == INV:
            psize = int(red_log[self.qp_header[self.qp_log_header.HSSNZ.value]])
        if not red_log[self.qp_header[self.qp_log_header.CNSNZ.value]] == INV:
            psize += int(red_log[self.qp_header[self.qp_log_header.CNSNZ.value]])
        if not red_log[self.qp_header[self.qp_log_header.INQNZ.value]] == INV:
            psize += int(red_log[self.qp_header[self.qp_log_header.INQNZ.value]])
        red_log[self.qp_header[self.qp_log_header.PROBS.value]] = psize
        return psize


    def find_problem_size(self, report_path, idx_eq, idx_min=qp_log_header.INVAL,
                          tol=1e-6, tol_str='1e-06', class_type=["All"]):
        prt_str = self.qp_header[self.qp_log_header.PERTB.value]
        ir_in = self.qp_header[self.qp_log_header.INIR.value]
        ir_out = self.qp_header[self.qp_log_header.OUTIR.value]
        ir_tol = self.qp_header[self.qp_log_header.TOLIR.value]

        sym_log = self.read_csv(report_path)
        #        tn = sym_log['Tool Name']
        splited_log = []
        if len(class_type) == 1:
            if class_type[0] == "" or class_type[0] == "All":
                splited_log = self.qp_lib_list_split_byname(sym_log)
            else:
                splited_log = self.qp_lib_list_split_bycat(sym_log, class_type)
        else:
            splited_log = self.qp_lib_list_split_bycat(sym_log, class_type)
        #        if tn == 'NASOQ-auto':
        #            splited_log = sym_log
        red_tol_log = []
        for lst in splited_log:
            tt = lst[0][self.qp_header[self.qp_log_header.TOOLN.value]]
            if (lst[0]['eps_abs'] == tol_str and not tt == 'QL') or (lst[0]['eps_abs'] == tol_str and tt == 'QL'):
                lst1 = self.qp_lib_list_reduction(lst, idx_eq, idx_min, tol)
                for l in lst1:
                    tool_name = l[self.qp_header[self.qp_log_header.TOOLN.value]]
                    prt_str_val = l[prt_str]
                    ir_in_val = l[ir_in]
                    ir_out_val = l[ir_out]
                    ir_tol_val = l[ir_tol]
                    # tn = tool_name + '_' + prt_str_val + '_' + ir_out_val + '_' + ir_in_val + '_' + ir_tol_val
                    # tn = tn.split('/')[0]
                    tn = tool_name
                    l[self.qp_header[self.qp_log_header.TOOLN.value]] = tn
                red_tol_log.append(lst1)

        for red_log in red_tol_log:
            if len(red_log) > 0:
                if "OSQP-polished" not in red_log[0][self.qp_header[self.qp_log_header.TOOLN.value]]:
                    continue
            for i in range(len(red_log)):
                if "OSQP-polished" in red_log[i][self.qp_header[self.qp_log_header.TOOLN.value]]:
                    full_name = red_log[i][self.qp_header[self.qp_log_header.PROBN.value]]
                    mat_full_name = full_name.split('/')
                    mat_name = mat_full_name[len(mat_full_name) - 1].split('.')[0]
                    # mat_name = mat_name.split('_')
                    us_pos = mat_name.rfind('_')
                    if us_pos > 0:
                        mat_name = mat_name[0:us_pos]
                    self.qp2size[mat_name] = self.compute_qp_size(red_log[i])



    def read_qp_log_split(self, report_path, idx_eq, idx_min=qp_log_header.INVAL,
                          tol=1e-6, tol_str='1e-06', class_type=["All"]):
        prt_str = self.qp_header[self.qp_log_header.PERTB.value]
        ir_in = self.qp_header[self.qp_log_header.INIR.value]
        ir_out = self.qp_header[self.qp_log_header.OUTIR.value]
        ir_tol = self.qp_header[self.qp_log_header.TOLIR.value]

        sym_log = self.read_csv(report_path)
#        tn = sym_log['Tool Name']
        splited_log = []
        if len(class_type)==1:
            if class_type[0] == "" or class_type[0] == "All":
                splited_log = self.qp_lib_list_split_byname(sym_log)
            else:
                splited_log = self.qp_lib_list_split_bycat(sym_log, class_type)
        else:
            splited_log = self.qp_lib_list_split_bycat(sym_log,class_type)
#        if tn == 'NASOQ-auto':
#            splited_log = sym_log
        red_tol_log = []
        for lst in splited_log:
            tt = lst[0][self.qp_header[self.qp_log_header.TOOLN.value]]
            if (lst[0]['eps_abs'] == tol_str and not tt == 'QL') or (lst[0]['eps_abs'] == tol_str and tt == 'QL'):
                lst1 = self.qp_lib_list_reduction(lst, idx_eq, idx_min, tol)
                for l in lst1:
                    tool_name = l[self.qp_header[self.qp_log_header.TOOLN.value]]
                    prt_str_val = l[prt_str]
                    ir_in_val = l[ir_in]
                    ir_out_val = l[ir_out]
                    ir_tol_val = l[ir_tol]
                    #tn = tool_name + '_' + prt_str_val + '_' + ir_out_val + '_' + ir_in_val + '_' + ir_tol_val
                    #tn = tn.split('/')[0]
                    tn = tool_name
                    l[self.qp_header[self.qp_log_header.TOOLN.value]] = tn
                red_tol_log.append(lst1)

        # extracting matrix name from the path and seting size

        lst_i = 0
        for red_log in red_tol_log:
            for i in range(len(red_log)):
                full_name = red_log[i][self.qp_header[self.qp_log_header.PROBN.value]]
                mat_full_name = full_name.split('/')
                mat_name = mat_full_name[len(mat_full_name) - 1].split('.')[0]
                #mat_name = mat_name.split('_')
                us_pos = mat_name.rfind('_')
                if us_pos>0:
                    mat_name = mat_name[0:us_pos]
                red_log[i][self.qp_header[self.qp_log_header.PROBN.value]] = mat_name
                psize = self.qp2size[mat_name]
                red_log[i][self.qp_header[self.qp_log_header.PROBS.value]] = psize
                # if not red_log[i][self.qp_header[self.qp_log_header.HSSNZ.value]] == INV:
                #     psize = int(red_log[i][self.qp_header[self.qp_log_header.HSSNZ.value]])
                #     # FIXME: remove this after the logs are fixed
                #     if "NASOQ" in red_log[i][self.qp_header[self.qp_log_header.TOOLN.value]] or \
                #             "MOSEK" in red_log[i][self.qp_header[self.qp_log_header.TOOLN.value]] or \
                #             "OSQP" in red_log[i][self.qp_header[self.qp_log_header.TOOLN.value]] or \
                #             "QL" in red_log[i][self.qp_header[self.qp_log_header.TOOLN.value]]:
                #         dim = int(red_log[i][self.qp_header[self.qp_log_header.HSSDM.value]])
                #         psize = 2 * psize - dim
                # if not red_log[i][self.qp_header[self.qp_log_header.CNSNZ.value]] == INV:
                #     psize += int(red_log[i][self.qp_header[self.qp_log_header.CNSNZ.value]])
                # if not red_log[i][self.qp_header[self.qp_log_header.INQNZ.value]] == INV:
                #     psize += int(red_log[i][self.qp_header[self.qp_log_header.INQNZ.value]])
                # red_log[i][self.qp_header[self.qp_log_header.PROBS.value]] = psize

            red_tol_log[lst_i] = sorted(red_log, key=itemgetter(self.qp_header[idx_eq.value]))
            lst_i = lst_i+1
        return red_tol_log

    def qp_lib_list_split(self, list_log):
        splitted = []
        prt_str = self.qp_header[self.qp_log_header.PERTB.value]
        eps = self.qp_header[self.qp_log_header.EPSAB.value]
        ir_in = self.qp_header[self.qp_log_header.INIR.value]
        ir_out = self.qp_header[self.qp_log_header.OUTIR.value]
        ir_tol = self.qp_header[self.qp_log_header.TOLIR.value]
        df = pd.DataFrame(list_log)
        gb = df.groupby([eps, prt_str, ir_in, ir_out, ir_tol])
        gbg = [gb.get_group(x) for x in gb.groups]
        for g in gbg:
            tmp = []
            for idx in g.index:
                tmp.append(list_log[idx])
            splitted.append(tmp)
        return splitted

    def qp_lib_list_split_byname(self, list_log):
        splitted = []
        eps = self.qp_header[self.qp_log_header.EPSAB.value]
        tool_str = self.qp_header[self.qp_log_header.TOOLN.value]
        df = pd.DataFrame(list_log)
        gb = df.groupby([eps, tool_str])
        gbg = [gb.get_group(x) for x in gb.groups]
        for g in gbg:
            tmp = []
            for idx in g.index:
                tmp.append(list_log[idx])
            splitted.append(tmp)
        return splitted

    def qp_lib_list_split_bycat(self, list_log, class_type):
        splitted = []
        reduced_logs = []
        eps = self.qp_header[self.qp_log_header.EPSAB.value]
        tool_str = self.qp_header[self.qp_log_header.TOOLN.value]
        prob_type = self.qp_header[self.qp_log_header.CLSTY.value]

        for s in list_log:
            tmp = s[prob_type]
            for ct in class_type:
                if tmp == ct:
                    reduced_logs.append(s)

        df = pd.DataFrame(reduced_logs)
        gb = df.groupby([eps, tool_str])
        gbg = [gb.get_group(x) for x in gb.groups]
        for g in gbg:
            tmp = []
            for idx in g.index:
                tmp.append(reduced_logs[idx])
            splitted.append(tmp)
        return splitted

    def qp_lib_list_reduction_2(self,list_log, idx_eq,idx_min, tol):
        min_time = 1e8
        red_sort_log = sorted(list_log, key=itemgetter(self.qp_header[idx_eq.value]))
        reduced_list = []
        if len(list_log) <= 1:
            return list_log
        tmp_list = []
        tmp_sort = np.array([])
        i = 0
        r_i = 0 # num of reduced elements
        while i < len(red_sort_log):
            num_rep = 0
            cur_prob = red_sort_log[i][self.qp_header[idx_eq.value]]
            sat_norm = float(red_sort_log[i][self.qp_header[self.qp_log_header.CSNRM.value]])
            lag_norm = float(red_sort_log[i][self.qp_header[self.qp_log_header.LRNRM.value]])
            min_time = float(red_sort_log[i][self.qp_header[idx_min.value]])
            psize = int(red_sort_log[i][self.qp_header[self.qp_log_header.HSSNZ.value]])
            psize = psize + int(red_sort_log[i][self.qp_header[self.qp_log_header.INQNZ.value]])
            psize = psize + int(red_sort_log[i][self.qp_header[self.qp_log_header.CNSNZ.value]])

            if sat_norm < tol and lag_norm < tol and min_time < 1e3:
                red_sort_log[i][self.qp_header[self.qp_log_header.STATS.value]] = 'optimal'
            else:
                min_time = 1e8
                red_sort_log[i][self.qp_header[self.qp_log_header.STATS.value]] = 'Not converged'
            tmp_list.append(red_sort_log[i])
            if i+1 < len(red_sort_log):
                nxt_prob = red_sort_log[i+1][self.qp_header[idx_eq.value]]
            else:
                break;
            i = i+1
            while cur_prob == nxt_prob:
                nxt_time = float(red_sort_log[i][self.qp_header[idx_min.value]])
                nxt_sat_norm = float(red_sort_log[i][self.qp_header[self.qp_log_header.CSNRM.value]])
                nxt_lag_norm = float(red_sort_log[i][self.qp_header[self.qp_log_header.LRNRM.value]])
                psize = int(red_sort_log[i][self.qp_header[self.qp_log_header.HSSNZ.value]])
                psize = psize + int(red_sort_log[i][self.qp_header[self.qp_log_header.INQNZ.value]])
                psize = psize + int(red_sort_log[i][self.qp_header[self.qp_log_header.CNSNZ.value]])
                #if nxt_sat_norm < tol and nxt_lag_norm < tol and nxt_time < 1000:
                if nxt_sat_norm < tol and nxt_lag_norm < tol and nxt_time < MAX_VAL :
                    min_time = nxt_time
                    red_sort_log[i][self.qp_header[self.qp_log_header.STATS.value]] = 'optimal'
                    tmp_list[r_i] = red_sort_log[i]
                i = i+1
                if i < len(red_sort_log):
                    nxt_prob = red_sort_log[i][self.qp_header[idx_eq.value]]
                else:
                    break
            r_i = r_i + 1
        rep = np.zeros((100, 17, 13))
        for irn in {0,1,2,3,4,5,6,7,8,9, 99}:
            for tol in range(15, 16):
                for pert in range(10, 11):
                    cur_tol = pow(10, -tol)
                    cur_pert = pow(10, -pert)
                    for i in range(len(red_sort_log)):
                        iter_noi = int(red_sort_log[i][self.qp_header[self.qp_log_header.INIR.value]])
                        diag_perti = float(red_sort_log[i][self.qp_header[self.qp_log_header.PERTB.value]])
                        toli = float(red_sort_log[i][self.qp_header[self.qp_log_header.TOLIR.value]])
                        statui = red_sort_log[i][self.qp_header[self.qp_log_header.STATS.value]]
                        if statui == "optimal" and iter_noi == irn and diag_perti == cur_pert and toli == cur_tol:
                            rep[irn, tol, pert] = rep[irn, tol, pert] + 1



        # for irn in {0,1,2,3,4,5,6,7,8,9, 99}:
        #     for tol in range(15, 16):
        #         for pert in range(7, 12):
        #             print(irn , ", ", tol, ", ", pert, ", ", rep[irn, tol, pert])

        rep2 = np.zeros(13)
        for pert in range(1, 12):
            cur_pert = pow(10, -pert)
            for i in range(len(red_sort_log)):
                diag_perti = float(red_sort_log[i][self.qp_header[self.qp_log_header.PERTB.value]])
                statui = red_sort_log[i][self.qp_header[self.qp_log_header.STATS.value]]
                if statui == "optimal"  and diag_perti == cur_pert:
                    rep2[ pert] = rep2[pert] + 1
        for pert in range(7, 12):
            print( pert, ", ", rep2[ pert])
        return red_sort_log

    def qp_log_export(self, log_list, problems_type):
        num_converged = 0
        avg_solve = 0
        perf_profile = []
        tool_name = log_list[0][self.qp_header[self.qp_log_header.TOOLN.value]]
        # file_name = os.path.join(".", "results", problems_type,tool_name,
        #                             "statistics.txt")
        file_name = "results_" + problems_type +"_" + tool_name + "_statistics.txt"
        f = open(file_name, "w")
        f.write(self.qp_header[self.qp_log_header.PROBN.value])
        f.write(",")
        f.write(self.qp_header[self.qp_log_header.TOOLN.value])
        f.write(",run_time")
        #f.write(self.qp_header[self.qp_log_header.STIME.value])
        f.write(",status")
        #f.write(self.qp_header[self.qp_log_header.STATS.value])
        f.write("\n")
        for el in log_list:
            f.write(el[self.qp_header[self.qp_log_header.PROBN.value]])
            f.write(",")
            f.write(el[self.qp_header[self.qp_log_header.TOOLN.value]])
            f.write(",")
            f.write(el[self.qp_header[self.qp_log_header.STIME.value]])
            f.write(",")
            f.write(el[self.qp_header[self.qp_log_header.STATS.value]])
            f.write("\n")
        f.close()
        return num_converged, avg_solve, perf_profile

    def qp_comp_export(self, log_list, problems_type):
        num_converged = 0
        avg_solve = 0
        perf_profile = []
        tool_name = log_list[0][self.qp_header[self.qp_log_header.TOOLN.value]]
        # file_name = os.path.join(".", "results", problems_type,tool_name,
        #                             "statistics.txt")
        file_name = "comp_" + problems_type +"_" + tool_name + "_statistics.txt"
        f = open(file_name, "w")
        f.write(self.qp_header[self.qp_log_header.PROBN.value])
        f.write(",")
        f.write(self.qp_header[self.qp_log_header.TOOLN.value])
        f.write(",Complementarity Inf")
        #f.write(self.qp_header[self.qp_log_header.STIME.value])
        #f.write(",status")
        #f.write(self.qp_header[self.qp_log_header.STATS.value])
        f.write("\n")
        for el in log_list:
            f.write(el[self.qp_header[self.qp_log_header.PROBN.value]])
            f.write(",")
            f.write(el[self.qp_header[self.qp_log_header.TOOLN.value]])
            f.write(",")
            f.write(el[self.qp_header[self.qp_log_header.CMINF.value]])
            f.write("\n")
        f.close()
        return num_converged, avg_solve, perf_profile


    def list_reduction(self, list_log, idx_eq, idx_sort, f_sort):
        '''
        :param list_log: The list of reports as a set of dictionaries
        :param idx_eq:
        :param idx_sort: The ind
        :param f_sort:
        :return: Remove duplicates that has the same idx_eq and take the one
        with best f_sort<idx_sort>()
        '''
        reduced_list = []
        if len(list_log) <= 1:
            return list_log
        tmp_list = []
        tmp_sort = np.array([])
        i = 0
        while i < len(list_log):
            num_rep = 0
            tmp_list.append(list_log[i])
            tmp_sort = np.concatenate((tmp_sort, [float(list_log[i][self.sympiler_header[idx_sort]])]))
            cur_mat_name = list_log[i][self.sympiler_header[idx_eq]]
            for j in range(i+1, len(list_log)):
                if cur_mat_name == list_log[j][self.sympiler_header[idx_eq]]:
                    # tmp_sort = np.concatenate((tmp_sort, list_log[j][self.sympiler_header[idx_sort]]))
                    tmp_sort = np.concatenate((tmp_sort, [float(list_log[j][self.sympiler_header[idx_sort]])]))
                    tmp_list.append(list_log[j])
                    num_rep += 1
                else:
                    break
            idx_fav = f_sort(tmp_sort)
            reduced_list.append(list_log[i + idx_fav])
            tmp_list.clear()
            tmp_sort = np.array([])
            i = i + num_rep + 1
        return reduced_list

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
            tool_name.append(lst[0][self.qp_header[self.qp_log_header.TOOLN.value]])
        num_tool = len(list_args)
        all_data = np.zeros((num_tool,n_groups))
        all_data_rel = np.zeros((num_tool, n_groups))
        mat_list = []
        lst_no = 0
        for lst in list_args:
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
            tool_name.append(lst[0][self.qp_header[self.qp_log_header.TOOLN.value]])
        num_tool = len(list_args)
        all_data = np.zeros((num_tool,n_groups))
        all_data_rel = np.zeros((num_tool, n_groups))
        mat_list = []
        lst_no = 0
        for lst in list_args:
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
        # for dim in range(num_tool):
        #     all_data_rel[dim, :] = all_data[dim, :]/all_data[0, :]
        # for v in range(n_groups):
        #     print(mat_list[v],',',all_data[0,v],',',all_data[1,v])

        ## Preparing the graph info
        # Other colors: https://matplotlib.org/gallery/color/named_colors.html
        colors = ('red', 'skyblue', 'indigo', 'olive', 'slateblue', 'magenta',
                  'slategray', 'limegreen', 'maroon', 'teal', 'khaki', 'purple',
                  'r', 'c', 'm', 'y', 'k')
        font_axis = self.font0.copy()
        fig, ax = plt.subplots(figsize=(12, 6))  # size=fig(10,20)
        width = 0.8
        ind = (num_tool + 1) * width * np.arange(n_groups)
        if y_label == '' and num_stacks == 1:
            y_label = self.sympiler_header[idx_y[0]]
        if file_name == '':
            file_name = self.sympiler_header[idx_y[0]]
        i = 0
        c = 0
        st = 0
        for linePlt in all_data_rel:
            st = 0
            bot_stacks = np.zeros(n_groups)
            for st_val in range(num_stacks):
                metric_label = self.sympiler_header[idx_y[st]]
                st_val_array = linePlt[:, st_val]
                ax.bar(ind + i*width,  st_val_array, width,
                       label=tool_name[i]+' '+metric_label, color=colors[c],
                       align='edge', bottom=bot_stacks)
                bot_stacks += st_val_array
                c += 1
                st += 1
            i += 1
        ax.legend(frameon=False)
        # ax.grid()
        font_axis.set_weight('black') # ['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']
        ax.set_ylabel(y_label,  fontsize=12, fontproperties=font_axis)
        ax.set_xticks(ind-6*width)
        ax.set_xticklabels(mat_list, fontsize=12)

        ax.tick_params(color='black', labelcolor='black')
        plt.yscale(y_scale)
        plt.xticks(rotation=55)
        # plt.set_title('Quality Progress for file:' + fileName)
        ax.tick_params(top=False, bottom=False, left=True, right=False,
                       labelleft=True, labelbottom=True, labeltop=False, labelright=False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.show()
        if self.path != '' and os.path.exists(self.path):
            fig.savefig(join(self.path, file_name+".pdf"))
        else:
            print("Output path is not valid, saving current!")
            fig.savefig(file_name+".pdf", bbox_inches='tight')

        return all_data_rel

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


def main(argv):
    output_path = ''
    input_path1 = ''
    input_path2 = ''
    input_path3 = ''
    input_path4 = ''
    report_type = 5
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
    print('Reading CSV files from the input directory: ', input_path1)
    #report_type = 2
    if report_type == reporting_type.LinSolve2.value:
        fact_or_acc = 0  # 1: acc   0: fact
        key_item = graphGenerator.sympiler_header[Visualizer.LogHeader.MATN.value]
        if (not os.path.exists(input_path1) or not os.path.exists(input_path2)):
            print('Please specify a path and right command arguments. '
                  'You need to have correct input paths to see a graph')
            print('visualize.py -i <sympiler report path> -f <mkl report path> -o <output Path for graphs> ')
            sys.exit(-1)

        del_items = ["boyd1", "boyd2", "cont-201", "cont-300", "dynamicSoaringProblem_5",
                     "dynamicSoaringProblem_6", "dynamicSoaringProblem_7", "dynamicSoaringProblem_8",
                     "lowThrust_1", "lowThrust_2", "lowThrust_3", "lowThrust_4", "lowThrust_5", "lowThrust_6",
                     "lowThrust_7", "lowThrust_8", "lowThrust_9", "lowThrust_10", "lowThrust_11","lowThrust_12",
                     "lowThrust_13", "reorientation_5", "reorientation_6", "reorientation_7", "reorientation_8",
                     "freeFlyingRobot_5", "freeFlyingRobot_6", "freeFlyingRobot_7", "freeFlyingRobot_8",
                     "freeFlyingRobot_9", "freeFlyingRobot_10", "freeFlyingRobot_11", "freeFlyingRobot_12",
                     "freeFlyingRobot_13", "freeFlyingRobot_14", "freeFlyingRobot_15", "freeFlyingRobot_16",
                     "spaceStation_5", "spaceStation_6", "spaceStation_7", "spaceStation_8", "spaceStation_9",
                     "spaceStation_10","spaceStation_11","spaceStation_12","spaceStation_13","spaceStation_14"]
        sym_log4 = graphGenerator.read_sympiler_log(input_path1, Visualizer.LogHeader.MATN.value,
                                                   Visualizer.LogHeader.FACT.value, min_idx)
        if input_path3 != '':
            sym_log3 = graphGenerator.read_sympiler_log(input_path3, Visualizer.LogHeader.MATN.value,
                                                   Visualizer.LogHeader.FACT.value, min_idx)
        #mkl_log = graphGenerator.fact_time_mkl(input_path2)
        mkl_log = graphGenerator.read_sympiler_log(input_path2, Visualizer.LogHeader.MATN.value,
                                                   Visualizer.LogHeader.FACT.value, min_idx)

        for del_item in del_items:
            sym_log4 = exclude_items(sym_log4, del_item, key_item)
            mkl_log = exclude_items(mkl_log, del_item, key_item)

        # Writing the CSV file if necessary
        graphGenerator.write_to_csv(mkl_log,'out_csv_mkl.csv')
        graphGenerator.write_to_csv(sym_log4, 'out_csv_sym.csv')
        print('Drawing a graph')
        # idx_y_val = np.array([
        #     Visualizer.LogHeader.ANAT.value,
        #     Visualizer.LogHeader.FACT.value,
        #     Visualizer.LogHeader.SOLT.value,
        #     Visualizer.LogHeader.ITERT.value
        # ])
        f_val = nop
        y_scale = "linear"
        file_name = "factorization"
        y_label = "MKL Pardiso/ LBL Fact Time"
        if fact_or_acc == 1:
            idx_y_val = np.array([Visualizer.LogHeader.REDA.value])
            y_scale = "log"
            file_name = "accuracy"
            y_label = "Residual Inf Norm"
        else:
            idx_y_val = np.array([Visualizer.LogHeader.FACT.value])
            f_val = scaling_array
        t_sym = []
        t_mkl = []
        for it in sym_log4:
            t_sym.append(it[graphGenerator.sympiler_header[idx_y_val[0]]])
        for it in mkl_log:
            t_mkl.append(it[graphGenerator.sympiler_header[idx_y_val[0]]])

        graphGenerator.visualize_list(sym_log4, mkl_log,
                                      idx_x=Visualizer.LogHeader.MATN.value,
                                      idx_y=idx_y_val,
                                      file_name=file_name,
                                      y_label=y_label,
                                      f_val=f_val,
                                      y_scale=y_scale
                                      )
    elif report_type == reporting_type.AllFile.value:
        if not os.path.isdir(input_dir):
            print('Please specify a path to the input directory.')
            print('visualize.py -i <sympiler report path> -f <mkl report path> -o <output Path for graphs> ')
            sys.exit(-1)

        def condition(dic_in):
            p_size = graphGenerator.qp_header[Visualizer.qp_log_header.PROBS.value]
            #if dic_in[p_size] <= 5000:
            #if dic_in[p_size] > 5000 and dic_in[p_size] <= 200000:
            #if dic_in[p_size] > 200000:

            if dic_in[p_size] <= 20000:
                return True
            else:
                return False

        def condition_ql(dic_in):
            p_size = graphGenerator.qp_header[Visualizer.qp_log_header.PROBD.value]
            psize = int(dic_in[p_size])
            if psize < 8000:
                # if dic_in[p_size] > 5000 and dic_in[p_size] <= 200000:
                # if dic_in[p_size] > 200000:
                return True
            else:
                return False

        def condition_ql(dic_in):
            p_size = graphGenerator.qp_header[Visualizer.qp_log_header.PROBD.value]
            psize = int(dic_in[p_size])
            if psize < 8000:
                # if dic_in[p_size] > 5000 and dic_in[p_size] <= 200000:
                # if dic_in[p_size] > 200000:
                return True
            else:
                return False

        def condition_large_var(dic_in):
            p_dim = graphGenerator.qp_header[Visualizer.qp_log_header.PROBD.value]
            p_ed = graphGenerator.qp_header[Visualizer.qp_log_header.CNSED.value]
            p_id = graphGenerator.qp_header[Visualizer.qp_log_header.CNSID.value]
            dim_val = float(dic_in[p_dim])
            ed_val = int(dic_in[p_ed])
            id_val = int(dic_in[p_id])

            if (ed_val+id_val / dim_val) >= 1.9:
                # if dic_in[p_size] > 5000 and dic_in[p_size] <= 200000:
                # if dic_in[p_size] > 200000:
                return True
            else:
                return False

        problem_group = ['All']
        #problem_group = ['Contact simulation']
        #problem_group = ['Maros Mezaros']
        #problem_group = ['MPC']
        #problem_group = ['Cloth simulation', 'Image deformation', 'Model construction', 'Mesh processing']
        if len(problem_group) > 1:
            problem_type = "Mixed"
        else:
            problem_type = problem_group[0]
        tol = 1e-03
        tol_str = '0.001'
        sorting_idx = Visualizer.qp_log_header.PROBN
        sorting_idx2 = Visualizer.qp_log_header.PROBS
        key_item = graphGenerator.qp_header[Visualizer.qp_log_header.PROBN.value]
        ql_wanted = 0
        del_items = ["bunny2_bunny2","wolf_falling6","knight_falling2",
                     "InjBnd_Mapping0","InjBnd_Mapping1",
                     "InjBnd_Mapping2","InjBnd_Mapping3",
                     #"InjBnd_MappingImp0","InjBnd_MappingImp1","InjBnd_MappingImp2",
                     "lamb_falling1","horse_horse5","horse_horse1",
                     "horse_falling1","garg3_floor1","bunny_falling1","bunny2_falling1","bunny_bunny2",
                     "brick0_brick","screwdriver_falling"]
        for i in range(1, 130):
            item = "Beam_floor"+str(i)
            del_items.append(item)
        for i in range(290, 338):
            item = "Beam_floor"+str(i)
            del_items.append(item)
        for i in range(400, 524):
            item = "Cube_floor"+str(i)
            #del_items.append(item)
        for i in range(4, 12):
            item = "lion"+str(i)+"_s"
            del_items.append(item)
        for i in range(5, 800):
            item = "Bar2k_floor"+str(i)
            del_items.append(item)
        for i in range(10, 19):
            item = "nail"+str(i)+"t"
            del_items.append(item)

        files = [f for f in glob.glob(input_dir + "*.csv", recursive=False)]
        all_logs = []
        tool_names = []
        for f in files:
            if 'ql' in f or 'quadprog' in f:
                continue
            graphGenerator.find_problem_size(f, sorting_idx,
                                                    Visualizer.qp_log_header.STIME,
                                                        tol, tol_str, problem_group)
        for f in files:
            if 'ql' in f or 'quadprog' in f:
                continue
            logs_tmp = graphGenerator.read_qp_log_split(f, sorting_idx,
                                                    Visualizer.qp_log_header.STIME,
                                                        tol, tol_str, problem_group)

            for log_tmp in logs_tmp:
                tn = log_tmp[0][graphGenerator.qp_header[Visualizer.qp_log_header.TOOLN.value]]
                # if not "1e-09_1_1_1e-15" in tn and not "1e-09_2_2_1e-15" in tn and \
                #         not "1e-09_3_3_1e-15" in tn and not "OSQP" in tn:
                #     continue
                for del_item in del_items:
                    log_tmp = exclude_items(log_tmp, del_item, key_item)
                #if ql_wanted:

                #log_tmp = exclude_items_general(log_tmp, condition)
                #log_tmp = exclude_items_general(log_tmp, condition_ql)
                #log_tmp = exclude_items_general(log_tmp, condition_large_var)
                #log_tmp = sorted(log_tmp, key=itemgetter(graphGenerator.qp_header[sorting_idx2.value]))
                log_tmp = sorted(log_tmp, key=itemgetter(graphGenerator.qp_header[sorting_idx.value]))
                all_logs.append(log_tmp)
                tool_names.append(tn)
                graphGenerator.qp_log_export(log_tmp, problem_type)
                graphGenerator.qp_comp_export(log_tmp, problem_type)
        compute_performance_profiles(tool_names, problem_type,tol_str)
        su = compute_speedup(tool_names,"NASOQ-Fixed",problem_type,tol_str)
        fr = compute_failure_rates(tool_names,
                              problem_type, tol_str)
        plot_performance_profiles(problem_type,
                                  tool_names, tol_str)
        min_fr = 100
        min_log = []
        all_reduced_logs = []
        for l in all_logs:
            ln = l[0][graphGenerator.qp_header[Visualizer.qp_log_header.TOOLN.value]]
            if "NASOQ" in ln:
                cur_fr = fr[ln]
                if cur_fr < min_fr:
                    min_fr = cur_fr
                    min_log = l
            else:
                all_reduced_logs.append(l)

        all_reduced_logs.append(min_log)

        graphGenerator.plot_bar_chart(all_reduced_logs, idx_x=Visualizer.qp_log_header.PROBN,
                                   idx_y=Visualizer.qp_log_header.STIME,
                                   idx_y2=Visualizer.qp_log_header.STATS,
                                   y_label='Solve Time with tol = %s (Sec) ' % tol_str,
                                      file_name='results_time_%s_%s.png' % (problem_type, tol_str) )
       # graphGenerator.write_to_csv(nasoq_log1, 'NASOQ.csv')


    else:
        if (not os.path.exists(input_path1) or not os.path.exists(input_path2)):
            print('Please specify a path and right command arguments. '
                  'You need to have correct input paths to see a graph')
            print('visualize.py -i <sympiler report path> -f <mkl report path> -o <output Path for graphs> ')
            sys.exit(-1)
        nasoqp_log2 = []
        nasoqp_log3 = []
        nasoqp_log4 = []
        problem_type = 'AllQPs'
        tol = 1e-9
        sorting_idx = Visualizer.qp_log_header.PROBN
        sorting_idx2 = Visualizer.qp_log_header.PROBS
        nasoq_log1 = graphGenerator.read_qp_log(input_path1, sorting_idx,
                                                    Visualizer.qp_log_header.STIME, tol)
        graphGenerator.write_to_csv(nasoq_log1, 'NASOQall.csv')

        nasoq_log1 = sorted(nasoq_log1, key=itemgetter(graphGenerator.qp_header[sorting_idx2.value]))
        graphGenerator.qp_log_export(nasoq_log1, problem_type)
        if input_path2 != '':
            nasoqp_log2 = graphGenerator.read_qp_log(input_path2, sorting_idx,
                                                    Visualizer.qp_log_header.STIME,tol)
            nasoqp_log2 = sorted(nasoqp_log2, key=itemgetter(graphGenerator.qp_header[sorting_idx2.value]))
            graphGenerator.qp_log_export(nasoqp_log2, problem_type)
        if input_path3 != '':
            nasoqp_log3 = graphGenerator.read_qp_log(input_path3, sorting_idx,
                                                    Visualizer.qp_log_header.STIME,tol)
            nasoqp_log3 = sorted(nasoqp_log3, key=itemgetter(graphGenerator.qp_header[sorting_idx2.value]))
            graphGenerator.qp_log_export(nasoqp_log3, problem_type)
        if input_path4 != '':
            nasoqp_log4 = graphGenerator.read_qp_log(input_path4, sorting_idx,
                                                    Visualizer.qp_log_header.STIME,tol)
            nasoqp_log4 = sorted(nasoqp_log4, key=itemgetter(graphGenerator.qp_header[sorting_idx2.value]))
            graphGenerator.qp_log_export(nasoqp_log4, problem_type)
        # Writing the CSV file if necessary
        #graphGenerator.write_to_csv(osqp_log, 'out_csv_mkl.csv')
        #graphGenerator.write_to_csv(nasoq_log1, 'out_csv_sym.csv')
        print('Drawing a graph')

        if report_type == reporting_type.SCATQP.value:
            graphGenerator.compute_scattered(nasoq_log1,nasoqp_log2,idx_x=Visualizer.qp_log_header.PROBN,
                                      idx_y=Visualizer.qp_log_header.STIME,
                                         idx_y2 = Visualizer.qp_log_header.STATS,
                                             y_label='Solve Time(sec)')
        elif report_type == reporting_type.BARQP.value:
            graphGenerator.compute_bar(nasoq_log1, nasoqp_log2, idx_x=Visualizer.qp_log_header.PROBN,
                                             idx_y=Visualizer.qp_log_header.STIME,
                                             idx_y2=Visualizer.qp_log_header.STATS,
                                                y_label='Solve Time(sec)')
        elif report_type == reporting_type.PerfProf.value:
            compute_performance_profiles(["NASOQ","OSQP","OSQP-polished"],problem_type)
            compute_failure_rates(["NASOQ", "OSQP","OSQP-polished"],
                                  problem_type)
            plot_performance_profiles(problem_type,
                                      ["NASOQ","OSQP","OSQP-polished"])
            graphGenerator.write_to_csv(nasoq_log1, 'NASOQ.csv')
            graphGenerator.write_to_csv(nasoqp_log2, 'OSQP.csv')
            graphGenerator.write_to_csv(nasoqp_log3, 'OSQP-polished.csv')
            #graphGenerator.write_to_csv(nasoqp_log3, 'NASOQ-updown.csv')
            #graphGenerator.write_to_csv(nasoqp_log4, 'OSQP-polished.csv')
            # idx_y_val = np.array([
            #     Visualizer.LogHeader.ANAT.value,
            #     Visualizer.LogHeader.FACT.value,
            #     Visualizer.LogHeader.SOLT.value,
            #     Visualizer.LogHeader.ITERT.value
            # ])
        else:
            print(" No report type is selected.")
            # idx_y_val = np.array([Visualizer.LogHeader.NRM1.value])
            # graphGenerator.visualize_list(nasoq_log1, osqp_log,
            #                               idx_x=Visualizer.LogHeader.MATN.value,
            #                               idx_y=idx_y_val,
            #                               # y_label="MKL/ParSy Fact time",
            #                               file_name="accuracy1",
            #                               # y_label="MKL Pardiso/ ParSy LDL Time",
            #                               f_val=nop,
            #                               y_scale="log"
            #                               )
    sys.exit(0)

if __name__ == "__main__":
   main(sys.argv[1:])