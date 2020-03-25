
import numpy as np
import os
import pandas as pd
import pandas as pd
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pylab as plt
import csv

# Solver constants
OPTIMAL = "optimal"
OPTIMAL_INACCURATE = "optimal inaccurate"
PRIMAL_INFEASIBLE = "primal infeasible"
PRIMAL_INFEASIBLE_INACCURATE = "primal infeasible inaccurate"
DUAL_INFEASIBLE = "dual infeasible"
DUAL_INFEASIBLE_INACCURATE = "dual infeasible inaccurate"
PRIMAL_OR_DUAL_INFEASIBLE = "primal or dual infeasible"
SOLVER_ERROR = "solver_error"
MAX_ITER_REACHED = "max_iter_reached"
TIME_LIMIT = "time_limit"

SOLUTION_PRESENT = [OPTIMAL]
MAX_TIMING = 1e8

MAX_VAL = 1600


def compute_performance_profiles(solvers, problems_type, tol=''):
    t = {}
    status = {}

    # Get time and status
    for solver in solvers:
        # path = os.path.join('.', 'results', problems_type,
        #                     solver, 'results.csv')
        file_name = "results_" + problems_type + "_" + solver + "_statistics.txt"
        path = os.path.join('.',file_name)
        df = pd.read_csv(path)

        # Get total number of problems
        n_problems = len(df)

        t[solver] = df['run_time'].values
        status[solver] = df['status'].values

        # Set maximum time for solvers that did not succeed
        for idx in range(n_problems):
            if status[solver][idx] not in SOLUTION_PRESENT:
                t[solver][idx] = MAX_TIMING

    r = {}  # Dictionary of relative times for each solver/problem
    for s in solvers:
        r[s] = np.zeros(n_problems)

    # Iterate over all problems to find best timing between solvers
    for p in range(n_problems):
        # Get minimum time
        min_time = np.min([t[s][p] for s in solvers])

        # Normalize t for minimum time
        if min_time == 0:
            print("division by zero")
        for s in solvers:
            r[s][p] = t[s][p]/min_time

    # Compute curve for all solvers
    n_tau = 1000
    tau_vec = np.logspace(0, 4, n_tau)
    rho = {'tau': tau_vec}  # Dictionary of all the curves

    for s in solvers:
        rho[s] = np.zeros(n_tau)
        for tau_idx in range(n_tau):
            count_problems = 0  # Count number of problems with t[p, s] <= tau
            for p in range(n_problems):
                if r[s][p] <= tau_vec[tau_idx]:
                    count_problems += 1
            rho[s][tau_idx] = count_problems / n_problems

    # Store final pandas dataframe
    df_performance_profiles = pd.DataFrame(rho)
    performance_profiles_file = os.path.join('.', 'results'+'_'+
                                             problems_type+'_'+
                                             'performance_profiles.csv')
    df_performance_profiles.to_csv(performance_profiles_file, index=False)

    # Plot performance profiles
    # import matplotlib.pylab as plt
    # for s in solvers:
    #     plt.plot(tau_vec, rho[s], label=s)
    # plt.legend(loc='best')
    # plt.ylabel(r'$\rho_{s}$')
    # plt.xlabel(r'$\tau$')
    # plt.grid()
    # plt.xscale('log')
    # plt.show(block=False)


def compute_failure_rates(solvers, problems_type, tol=''):
    """
    Compute and show failure rates
    """

    # Check if results file already exists
    results_file = os.path.join(".", "failure_" + problems_type + "_"+tol+
                                "_results.txt")

    f_rate = {}
    # Always overwrite file
    f = open(results_file, "w")

    f.write('[Failure rates]\n')
    for solver in solvers:
        # results_file = os.path.join('.', 'results', problems_type,
        #                             solver, 'results.csv')
        file_name = "results_" + problems_type + "_" + solver + "_statistics.txt"
        results_file = os.path.join('.', file_name)
        df = pd.read_csv(results_file)

        n_problems = len(df)

        failed_statuses = np.logical_and([df['status'].values != s
                                           for s in SOLUTION_PRESENT],True)
        n_failed_problems = np.sum(failed_statuses)
        failure_rate = n_failed_problems / n_problems
        f_rate[solver] = failure_rate
        f.write(" - %s = %.4f %%\n" % (solver, 100 * failure_rate))
    f.write("\n")

    f.close()
    return f_rate

def plot_performance_profiles(problems, solvers, tol=''):
    """
    Plot performance profiles in matplotlib for specified problems and solvers
    """
    df = pd.read_csv('./results_%s_performance_profiles.csv' % problems)
    plt.figure(0)
    for solver in solvers:
        plt.plot(df["tau"], df[solver], label=solver)
    plt.xlim(1., 10000.)
    plt.ylim(0., 1.)
    plt.xlabel(r'Performance ratio $\tau$')
    plt.ylabel('Ratio of problems solved (tol = %s)' %tol)
    plt.xscale('log')
    plt.legend()
    plt.grid()
    results_file = './results_perfprofile_%s_%s.png' % (problems, tol)
    print("Saving plots to %s" % results_file)
    plt.savefig(results_file, dpi=100)
    plt.show(block=False)


def exclude_items(log_lst, del_items, key_items):
    res = [i for i in log_lst if (not del_items  == i[key_items] )]
    return res


def compare_size(l_item, thresh, idx):
    if l_item(idx) > thresh:
        return True
    return False


def exclude_items_general(log_lst, condition):
    res = [i for i in log_lst if (condition(i))]
    return res


def compute_speedup(solvers, ref_solver, problems_type, tol=''):
    t = {}
    status = {}
    speedup = {}
    if ref_solver == "":
        ref_solver = solvers[0]
    # Get time and status
    for solver in solvers:
        # path = os.path.join('.', 'results', problems_type,
        #                     solver, 'results.csv')
        file_name = "results_" + problems_type + "_" + solver + "_statistics.txt"
        path = os.path.join('.', file_name)
        df = pd.read_csv(path)

        # Get total number of problems
        n_problems = len(df)

        t[solver] = df['run_time'].values
        status[solver] = df['status'].values

        # Set maximum time for solvers that did not succeed
        for idx in range(n_problems):
            if status[solver][idx] not in SOLUTION_PRESENT:
                t[solver][idx] = MAX_TIMING

    r = {}  # Dictionary of relative times for each solver/problem
    for s in solvers:
        r[s] = np.zeros(n_problems)
        speedup[s] = np.zeros(1)

    # Iterate over all problems to find best timing between solvers
    for p in range(n_problems):
        # Get minimum time
        ref_time = t[ref_solver][p]

        # Normalize t for minimum time
        if ref_time == 0:
            print("division by zero")
        for s in solvers:
            if t[s][p] < MAX_VAL and t[ref_solver][p] < MAX_VAL:
                r[s][p] = ref_time / t[s][p]
                speedup[s] = speedup[s] + 1
            else:
                r[s][p] = 0 # make it zero so it does not affect avg


    # Store final pandas dataframe
    speedup_new = {}
    speedup_file = os.path.join('.', 'speedup' + '_' +
                                problems_type + '_' + tol + '_'
                                                            'results.csv')

    header = ['Tool name', 'Average speedup', 'Max', 'Min', 'SDE', 'Median', '# of slower', '# of 2X slower',
              'Tol', 'Problem type', '# of problems']
    tmp_dic = {}
    with open(speedup_file, 'w') as csv_file:
        writer = csv.DictWriter(csv_file,fieldnames=header)
        writer.writeheader()
        for s in solvers:
            tmp_dic[header[0]] = s
            sum_time = np.sum(r[s])
            max_time = np.max(r[s])
            median_array = []
            slower_cases_no = 0
            slower2_cases_no = 0
            median_time = 0
            min_time = 1e10
            for t in r[s]:
                if min_time > t and t > 0:
                    min_time = t

            for t in r[s]:
                if t>0:
                    median_array.append(t)
                    if t < 0.9:
                        slower_cases_no = slower_cases_no + 1
                    if t < 0.5:
                        slower2_cases_no = slower2_cases_no + 1
            median_time = np.median(median_array)
            avg = sum_time/speedup[s][0]
            tmp_dic[header[1]] = avg
            var = 0.0
            for v in r[s]:
                if v > 0:
                    var += np.sum(pow((v - avg), 2))
            variance = var/speedup[s][0]

            tmp_dic[header[2]] = max_time
            tmp_dic[header[3]] = min_time
            tmp_dic[header[4]] = np.sqrt(variance)
            tmp_dic[header[5]] = median_time
            tmp_dic[header[6]] = slower_cases_no
            tmp_dic[header[7]] = slower2_cases_no
            tmp_dic[header[8]] = tol
            tmp_dic[header[9]] = problems_type
            tmp_dic[header[10]] = n_problems
            writer.writerow(tmp_dic)

def get_name(s):
    full_name = s
    mat_full_name = full_name.split('/')
    mat_name = mat_full_name[len(mat_full_name) - 1].split('.')[0]
    us_pos = mat_name.rfind('_')
    if us_pos > 0:
        mat_name = mat_name[0:us_pos]
    return mat_name

def compute_speedup_general(list, ref_idx):
    num_tools = len(list)
    num_entry = len(list[0])
    speed_up = np.zeros((num_tools, num_entry))
    for i in range(num_tools):
        for j in range(num_entry):
            speed_up[i, j] = list[ref_idx, j] / list[i, j]
    return speed_up


def plot_bar_chart(all_data, tool_name, y_label, file_name,font0):
    colors = ('red', 'skyblue', 'indigo', 'olive', 'slateblue', 'magenta',
              'slategray', 'limegreen', 'maroon', 'teal', 'khaki', 'purple',
              'r', 'c', 'm', 'y', 'k')
    markers = ('s', 'o', '+', ".", ",", "v", "^", "<", ">",
               "1", "2", "3", "4", "8", "s", "p", "P")
    num_tool = len(tool_name)
    mat_list = len(len(all_data[0,:]))
    width = 0.3
    ind = (num_tool + 1) * width * np.arange(len(mat_list))
    x_label = np.arange(1, len(mat_list) + 1)
    font_axis = font0.copy()
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