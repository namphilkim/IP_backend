import matplotlib.pyplot as plt
import os
from scipy.stats import ttest_ind
import random as rd

def get_repertoire(workdir, isotype_dict):
                   
    file_list = os.listdir(workdir)
    
    isotype_list = []
    gini_list = [[]]
    shm_list = [[]]
    all_freq_list = []
    VJ_dict = {}
    V_dict = {}
    J_dict = {}
    for i in range(len(isotype_dict)):
        isotype_list.append([])
        gini_list.append([])
        shm_list.append([])

    scatter_shm = []
    scatter_freq = []
    for f,file in enumerate(file_list):
        total_reads = 0
        file_pass = open(workdir+file,'r')
        for l, line in enumerate(file_pass):
            split_line_pass = line.split('\t')
            list_line_pass = split_line_pass[:-1]+[split_line_pass[-1][:-1]]
            if l == 0:
                read_index = list_line_pass.index('duplicate_count')
            else:
                total_reads += int(list_line_pass[read_index])
                
        print(file, 'Total reads : ' + str(total_reads))
        file_pass = open(workdir+file,'r')
        isotype = []
        shm = [0]
        gini = [0]
        freq_list = [[]]
        unique_seqs = [0]
        for i in range(len(isotype_dict)):
            isotype.append(0)
            shm.append(0)
            gini.append(0)
            freq_list.append([])
            unique_seqs.append(0)
        for l, line in enumerate(file_pass):
            split_line_pass = line.split('\t')
            list_line_pass = split_line_pass[:-1]+[split_line_pass[-1][:-1]]
            if l == 0:
                read_index = list_line_pass.index('duplicate_count')
                isotype_index = list_line_pass.index('c_call')
                shm_index = list_line_pass.index('v_alignment_mutation')
                v_index = list_line_pass.index('v_call')
                j_index = list_line_pass.index('j_call')
            else:
                unique_seqs[-1] += 1
                this_frequency = int(list_line_pass[read_index])/total_reads
                freq_list[-1].append(this_frequency)
                shm[-1] += int(list_line_pass[shm_index])*this_frequency
                try:
                    this_isotype = isotype_dict[list_line_pass[isotype_index].split('*')[0]]
                    isotype[this_isotype] += this_frequency
                    shm[this_isotype] += int(list_line_pass[shm_index])*this_frequency
                    freq_list[this_isotype].append(this_frequency)
                    unique_seqs[this_isotype] += 1
                except KeyError:
                    pass
                if VJ_dict.get(list_line_pass[v_index].split('*')[0]+'_'+list_line_pass[j_index].split('*')[0]) == None:
                    if f == 0:
                        VJ_dict[list_line_pass[v_index].split('*')[0]+'_'+list_line_pass[j_index].split('*')[0]] = [this_frequency]
                    else:
                        temp = []
                        for i in range(f):
                            temp.append(0)
                        VJ_dict[list_line_pass[v_index].split('*')[0]+'_'+list_line_pass[j_index].split('*')[0]] = temp + [this_frequency]
                else: 
                    try:
                        VJ_dict[list_line_pass[v_index].split('*')[0]+'_'+list_line_pass[j_index].split('*')[0]][f] += this_frequency
                    except IndexError:
                        VJ_dict[list_line_pass[v_index].split('*')[0]+'_'+list_line_pass[j_index].split('*')[0]].append(this_frequency) 
                   
                if V_dict.get(list_line_pass[v_index].split('*')[0]) == None:
                    if f == 0:
                        V_dict[list_line_pass[v_index].split('*')[0]] = [this_frequency]
                    else:
                        temp = []
                        for i in range(f):
                            temp.append(0)
                        V_dict[list_line_pass[v_index].split('*')[0]] = temp + [this_frequency]
                else: 
                    try:
                        V_dict[list_line_pass[v_index].split('*')[0]][f] += this_frequency
                    except IndexError:
                        V_dict[list_line_pass[v_index].split('*')[0]].append(this_frequency) 
                   
                if J_dict.get(list_line_pass[j_index].split('*')[0]) == None:
                    if f == 0:
                        J_dict[list_line_pass[j_index].split('*')[0]] = [this_frequency]
                    else:
                        temp = []
                        for i in range(f):
                            temp.append(0)
                        J_dict[list_line_pass[j_index].split('*')[0]] = temp + [this_frequency]
                else: 
                    try:
                        J_dict[list_line_pass[j_index].split('*')[0]][f] += this_frequency
                    except IndexError:
                        J_dict[list_line_pass[j_index].split('*')[0]].append(this_frequency) 
    
        for i in range(len(isotype_dict)+1):
            length = len(freq_list[i])
            temp = sorted(freq_list[i])
            all_freq_list.append(temp)
            total = sum(freq_list[i])
            if total != 0:
                gini_index = 0.0
                cumulative = 0.0
                for f,freq in enumerate(temp):
                    a = total*f/length-cumulative
                    cumulative += freq 
                    b = total*(f+1)/length-cumulative
                    gini_index += (a+b)/length
                gini[i] = gini_index/total     
            elif total == 0:
                gini[i] = 0

        for i in range(len(isotype_dict)):
            isotype_list[i].append(isotype[i])
            gini_list[i].append(gini[i])
            if isotype[i] != 0:
                shm_list[i].append(shm[i]/isotype[i])
            elif isotype[i] == 0:
                shm_list[i].append(0)

        gini_list[-1].append(gini[-1])
        shm_list[-1].append(shm[-1])
    return isotype_list, gini_list, shm_list, VJ_dict, V_dict, J_dict
             
def plot_significance(list_1,list_2,x1,x2,y,y_offset,fontsize,linewidth):
    t,p = ttest_ind(list_1,list_2)
    print(p)
    if p <= 0.05:
        plt.plot([x1,x2], [y, y], linewidth=linewidth, color='k')
    if p <= 0.05 and p > 0.01:
        plt.text((x1+x2)/2,y+y_offset,'*',fontsize=fontsize,horizontalalignment='center')
    elif p <= 0.01 and p > 0.001:
        plt.text((x1+x2)/2,y+y_offset,'**',fontsize=fontsize,horizontalalignment='center')
    elif p <= 0.001 and p > 0.0001:
        plt.text((x1+x2)/2,y+y_offset,'***',fontsize=fontsize,horizontalalignment='center')
    elif p <= 0.0001:
        plt.text((x1+x2)/2,y+y_offset,'****',fontsize=fontsize,horizontalalignment='center')
    return
                   
def basic_box_plot(data_list,x_labels,colors,size_x,size_y,y_lim,y_label='Frequency',has_jitter=True):
    num_cohort = len(data_list)
    num_x_val = len(data_list[0])
    bias = 1/(num_cohort+1)
    line_width = 2
    plt.rcParams["figure.figsize"] = (size_x,size_y) 
    fig,ax = plt.subplots(1,1)
    positions = []
    widths = []
    for i in range(num_cohort):
        positions.append([])
        widths.append([])
        for j in range(num_x_val):
            positions[-1].append(j+(i+1)*bias)
            widths[-1].append(bias*0.8)
    
    for i in range(num_cohort):
        plt.boxplot(data_list[i],sym='',positions = positions[i],widths = widths[i],patch_artist=True,
                         boxprops = {'facecolor':colors[i],'edgecolor':'#555555','linewidth':line_width},
                        whiskerprops = {'color':'#555555','linewidth':line_width},
                        capprops = {'color':'#555555','linewidth':line_width},
                        medianprops = {'color':'#555555','linewidth':line_width})

    ticks = x_labels
    x_ticks = []
    for j in range(num_x_val):
        x_ticks.append(0.5+j)
    plt.xticks(x_ticks, ticks,font='serif',fontsize = 20,fontweight='bold')
    for j in range(num_x_val):
        for i in range(1,num_cohort):
            plot_significance(data_list[0][j],data_list[i][j],j+bias,j+(i+1)*bias,max(max(data_list[0][j]),max(data_list[i][j]))*1.1,0.005,24,4)
        for i in range(num_cohort):
            for k in range(len(data_list[i][j])):
                random_x = rd.random()-0.5
                plt.scatter(j+(i+1)*bias+0.5*bias*random_x,data_list[i][j][k],s=8,color = 'k',zorder=10)

    plt.yticks(font='serif',fontsize = 20,fontweight='bold')
    plt.xlim(0, num_x_val)
    plt.ylim(0, y_lim)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.ylabel(y_label,font='serif',fontsize = 24,fontweight='bold')
    plt.xlabel('Isotype',font='serif',fontsize = 24,fontweight='bold')
    #plt.title('Mutation count',font='serif',fontsize = 48,fontweight='bold')
    plt.tight_layout()
    #plt.savefig(work_dir+'figures/in-house_DB_divergence_box_and_violin.png')
    plt.show()