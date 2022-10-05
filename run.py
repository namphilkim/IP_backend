import IP_functionality as iplib
import time
import sys

def main(workdir_list):
  start_time = time.time()
  n = len(workdir_list)
  isotype_dict = {'IGHM':0,'IGHD':1,'IGHG3':2,'IGHG1':3,'IGHA1':4,'IGHG2':5,'IGHG4':6,'IGHA2':7}

  isotpye_list = []
  gini_list = []
  shm_list = []
  VJ_dict_list = []
  V_dict_list = []
  J_dict_list = []
  
  for i in range(n):
    isotpye_list.append([])
    gini_list.append([])
    shm_list.append([])
    VJ_dict_list.append([])
    V_dict_list.append([])
    J_dict_list.append([])
    isotype_list[-1], gini_list[-1], shm_list[-1], VJ_dict[-1], V_dict[-1], J_dict[-1] = iplib.get_repertoire(workdir_list[i],isotype_dict)

  print("--- %s seconds ---" % (time.time() - start_time))
  
  iplib.basic_box_plot(data_list=isotype_list,
               x_labels=['M','D','G3','G1','A1','G2','G4','A2'],colors = ['b','r','g','orange','purple','brown'],
               size_x=12,size_y=8,y_lim=1.0,y_label='Frequency',has_jitter=True,save='Isotype.png')

  iplib.basic_box_plot(data_list=gini_list,
               x_labels=['M','D','G3','G1','A1','G2','G4','A2','All'],colors = ['b','r','g','orange','purple','brown'],
               size_x=12,size_y=8,y_lim=1.0,y_label='Clonality',has_jitter=True,save='Clonality.png')

  iplib.basic_box_plot(data_list=shm_list,
               x_labels=['M','D','G3','G1','A1','G2','G4','A2','All'],colors = ['b','r','g','orange','purple','brown'],
               size_x=12,size_y=8,y_lim=50.0,y_label='Somatic mutation',has_jitter=True,save='Mutation.png')

main(list(sys.argv))
