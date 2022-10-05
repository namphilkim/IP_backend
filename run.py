import IP_functionality as iplib
import time

def main():
  start_time = time.time()
  isotype_dict = {'IGHM':0,'IGHD':1,'IGHG3':2,'IGHG1':3,'IGHA1':4,'IGHG2':5,'IGHG4':6,'IGHA2':7}

  control_workdir = 'sjogren/in_house_healthy/'
  HC_isotype_list, HC_gini_list, HC_shm_list, HC_VJ_dict, HC_V_dict, HC_J_dict = iplib.get_repertoire(control_workdir,isotype_dict)

  SS_workdir = 'sjogren/sjogren_file/'
  SS_isotype_list, SS_gini_list, SS_shm_list, SS_VJ_dict, SS_V_dict, SS_J_dict = iplib.get_repertoire(SS_workdir,isotype_dict)

  AD_workdir = 'AD/AD_Second/'
  AD_isotype_list, AD_gini_list, AD_shm_list, AD_VJ_dict, AD_V_dict, AD_J_dict = iplib.get_repertoire(AD_workdir,isotype_dict)
  print("--- %s seconds ---" % (time.time() - start_time))
  
  iplib.basic_box_plot(data_list=[HC_isotype_list,SS_isotype_list,AD_isotype_list],
               x_labels=['M','D','G3','G1','A1','G2','G4','A2'],colors = ['b','r','g'],
               size_x=12,size_y=8,y_lim=1.0,y_label='Frequency',has_jitter=True)

iplib.basic_box_plot(data_list=[HC_gini_list,SS_gini_list,AD_gini_list],
               x_labels=['M','D','G3','G1','A1','G2','G4','A2','All'],colors = ['b','r','g'],
               size_x=12,size_y=8,y_lim=1.0,y_label='Clonality',has_jitter=True)

iplib.basic_box_plot(data_list=[HC_shm_list,SS_shm_list,AD_shm_list],
               x_labels=['M','D','G3','G1','A1','G2','G4','A2','All'],colors = ['b','r','g'],
               size_x=12,size_y=8,y_lim=40.0,y_label='Somatic mutation',has_jitter=True)
