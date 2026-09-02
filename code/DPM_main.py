from DPM_run import *
from DPM_plot import DPM_plot_pho_cumulative

# Run flags
RUN_GENERATE_R1R2 = False
RUN_SIMULATION_WITHOUT_MISSPECIFICATION = False
RUN_SIMULATION_WITH_MISSPECIFICATION = True
RUN_PREPROCESSING_WITHOUT_MISSPECIFICATION = False
RUN_PREPROCESSING_WITH_MISSPECIFICATION = False
RUN_ANALYSIS_LOD = False
RUN_RUN_PLOT_1PAR = False
RUN_FIGURE1 = False
RUN_PROCESSING_MISSPECIFY_POP = False


# Generate R1 and R2 drug-resistant subclones for virtual patients based on the probability density function (PDF)
if RUN_GENERATE_R1R2:
    DPM_run_generate_R1R2_pdftrue()

# Run simulations without misspecification
if RUN_SIMULATION_WITHOUT_MISSPECIFICATION:
    parameterset = "data_ori"  # "data_pdf"   #
    filename_pattern = '_para'
    pathload = f'./{parameterset}/'
    pathsave = f'./result/{parameterset}'
    DPM_run_par_csv_folder(pathload=pathload,
                           pathsave=pathsave,
                           filename_pattern=filename_pattern,
                           run_sim=True,
                           Strategy_name=['strategy0', 'strategy2.2'],
                           fullinput=True,
                           save_filename_param=True,
                           save_filename_stopt=True,
                           save_filename_dosage=True,
                           save_filename_pop=True,
                           save_filename_eachtimepoint=False,
                           misspecification_sim_only=False,
                           use_parallel=False)

# Run simulations with misspecification
if RUN_SIMULATION_WITH_MISSPECIFICATION:
    parameterset = "data_ori"   # "data_pdf"
    filename_pattern = '_para'
    subclone_LOD = [1e-06, 1e-05, 1e-04, 1e-03, 1e-02, 1e-01]
    misspecification_LOD = [0, 'max', 'loguni', 'pdf']
    pathload = f'./{parameterset}/'
    for i_subclone_LOD in subclone_LOD:
        for i_misspecification_LOD in misspecification_LOD:
            pathsave = f'./result/{parameterset}_mis_LOD/{i_subclone_LOD:.0e}/{i_misspecification_LOD}/'
            DPM_run_par_csv_folder(pathload=pathload,
                                   pathsave=pathsave,
                                   filename_pattern=filename_pattern,
                                   subclone_LOD=i_subclone_LOD,
                                   misspecification_LOD=i_misspecification_LOD,
                                   Strategy_name=['strategy0', 'strategy2.2'],
                                   save_filename_param=True,
                                   save_filename_stopt=True,
                                   save_filename_dosage=True,
                                   save_filename_pop=False,
                                   save_filename_eachtimepoint=False,
                                   misspecification_sim_only=True,
                                   use_parallel=False)

# Run preprocessing without misspecification
if RUN_PREPROCESSING_WITHOUT_MISSPECIFICATION:
    parameterset = 'data_ori'  # 'data_pdf'
    pathload = f'./result/{parameterset}/'
    pathsave = f'./result/{parameterset}/figure'
    DPM_run_processing(Num_drug=2,
                       pathload=pathload,
                       pathsave=pathsave,
                       Num_stepdiff=2,
                       Strategy_name=['strategy0', 'strategy2.2'],
                       pathloadmis='',
                       use_parallel=False)

# Run preprocessing with misspecification
if RUN_PREPROCESSING_WITH_MISSPECIFICATION:
    parameterset = 'data_ori'  # 'data_pdf'
    pathload = f'./result/{parameterset}/'
    pathsave = f'./result/{parameterset}/figure'
    for i_LOD in ['1e-06', '1e-05', '1e-04', '1e-03', '1e-02', '1e-01']:
        pathloadmis = f'./result/{parameterset}_mis_LOD/' + i_LOD
        DPM_run_processing(Num_drug=2,
                           pathload=pathload,
                           pathsave=pathsave,
                           Num_stepdiff=2,
                           Strategy_name=['strategy0', 'strategy2.2'],
                           pathloadmis=pathloadmis,
                           use_parallel=False)


# Run the analysis for the simualtion results under both misspecification and non-misspecification scenarios
if RUN_ANALYSIS_LOD:
    filename_pattern = 'result'
    LOD = ['1e-06', '1e-05', '1e-04', '1e-03', '1e-02', '1e-01']
    for i_parameterset in ['data_ori', 'data_pdf']:
        pathload = f'./result/{i_parameterset}/'
        pathsave = f'./result/{i_parameterset}/figure'
        pathload_misLOD = f'./result/{i_parameterset}_mis_LOD/'
        DPM_run_analysis_LOD(Num_drug=2,
                             LOD=LOD,
                             filename_pattern=filename_pattern,
                             pathload=pathload,
                             pathsave=pathsave,
                             pathloadmis=pathload_misLOD,
                             Strategy_name=['strategy0', 'strategy2.2'],
                             plot=True)

# Run a single virtual patient example to demonstrate how subclone misspecification can worsen the DPM stratey. The required para_sel.pckl file
# is generated by DPM_run_analysis_LOD, which must be run first. Used for FigS4 and FigS5.
if RUN_RUN_PLOT_1PAR:
    subclone_LOD = 1e-2
    for i_parameterset in ['data_ori', 'data_pdf']:
        pathload = f'./result/{i_parameterset}/figure/example/para_sel.pckl'
        with bz2.BZ2File(pathload, 'rb') as f:
            par = pickle.load(f)
        subclone_LOD_str = f'{subclone_LOD:.0e}'
        subclone_LOD_str = subclone_LOD_str.replace('e-0', 'e-')
        pathsave = f'./result/{i_parameterset}/figure/example/{subclone_LOD_str}/'
        DPM_run_plot_1PAR(pathsave=pathsave,
                          par=par,
                          Strategy_name=['strategy0', 'strategy2.2'],
                          subclone_LOD=subclone_LOD,
                          misspecification_LOD=0)


# Run the cell population analysis.
if RUN_PROCESSING_MISSPECIFY_POP:
    parameterset = 'data_ori'   # 'data_pdf'
    pathload = f'./result/{parameterset}/'
    pathsave = f'./result/{parameterset}/figure'
    ylim = [500, 1000, 1000, 1000, 10000, 100000]
    for i, i_LOD in enumerate(['1e-06', '1e-05', '1e-04', '1e-03', '1e-02', '1e-01']):
        pathloadmis = f'./result/{parameterset}_mis_LOD/' + i_LOD
        DPM_run_processing_misspecfiy_pop(pathload=pathload,
                                          pathsave=pathsave,
                                          Num_stepdiff=2,
                                          Strategy_name=['strategy0', 'strategy2.2'],
                                          pathloadmis=pathloadmis,
                                          ylim=ylim[i])

# For generating Figure 1 in the manuscript.
if RUN_FIGURE1:
    pathsave = './figure'
    DPM_plot_pho_cumulative(pathsave=pathsave)