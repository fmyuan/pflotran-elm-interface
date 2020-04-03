# swap_comment_character.py
import sys
import shutil
import os
import fnmatch

timestepper_flow_keywords = [
        'PRESSURE_CHANGE_GOVERNOR',
        'TEMPERATURE_CHANGE_GOVERNOR',
        'CFL_GOVERNOR',
        'CONCENTRATION_CHANGE_GOVERNOR',
        'SATURATION_CHANGE_GOVERNOR',
        'GAS_SAT_CHANGE_TS_GOVERNOR',
        'LIQ_PRES_CHANGE_TS_GOVERNOR',
        'GAS_SAT_GOV_SWITCH_ABS_TO_REL',
        'MINIMUM_TIMESTEP_SIZE'
        ]

timestepper_trans_keywords = [
        'CFL_GOVERNOR',
        'VOLUME_FRACTION_CHANGE_GOVERNOR'
        ]  

newton_flow_keywords = [
        'NUMERICAL_JACOBIAN',
        'ANALYTICAL_JACOBIAN',
        'MAX_NEWTON_ITERATIONS',
        'RESIDUAL_INF_TOL',
        'RESIDUAL_ABS_INF_TOL',
        'LIQUID_RESIDUAL_ABS_INF_TOL',
        'GAS_RESIDUAL_ABS_INF_TOL',
        'ENERGY_RESIDUAL_ABS_INF_TOL',
        'RESIDUAL_SCALED_INF_TOL',
        'ITOL_SCALED_RESIDUAL',
        'LIQUID_RESIDUAL_SCALED_INF_TOL',
        'GAS_RESIDUAL_SCALED_INF_TOL',
        'ENERGY_RESIDUAL_SCALED_INF_TOL',
        'UPDATE_INF_TOL',
        'ABS_UPDATE_INF_TOL',
        'PRES_ABS_UPDATE_INF_TOL',
        'TEMP_ABS_UPDATE_INF_TOL',
        'SAT_ABS_UPDATE_INF_TOL',
        'XMOL_ABS_UPDATE_INF_TOL',
        'LIQUID_PRES_ABS_UPDATE_INF_TOL',
        'GAS_PRES_ABS_UPDATE_INF_TOL',
        'AIR_PRES_REL_UPDATE_INF_TOL',
        'PRESSURE_DAMPENING_FACTOR',
        'MAX_ITERATION_BEFORE_DAMPING',
        'DAMPING_FACTOR',
        'ITOL_RELATIVE_UPDATE',
        'PRESSURE_CHANGE_LIMIT',
        'SATURATION_CHANGE_LIMIT',
        'TEMPERATURE_CHANGE_LIMIT',
        'USE_INFINITY_NORM_CONVERGENCE',
        'LIQUID_RESIDUAL_INFINITY_TOL',
        'GAS_RESIDUAL_INFINITY_TOL',
        'MAX_ALLOW_REL_LIQ_PRES_CHANG_NI',
        'MAX_ALLOW_REL_GAS_SAT_CHANGE_NI',
        'GAS_SAT_THRESH_FORCE_TS_CUT',
        'MIN_LIQ_PRES_FORCE_TS_CUT',
        'MAX_ALLOW_GAS_SAT_CHANGE_TS',
        'MAX_ALLOW_LIQ_PRES_CHANGE_TS',
        'REL_LIQ_PRESSURE_PERTURBATION',
        'MIN_LIQ_PRESSURE_PERTURBATION',
        'REL_GAS_SATURATION_PERTURBATION',
        'MIN_GAS_SATURATION_PERTURBATION',
        'CONVERGENCE_TEST',
        'JACOBIAN_PRESSURE_DERIV_SCALE',
        'SCALE_JACOBIAN',
        'DO_NOT_SCALE_JACOBIAN',
        'GAS_SAT_THRESH_FORCE_EXTRA_NI',
        'MAXIMUM_PRESSURE_CHANGE'
        ]

newton_trans_keywords = [
        'NUMERICAL_JACOBIAN',
        'ANALYTICAL_JACOBIAN',
        'ITOL_RELATIVE_UPDATE'
        ]
  
def write_list(f,list_):
    if len(list_) > 0:
        for line in list_:
            f.write('  '+line.rstrip()+'\n')
        f.write('\n') 
        
def write_block(f,process,list_):    
    if list_[0].strip() != process:
        list_.insert(0,'{}\n'.format(process))
    pop_indexes = []
    for i in range(1,len(list_)):
        if list_[i].strip() == process:
            pop_indexes.append(i)
        else:
            list_[i] = '  ' + list_[i]
    for n in pop_indexes:
        list_.pop(n)
        
    list_.append('/\n')
    # do not write if empty
    if len(list_) > 2:
        write_list(f,list_)

def refactor_file(filename,replace_file_flag):

    f = open(filename,'r')
    
    linear_solver_flow = []
    newton_solver_flow = []
    timestepper_flow = []
    
    linear_solver_tran = []
    newton_solver_tran = []
    timestepper_tran = []
    
    tupl = ('TIMESTEPPER','NEWTON_SOLVER','LINEAR_SOLVER')
    
    store_mode = 0
    subsurface_card_found = False
    flow_process_model_found = False
    tran_process_model_found = False
    skip_count = 0
    while True:
        line = f.readline()

        if len(line) == 0:
            break
        w = line.split()

        card = ''
        if len(w) > 0:
          card = w[0].strip().upper()

        if card == 'SKIP':
            skip_count += 1

        if card == 'NOSKIP':
            skip_count -= 1

        if skip_count > 0:
            continue

        if card.startswith('SUBSURFACE_FLOW'):
            flow_process_model_found = True
            store_mode = 1

        elif card.startswith('SUBSURFACE_TRANSPORT'):
            tran_process_model_found = True
            store_mode = 2

        elif card.startswith('SUBSURFACE'):
            subsurface_card_found = True

        elif card.startswith('CPR_OPTIONS'):
                # skip file
                return

        elif card.startswith('EXTERNAL_FILE'):
            if not subsurface_card_found:
                # skip file
                return

        elif card.startswith(tupl):
            store_mode = 3

            list_ = []
            block_card = card
            
            if len(w) > 1:
                card2 = w[1].strip().upper()
            else:
                card2 = ""
            if card2 == 'FLOW':
                process_model = 1
                line = '{} \n'.format(card)

            elif card2 == 'TRANSPORT':
                process_model = 2
                line = '{} \n'.format(card)
            else:
                process_model = 0
        
        if store_mode == 1:
            
            if card in timestepper_flow_keywords:
                timestepper_flow.append('{} \n'.format(line.strip()))
                    
            if card in newton_flow_keywords:
                newton_solver_flow.append('{} \n'.format(line.strip()))
                    
        elif store_mode == 2:
            if card in timestepper_trans_keywords:
                timestepper_tran.append('{} \n'.format(line.strip()))

            if card in newton_trans_keywords:
                newton_solver_tran.append('{} \n'.format(line.strip()))                
        elif store_mode > 0 and \
             (card.startswith('END') or card.startswith('/')):
            store_mode = 0
            if block_card.startswith('TIMESTEPPER'):
                if process_model == 1:
                    timestepper_flow = list_ + timestepper_flow
                elif process_model == 2:
                    timestepper_tran = list_ + timestepper_tran
            elif block_card.startswith('NEWTON_SOLVER'):
                if process_model == 1:
                    newton_solver_flow = list_ + newton_solver_flow
                elif process_model == 2:
                    newton_solver_tran = list_ + newton_solver_tran
            elif block_card.startswith('LINEAR_SOLVER'):
                if process_model == 1:
                    linear_solver_flow = list_ + linear_solver_flow
                elif process_model == 2:
                    linear_solver_tran = list_ + linear_solver_tran
        
        if store_mode == 3:
            list_.append('{} \n'.format(line.strip()))
            

    f.seek(0)
    f2 = open(filename+'.tmp','w')

    skip_mode = 0
    skip_count = 0   # this integer applies to the SKIP/NOSKIP card stack
    insert_mode = False
    while True:
        line = f.readline()
        if len(line) == 0:
            break
        w = line.split()
        card = ''
        if len(w) > 0:
            card = w[0].strip().upper()

        if card == 'SKIP':
            skip_count += 1

        if card == 'NOSKIP':
            skip_count -= 1

        if skip_count > 0:
            f2.write(line)
            continue

        if card.startswith(tupl):
            skip_mode = 1
        elif card.startswith('SUBSURFACE_FLOW') or \
             card.startswith('SUBSURFACE_TRANSPORT'):
            skip_mode = 2
        if skip_mode == 0:
            f2.write(line)
        if skip_mode == 2:
            if card not in newton_flow_keywords and \
               card not in timestepper_flow_keywords and \
               card not in newton_trans_keywords and \
               card not in timestepper_trans_keywords:
                   f2.write(line)
                   
        if card.startswith('END') or card.startswith('/'):
            skip_mode = 0

        if card == 'SUBSURFACE':
            insert_mode = True

        if insert_mode:
            write_flow = flow_process_model_found and \
                         (timestepper_flow or newton_solver_flow or \
                          linear_solver_flow)
            write_tran = tran_process_model_found and \
                         (timestepper_tran or newton_solver_tran or \
                          linear_solver_tran)
            if write_flow or write_tran:
                    
                f2.write('\n#=========================== numerical methods '
                     '================================\n')
                if write_flow:
                    f2.write('NUMERICAL_METHODS FLOW\n\n')

                    if timestepper_flow:
                        write_block(f2,'TIMESTEPPER',timestepper_flow)
                    if newton_solver_flow:
                        write_block(f2,'NEWTON_SOLVER',newton_solver_flow)   
                    if linear_solver_flow:
                        write_block(f2,'LINEAR_SOLVER',linear_solver_flow)
                    f2.write('END\n\n')

                if write_tran:
                    f2.write('NUMERICAL_METHODS TRANSPORT\n\n')
                    if timestepper_tran:
                        write_block(f2,'TIMESTEPPER',timestepper_tran)
                    if newton_solver_tran:
                        write_block(f2,'NEWTON_SOLVER',newton_solver_tran)
                    if linear_solver_tran:
                        write_block(f2,'LINEAR_SOLVER',linear_solver_tran)
                    f2.write('END\n\n')
            insert_mode = False
    f.close()
    f2.close()
    
    if replace_file_flag:
        print('File {} has been updated.'.format(filename))
        # using shutil.move adds ^M to end of lines.
        os.remove(filename)
        shutil.copy(filename+'.tmp',filename)
    os.remove(filename+'.tmp')

replace_file_flag = True


def main():
#    filename = '1d_flux.in'
    suffix = '*.in'
    for root, dirnames, filenames in os.walk('.'):
        for filename in fnmatch.filter(filenames,suffix):
            filename = os.path.join(root,filename)
            print(filename)
            refactor_file(filename,replace_file_flag)


if __name__ == "__main__":
    try:
        suite_status = main()
        print("success")
        sys.exit(suite_status)
    except Exception as error:
        print(str(error))
#        if cmdl_options.backtrace:
#            traceback.print_exc()
#        traceback.print_exc()
        print("failure")
        sys.exit(1)

print('done')
