# swap_comment_character.py
import sys
import shutil
import os
import fnmatch

cards = {
         'LIQUID_EQUATION_TOLERANCE' : 'LIQUID_RESIDUAL_INFINITY_TOL',
         'GAS_EQUATION_TOLERANCE' : 'GAS_RESIDUAL_INFINITY_TOL',
         'LIQUID_PRESSURE_TOLERANCE' : 'MAX_ALLOW_REL_LIQ_PRES_CHANG_NI',
         'GAS_SATURATION_TOLERANCE' : 'MAX_ALLOW_REL_GAS_SAT_CHANGE_NI',
         'SATURATION_REL_PERTURBATION' : 'REL_GAS_SATURATION_PERTURBATION',
         'PRESSURE_REL_PERTURBATION' : 'REL_LIQ_PRESSURE_PERTURBATION',
         'SATURATION_MIN_PERTURBATION' : 'MIN_GAS_SATURATION_PERTURBATION',
         'PRESSURE_MIN_PERTURBATION' : 'MIN_LIQ_PRESSURE_PERTURBATION',
         'SATLIMIT' : 'GAS_SAT_THRESH_FORCE_EXTRA_NI',
         'DSATLIM' : 'GAS_SAT_THRESH_FORCE_TS_CUT',
         'DPRELIM' : 'MIN_LIQ_PRES_FORCE_TS_CUT',
         'DSAT_MAX' : 'MAX_ALLOW_GAS_SAT_CHANGE_TS',
         'DPRES_MAX' : 'MAX_ALLOW_LIQ_PRES_CHANGE_TS',
         'SATNORM' : 'GAS_SAT_CHANGE_TS_GOVERNOR',
         'PRESNORM' : 'LIQ_PRES_CHANGE_TS_GOVERNOR',
         'TSWITCH' : 'GAS_SAT_GOV_SWITCH_ABS_TO_REL',
         'LSCALE' : 'SCALE_JACOBIAN',
         'DO_NOT_LSCALE' : 'DO_NOT_SCALE_JACOBIAN',
         'P_SCALE' : 'JACOBIAN_PRESSURE_DERIV_SCALE',
         'MAXIT' : 'MAXIMUM_NUMBER_OF_ITERATIONS',
         'MAX_TS_CUTS' : 'MAXIMUM_CONSECUTIVE_TS_CUTS',
         'MAX_STEPS' : 'MAXIMUM_NUMBER_OF_TIMESTEPS',
         'DELTMIN' : 'MINIMUM_TIMESTEP_SIZE',
         'TIMESTEP_OVERSTEP_TOLERANCE' : 'TIMESTEP_OVERSTEP_REL_TOLERANCE',
         'ICONVTEST' : 'CONVERGENCE_TEST BOTH'}

remove_card_list = {'FIX_UPWIND_DIRECTION',
                    'DTIMEMAX',
                    'EPS_SAT',
                    'EPS_PRES',
                    '!NO_CREEP_CLOSURE',
                    '!NO_FRACTURE'}


def swap_file(filename):
  f = open(filename,'r')
  f2 = open(filename+'.tmp','w')
  # initialize local variables
  time_step_growth_factor = '' 
  for line in f:
    line2 = line
    # replace cards
    if len(line.rstrip()) > 1:
      card = line.split()[0]
      if card in cards:
        line2 = line.replace(card,cards[card])
      if card in remove_card_list:
        if card.startswith('DTIMEMAX'):
          time_step_growth_factor = line.split()[1]
        continue
    # custom replacement
    line_strip = line.lstrip()
    if line_strip.startswith('MODE BRAGFLO'):
      line2 = line.replace('BRAGFLO','WIPP_FLOW')
    if line_strip.startswith('ICONVTEST'):
      line2 = '        CONVERGENCE_TEST BOTH\n'
    if line2.strip().startswith('END_SUBSURFACE') and \
      len(time_step_growth_factor) > 0:
      # we have reached the end of the SUBSURFACE block without a 
      # TIMESTEPPER FLOW block, and we need one to properly set 
      # TIMESTEP_MAXIMUM_GROWTH_FACTOR
      f2.write('TIMESTEPPER FLOW\n')
      f2.write('  TIMESTEP_MAXIMUM_GROWTH_FACTOR %s\n' % \
               time_step_growth_factor)
      f2.write('END\n')
      # set back to blank
      time_step_growth_factor = ''
    f2.write(line2)
    # add select cards
    if line2.strip().startswith('TIMESTEPPER') and \
       line.find('FLOW') > -1:
      if len(time_step_growth_factor) > 0:
        f2.write('  TIMESTEP_MAXIMUM_GROWTH_FACTOR %s\n' % \
                 time_step_growth_factor)
        # set back to blank
        time_step_growth_factor = ''
  f.close()
  f2.close()
  # using shutil.move adds ^M to end of lines.
  os.remove(filename)
  shutil.copy(filename+'.tmp',filename)
  os.remove(filename+'.tmp')

suffix = '*.in'
for root, dirnames, filenames in os.walk('.'):
  for filename in fnmatch.filter(filenames,suffix):
    filename = os.path.join(root,filename)
    print(filename)
    swap_file(filename)

print('done')
