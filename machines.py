import os
import sys
machines = {}

def add_machine(name, compiler, c_flags, l_flags, gsl_dir):
  machine = {}
  machine['NAME'] = name
  machine['COMPILER'] = compiler
  machine['COMPILER_FLAGS'] = c_flags
  machine['LIB_FLAGS'] = l_flags
  machine['GSL_DIR'] = gsl_dir
  machines[name] = machine

add_machine(name='meade', 
            compiler='h5pcc',  
            c_flags='-Ofast -Wall -Werror -fdiagnostics-color -fopenmp',
            l_flags='',
            gsl_dir='/home/brryan/Software/gsl')

def get_machine():
  for key in machines:
    print key
    if os.uname()[1] == machines[key]['NAME']:
      return machines[os.uname()[1]]
  print 'UNKNOWN MACHINE ' + os.uname()[1]
  sys.exit()

