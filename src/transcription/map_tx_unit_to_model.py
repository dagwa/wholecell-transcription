import argparse
import csv
#import antimony # 2.7
import libsbml

parser = argparse.ArgumentParser(description='Map genes to model.')
parser.add_argument('input', metavar='I', type=str, nargs='+',
                    help='Input files')

args = parser.parse_args()

# Generate states per base pair?
perLocus = False

class MakeStateName(object):
  def __call__(self, desc, tx_unit, direction, active):
    '''
    Make ID name for given state of polymerase

    :param tx_unit: the transcription unit which the RNAp is bound to
    :param direction: the direction (which strand)
    :param active: whether the polymerase is active or not

    :returns: the identifier for the given RNAp state
    '''

    return '{}_tu{}_dir{}_act{}'.format(desc, tx_unit, direction, active)

class MakeStateNamePerLocus(MakeStateName):
  def __call__(self, desc, tx_unit, locus, direction, active):
    '''
    Make ID name for given state of polymerase

    :param tx_unit: the transcription unit which the RNAp is bound to
    :param direction: the direction (which strand)
    :param active: whether the polymerase is active or not

    :returns: the identifier for the given RNAp state
    '''

    return '{}_tu{}_loc_{}_dir{}_act{}'.format(desc, tx_unit, locus, direction, active)

class MakeStateNameFac:
  def __call__(self):
    if perLocus:
      return MakeStateName()
    else:
      return MakeStateNamePerLocus()

class RNApState(object):
  def __init__(self, tx_unit, direction, active, desc):
    self.tx_unit = tx_unit
    self.direction = direction
    self.active = active
    self.desc = desc

  def __repr__(self):
    return MakeStateName()(self.desc, self.tx_unit, self.direction, self.active)

  def __str__(self):
    return self.__repr__()

class RNApStatePerLocus(RNApState):
  def __init__(self, tx_unit, locus, direction, active, desc):
    super(RNApStatePerLocus, self).__init__(tx_unit, direction, active, desc)
    self.locus = locus

  def __repr__(self):
    return MakeStateNamePerLocus()(self.desc, self.tx_unit, self.locus, self.direction, self.active)


def direc_range():
  '''
  Range of values for polymerase direction
  '''
  return [0,1]

for i in args.input:
  with open(i) as f:
    # Read CSV
    reader = csv.reader(f)
    
    # Data columns
    namecol = 1
    lengthcol = 5

    # States for RNAp
    active_RNAp = []
    spec_bound_RNAp = []
    ns_bound_RNAp = []
    free_RNAp = []

    # Antimony model str
    antmdl = ''

    # Get gene names, properties
    for row in reader:
      try:
        tu_name = row[namecol].replace(' ', '_').replace('-', 'X')
        tu_length = int(row[lengthcol])
        #print(', '.join([row[0], tu_name, str(tu_length)]))

        def make_states(tx_unit, locus, direction, active):
          if locus is not None:
            active_state = RNApStatePerLocus(tx_unit = tu_name, locus = k, direction = direc, active = act, desc = 'active')
            spec_bound_state = RNApStatePerLocus(tx_unit = tu_name, locus = k, direction = direc, active = act, desc = 'spec_bound')
            ns_bound_state = RNApStatePerLocus(tx_unit = tu_name, locus = k, direction = direc, active = act, desc = 'ns_bound')
            free_state = RNApStatePerLocus(tx_unit = tu_name, locus = k, direction = direc, active = act, desc = 'free')
          else:
            active_state = RNApState(tx_unit = tu_name, direction = direc, active = act, desc = 'active')
            spec_bound_state = RNApState(tx_unit = tu_name, direction = direc, active = act, desc = 'spec_bound')
            ns_bound_state = RNApState(tx_unit = tu_name, direction = direc, active = act, desc = 'ns_bound')
            free_state = RNApState(tx_unit = tu_name, direction = direc, active = act, desc = 'free')

          active_RNAp.append(active_state)
          spec_bound_RNAp.append(spec_bound_state)
          ns_bound_RNAp.append(ns_bound_state)
          free_RNAp.append(free_state)

        for direc in direc_range():
          for act in [True, False]:
            if perLocus:
              for k in range(tu_length):
                make_states(tu_name, k, direc, act)
            else:
              make_states(tu_name, None, direc, act)
      except ValueError:
        # Discard header row etc.
        pass

    print('active states {}'.format(active_RNAp))