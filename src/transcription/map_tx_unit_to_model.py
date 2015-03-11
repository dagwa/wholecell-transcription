import argparse
import csv
#import antimony # 2.7
import libsbml

parser = argparse.ArgumentParser(description='Map genes to model.')
#parser.add_argument('input', metavar='I', type=str, nargs='+',
                    #help='Input files')
parser.add_argument('--tx-units', dest='tx_units', nargs=1, help='Transcription units CSV')
parser.add_argument('--genes', dest='genes', nargs=1, help='Genes units CSV')
parser.add_argument('--locus', dest='locus', action='store_const',
                   const=True, default=False,
                   help='Generate per-locus')

args = parser.parse_args()

tu_len_cutoff = 3

# Generate states per base pair?
perLocus = args.locus

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
  def __init__(self, desc):
    self.desc = desc

  def __repr__(self):
    return self.id()

  def __str__(self):
    return self.__repr__()

  def id(self):
    return self.desc

class RNApBoundState(RNApState):
  pass

class RNApSpecBoundState(RNApBoundState):
  def __init__(self, tx_unit, direction, active, desc):
    self.tx_unit = tx_unit
    self.direction = direction
    self.active = active
    self.desc = desc

  def id(self):
    return MakeStateName()(self.desc, self.tx_unit, self.direction, self.active)

  def __str__(self):
    return self.__repr__()

  def global_locus(self):
    return self.tx_unit

class RNApSpecBoundStatePerLocus(RNApSpecBoundState):
  def __init__(self, tx_unit, locus, direction, active, desc):
    super(RNApSpecBoundStatePerLocus, self).__init__(tx_unit, direction, active, desc)
    self.locus = locus

  def id(self):
    return MakeStateNamePerLocus()(self.desc, self.tx_unit, self.locus, self.direction, self.active)

  def global_locus(self):
    return '_'.join(self.tx_unit, str(self.locus))


def direc_range():
  '''
  Range of values for polymerase direction
  '''
  return [0,1]

gene_seq = {}
# Construct gene sequences
with open(args.genes[0]) as gene_f:
  gene_reader = csv.reader(gene_f)

  # Columns for data
  name_col = 0
  seq_col = 22

  # Read rows
  for row in gene_reader:
    # Get name and sequence of gene
    name = row[name_col]
    seq = row[seq_col]

    gene_seq[name] = seq

# States for RNAp

class RNAp_states:
  def __init__(self):
    self.active_RNAp = []
    self.active_sigma_bound_RNAp = []
    self.spec_bound_RNAp = []
    self.ns_bound_RNAp = RNApBoundState('RNAp_NonSpecBound')
    self.free_RNAp = RNApState('RNAp_Free')

  def index(self):
    self.active_state_for_unit = {}
    self.active_state_sigma_bound_for_unit = {}
    self.spec_bound_for_unit = {}

    for state in self.active_RNAp:
      self.active_state_for_unit[state.global_locus()] = state
    for state in self.active_sigma_bound_RNAp:
      self.active_state_sigma_bound_for_unit[state.global_locus()] = state
    for state in self.spec_bound_RNAp:
      self.spec_bound_for_unit[state.tx_unit] = state

  def map_spec_bound_to_active_state(self, state):
    return self.active_state_for_unit[state.global_locus()]

  def add_active_RNAp_state(self, state):
    self.active_RNAp.append(state)

  def add_active_sigma_bound_RNAp_state(self, state):
    self.active_sigma_bound_RNAp.append(state)

  def add_spec_bound_RNAp_state(self, state):
    self.spec_bound_RNAp.append(state)

  # Allow iterating over states
  def __iter__(self):
    return iter(self.active_RNAp + self.spec_bound_RNAp + self.active_sigma_bound_RNAp + [self.ns_bound_RNAp] + [self.free_RNAp])

# Make specifically bound and active states for RNAp on each operon
def make_states(states, tx_unit, locus, direction, active):
  #states: RNAp_states

  if locus is not None:
    # Locus granularity
    active_state = RNApSpecBoundStatePerLocus(tx_unit = tu_name, locus = k, direction = direc, active = act, desc = 'active')
    active_state_sigma_bound = RNApSpecBoundStatePerLocus(tx_unit = tu_name, locus = k, direction = direc, active = act, desc = 'active_sigma_bound')
    spec_bound_state = RNApSpecBoundStatePerLocus(tx_unit = tu_name, locus = k, direction = direc, active = act, desc = 'spec_bound')
  else:
    # Operon granularity
    active_state = RNApSpecBoundState(tx_unit = tu_name, direction = direc, active = act, desc = 'active')
    active_state_sigma_bound = RNApSpecBoundState(tx_unit = tu_name, direction = direc, active = act, desc = 'active_sigma_bound')
    spec_bound_state = RNApSpecBoundState(tx_unit = tu_name, direction = direc, active = act, desc = 'spec_bound')

  # Append to lists of all states
  states.add_active_RNAp_state(active_state)
  states.add_active_sigma_bound_RNAp_state(active_state_sigma_bound)
  states.add_spec_bound_RNAp_state(spec_bound_state)


#def RNAp_states():
  #return active_RNAp + spec_bound_RNAp + active_sigma_bound_RNAp + [ns_bound_RNAp] + [free_RNAp]

states = RNAp_states()

#for i in args.input:
with open(args.tx_units[0]) as tx_f:
  # Read CSV
  tu_reader = csv.reader(tx_f)

  # Data columns
  namecol = 1
  genescol = 3
  lengthcol = 5

  # Antimony model str
  antmdl = ''

  # Get gene names, properties
  for row in tu_reader:
    try:
      tu_name = row[namecol].replace(' ', '_').replace('-', 'X').replace(',', '_')
      # Genes transcribed
      tu_genes = [x.strip() for x in row[genescol].split(',')]
      # Length of transcription unit
      tu_length = int(row[lengthcol])

      # Join genes to form transcription unit sequence
      tu_seq = ''.join([x for x in [gene_seq[g] for g in tu_genes]])
      tu_seq = tu_seq.replace('T', 'U') # DNA seq -> RNA seq
      #print(tu_seq)

      # Combinatorics for each param
      for direc in direc_range():
        for act in [True, False]:
          if perLocus:
            for k in range(min(tu_length, tu_len_cutoff)):
              make_states(states, tu_name, k, direc, act)
          else:
            make_states(states, tu_name, None, direc, act)
    except ValueError:
      # Discard header row etc.
      pass

  #print('active states {}'.format(active_RNAp))

# create SBML model

def check(value, message):
  """If 'value' is None, prints an error message constructed using
  'message' and then exits with status code 1. If 'value' is an integer,
  it assumes it is a libSBML return status code. If the code value is
  LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
  prints an error message constructed using 'message' along with text from
  libSBML explaining the meaning of the code, and exits with status code 1.
  """
  if value == None:
    raise RuntimeError('LibSBML returned a null value trying to ' + message + '.')
  elif type(value) is int:
    if value == libsbml.LIBSBML_OPERATION_SUCCESS:
      return
    else:
      err_msg = 'Error encountered trying to ' + message + '.' \
        + 'LibSBML returned error code ' + str(value) + ': "' \
        + libsbml.OperationReturnValue_toString(value).strip() + '"'
      raise RuntimeError(err_msg)
  else:
    return

states.index()

# http://sbml.org/Software/libSBML/docs/python-api/libsbml-python-creating-model.html

document = libsbml.SBMLDocument(3, 1)
model = document.createModel()
check(model, 'create model')
check(model.setTimeUnits("second"), 'set model-wide time units')
check(model.setExtentUnits("mole"), 'set model units of extent')
check(model.setSubstanceUnits('mole'), 'set model substance units')

per_second = model.createUnitDefinition()
check(per_second, 'create unit definition')
check(per_second.setId('per_second'), 'set unit definition id')
unit = per_second.createUnit()
check(unit, 'create unit on per_second')
check(unit.setKind(libsbml.UNIT_KIND_SECOND), 'set unit kind')
check(unit.setExponent(-1), 'set unit exponent')
check(unit.setScale(0), 'set unit scale')
check(unit.setMultiplier(1), 'set unit multiplier')

c1 = model.createCompartment()
check(c1, 'create compartment')
check(c1.setId('c1'), 'set compartment id')
check(c1.setConstant(True), 'set compartment "constant"')
check(c1.setSize(1), 'set compartment "size"')
check(c1.setSpatialDimensions(3), 'set compartment dimensions')
check(c1.setUnits('litre'), 'set compartment size units')

# Create states
for RNAp_state in states:
  spec = model.createSpecies()
  check(spec, 'create species spec')
  idstr = RNAp_state.id()
  #print(idstr)
  check(spec.setId(idstr), 'set species spec id')
  check(spec.setCompartment('c1'), 'set species spec compartment')
  check(spec.setConstant(False), 'set "constant" attribute on spec')
  check(spec.setInitialAmount(0), 'set initial amount for spec')
  check(spec.setSubstanceUnits('item'), 'set substance units for spec')
  check(spec.setBoundaryCondition(False), 'set "boundaryCondition" on spec')
  check(spec.setHasOnlySubstanceUnits(False), 'set "hasOnlySubstanceUnits" on spec')

# Create reactions

# Create states spec_bound -> active
counter = 0
for state in states.spec_bound_RNAp:
  active_state = states.map_spec_bound_to_active_state(state)

  r = model.createReaction()
  check(r, 'create reaction')
  check(r.setId('r_spec_to_active_{}'.format(counter)), 'set reaction id')
  check(r.setReversible(False), 'set reaction reversibility flag')
  check(r.setFast(False), 'set reaction "fast" attribute')
  counter += 1

  reactant = r.createReactant()
  check(reactant, 'create reactant')
  check(reactant.setSpecies(state.id()), 'assign reactant species')
  check(reactant.setConstant(False), 'set "constant" on species ref 1')

  product = r.createProduct()
  check(product, 'create product')
  check(product.setSpecies(active_state.id()), 'assign product species')
  check(product.setConstant(False), 'set "constant" on species ref 2')

  math_ast = libsbml.parseL3Formula('1')
  check(math_ast, 'create AST for rate expression')

  kinetic_law = r.createKineticLaw()
  check(kinetic_law, 'create kinetic law')
  check(kinetic_law.setMath(math_ast), 'set math on kinetic law')

sbmlstr = libsbml.writeSBMLToString(document)
with open('/tmp/tx.sbml', 'w') as f:
  f.write(sbmlstr)