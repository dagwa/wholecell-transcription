from __future__ import print_function

import argparse
import csv
#import antimony # 2.7
import libsbml
from sets import Set as set

class MissingItemError(LookupError):
  pass

parser = argparse.ArgumentParser(description='Map genes to model.')
#parser.add_argument('input', metavar='I', type=str, nargs='+',
                    #help='Input files')
parser.add_argument('--tx-units', dest='tx_units', nargs=1, help='Transcription units CSV')
parser.add_argument('--genes', dest='genes', nargs=1, help='Genes units CSV')
parser.add_argument('--bp-cutoff', dest='bp_cutoff', type=int, nargs=1, default=[3], help='Maximum number of base pairs per transcription unit')
parser.add_argument('--tu-cutoff', dest='tu_cutoff', type=int, nargs=1, default=[None], help='Maximum number of transcription units to use')
parser.add_argument('--locus', dest='locus', action='store_const',
                   const=True, default=False,
                   help='Generate per-locus')

args = parser.parse_args()

# Cutoff length of a transcription unit
tu_len_cutoff = args.bp_cutoff[0]

# Max number of transcription units to include (None for all)
tu_cutoff = args.tu_cutoff[0]

# Generate states per base pair?
perLocus = args.locus

def make_global_locus(tx_unit, locus=None):
  if locus is not None:
    return '_'.join([tx_unit, str(locus)])
  else:
    return tx_unit

# Find an element in a container
def find(cont, cond):
  result = None
  for elt in cont:
    if cond(elt):
      if result is not None:
        print('Warning: multiply defined elt')
      result = elt
      #return elt
  #print('container: {}'.format(cont))
  if result is None:
    raise MissingItemError('No such element found')
  return result

# Make ID for RNAp state
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
    return make_global_locus(self.tx_unit)

  def matches_props(self, x):
    '''
    Return true if polymerase x is identical except for Markov state
    '''
    if isinstance(x, RNApSpecBoundStatePerLocus) and \
      self.tx_unit == x.tx_unit and self.direction == x.direction:
      #self.desc == x.desc:
        return True
    else:
      return False

class RNApSpecBoundStatePerLocus(RNApSpecBoundState):
  def __init__(self, tx_unit, locus, direction, active, desc):
    super(RNApSpecBoundStatePerLocus, self).__init__(tx_unit, direction, active, desc)
    self.locus = locus

  def id(self):
    return MakeStateNamePerLocus()(self.desc, self.tx_unit, self.locus, self.direction, self.active)

  def global_locus(self):
    return make_global_locus(self.tx_unit, self.locus)

  def matches_props(self, x):
    '''
    Return true if polymerase x is identical except for Markov state
    '''
    if isinstance(x, RNApSpecBoundStatePerLocus) and \
      self.global_locus() == x.global_locus() and self.direction == x.direction:
      #self.desc == x.desc:
        return True
    else:
      return False


# --------------
# Transcription Factors
# --------------

class SigmaState(object):
  def __init__(self, desc):
    self.desc = desc

  def id(self):
    return self.desc

class SigmaFactorBoundState(SigmaState):
  def __init__(self, tx_unit, direction, desc):
    self.tx_unit = tx_unit
    self.direction = direction
    self.desc = desc

  def id(self):
    return '{}_tu{}_dir{}'.format(self.desc, self.tx_unit, self.direction)

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

def add_state_to_map(map, locus_key, state):
  if not locus_key in map:
    map[locus_key] = set()
  map[locus_key].add(state)

# Collection of polymerase states
class RNAp_states:
  def __init__(self):
    self.active_RNAp = []
    #self.active_sigma_bound_RNAp = []
    self.spec_bound_RNAp = []
    self.ns_bound_RNAp = RNApBoundState('RNAp_NonSpecBound')
    self.free_RNAp = RNApState('RNAp_Free')

  def index(self):
    self.active_state_for_locus = {}
    self.active_state_sigma_bound_for_locus = {}
    self.spec_bound_for_unit = {}
    self.tx_units = set()
    # Initial active states of RNAp - transcript length 0
    self.init_active_states = set()

    for state in self.spec_bound_RNAp:
      self.tx_units.add(state.tx_unit)
      add_state_to_map(self.spec_bound_for_unit, state.tx_unit, state)
    for state in self.active_RNAp:
      add_state_to_map(self.active_state_for_locus, state.global_locus(), state)
      if state.locus == 0:
        self.init_active_states.add(state)
    #for state in self.active_sigma_bound_RNAp:
      #add_state_to_map(self.active_state_sigma_bound_for_locus, state.global_locus(), state)

  def map_spec_bound_to_active_state(self, state):
    return find(self.active_state_for_locus[make_global_locus(state.tx_unit, 0)], lambda x: state.matches_props(x) and x.locus == 0)

  #def map_tx_unit_to_active_state(self, tx_unit):
    #return find(self.active_state_for_locus[make_global_locus(tx_unit, 0)], lambda x: state.matches_props(x) and x.locus == 0)

  def get_next_elongation_state(self, state):
    #print('get_next_elongation_state: start {}'.format(state))
    #print('get_next_elongation_state: {}'.format(self.active_state_for_locus[make_global_locus(state.tx_unit, state.locus+1)]))
    try:
      return find(self.active_state_for_locus[make_global_locus(state.tx_unit, state.locus+1)], lambda x: state.tx_unit == x.tx_unit and state.direction == x.direction and x.locus == state.locus+1)
    except KeyError:
      raise MissingItemError('No further elongation')

  def add_active_RNAp_state(self, state):
    self.active_RNAp.append(state)

  #def add_active_sigma_bound_RNAp_state(self, state):
    #self.active_sigma_bound_RNAp.append(state)

  def add_spec_bound_RNAp_state(self, state):
    self.spec_bound_RNAp.append(state)

  # Allow iterating over states
  def __iter__(self):
    #return iter(self.active_RNAp + self.spec_bound_RNAp + self.active_sigma_bound_RNAp + [self.ns_bound_RNAp] + [self.free_RNAp])
    return iter(self.active_RNAp + self.spec_bound_RNAp + [self.ns_bound_RNAp] + [self.free_RNAp])

# Collection of transcription factor states
class TF_states:
  def __init__(self):
    self.bound_sigma = []
    self.free_sigma = SigmaState('Sigma_Free')

    self.tx_unit_to_simga_bound_state = {}

  def add_sigma_factor_bound(self, state):
    self.bound_sigma.append(state)
    add_state_to_map(self.tx_unit_to_simga_bound_state, state.tx_unit, state)

  def map_spec_bound_to_simga_bound_state(self, state):
    return find(self.tx_unit_to_simga_bound_state[state.tx_unit], lambda x: state.tx_unit == x.tx_unit and state.direction == x.direction)

  #def map_tx_unit_to_simga_bound_state(self, tx_unit):
    #return self.tx_unit_to_simga_bound_state[tx_unit]

  def index(self):
    self.bound_sigma_for_locus = {}
    for state in bound_sigma:
      self.bound_sigma_for_locus[state.tx_unit] = state

  # Allow iterating over states
  def __iter__(self):
    return iter(self.bound_sigma + [self.free_sigma])


# Make specifically bound and active states for RNAp on each operon
def make_states(RNAp_states, tf_states, tx_unit, tu_length, direction):
  #RNAp_states: RNAp_RNAp_states

  # Polymerase has active state per locus
  if tu_len_cutoff is not None:
    for k in range(min(tu_length, tu_len_cutoff)):
      # Locus granularity
      active_state = RNApSpecBoundStatePerLocus(tx_unit = tu_name, locus = k, direction = direc, active = True, desc = 'RNAp_active')
      #active_state_sigma_bound = RNApSpecBoundStatePerLocus(tx_unit = tu_name, locus = k, direction = direc, active = True, desc = 'active_sigma_bound')

      # Add active states
      RNAp_states.add_active_RNAp_state(active_state)
      #RNAp_states.add_active_sigma_bound_RNAp_state(active_state_sigma_bound)
  else:
    # Operon granularity
    active_state = RNApSpecBoundState(tx_unit = tu_name, direction = direc, active = True, desc = 'RNAp_active')
    #active_state_sigma_bound = RNApSpecBoundState(tx_unit = tu_name, direction = direc, active = True, desc = 'active_sigma_bound')

    # Add active states
    RNAp_states.add_active_RNAp_state(active_state)
    #RNAp_states.add_active_sigma_bound_RNAp_state(active_state_sigma_bound)

  spec_bound_state = RNApSpecBoundState(tx_unit = tu_name, direction = direc, active = False, desc = 'RNAp_spec_bound')

  # Append specifically bound statest to lists of all RNAp_states
  RNAp_states.add_spec_bound_RNAp_state(spec_bound_state)

  # Append transcription factor states
  tf_states.add_sigma_factor_bound(SigmaFactorBoundState(tx_unit = tu_name, direction = direc, desc = 'Sigma_Bound'))

rnap_states = RNAp_states()

tf_states = TF_states()

# Read transcription units and create states
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
  count = 0
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
          make_states(rnap_states, tf_states, tu_name, tu_length, direc)

      #print('read tu {}: {}'.format(count, tu_name))
      count += 1
      if tu_cutoff and count == tu_cutoff:
        break
    except ValueError:
      # Discard header row etc.
      pass

  #print('active rnap_states {}'.format(active_RNAp))

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

rnap_states.index()

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

# Create rnap_states
for RNAp_state in rnap_states:
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

# Create TF states
for tf_state in tf_states:
  spec = model.createSpecies()
  check(spec, 'create species spec')
  idstr = tf_state.id()
  #print(idstr)
  check(spec.setId(idstr), 'set species spec id')
  check(spec.setCompartment('c1'), 'set species spec compartment')
  check(spec.setConstant(False), 'set "constant" attribute on spec')
  check(spec.setInitialAmount(0), 'set initial amount for spec')
  check(spec.setSubstanceUnits('item'), 'set substance units for spec')
  check(spec.setBoundaryCondition(False), 'set "boundaryCondition" on spec')
  check(spec.setHasOnlySubstanceUnits(False), 'set "hasOnlySubstanceUnits" on spec')

# Create reactions

# Reaction: Activate specifically bound polymerase
counter = 0
for state in rnap_states.spec_bound_RNAp:
  active_state = rnap_states.map_spec_bound_to_active_state(state)
  sigma_bound_state = tf_states.map_spec_bound_to_simga_bound_state(state)

  r = model.createReaction()
  check(r, 'create reaction')
  check(r.setId('r_RNAp_spec_to_active_{}'.format(counter)), 'set reaction id')
  check(r.setReversible(False), 'set reaction reversibility flag')
  check(r.setFast(False), 'set reaction "fast" attribute')
  counter += 1

  # Create reactant: specifically bound polymerase
  reactant = r.createReactant()
  check(reactant, 'create reactant')
  check(reactant.setSpecies(state.id()), 'assign reactant species')
  check(reactant.setConstant(False), 'set "constant" on species ref 1')

  # Create product: active polymerase
  product = r.createProduct()
  check(product, 'create product')
  check(product.setSpecies(active_state.id()), 'assign product species')
  check(product.setConstant(False), 'set "constant" on species ref 2')

  # Create modifier: sigma factor
  sigma = r.createModifier()
  check(sigma, 'create modifier')
  check(sigma.setSpecies(sigma_bound_state.id()), 'assign product species')

  math_ast = libsbml.parseL3Formula(sigma_bound_state.id())
  check(math_ast, 'create AST for rate expression')

  kinetic_law = r.createKineticLaw()
  check(kinetic_law, 'create kinetic law')
  check(kinetic_law.setMath(math_ast), 'set math on kinetic law')

# Elongation reactions
for init_state in rnap_states.init_active_states:
  state = init_state
  counter = 0
  try:
    while True:
      nxt = rnap_states.get_next_elongation_state(state)
      #print('nxt: {}'.format(nxt))
      r = model.createReaction()
      check(r, 'create reaction')
      check(r.setId('r_RNAp_elongation_tu{}_d{}_{}'.format(init_state.tx_unit, init_state.direction, counter)), 'set reaction id')
      check(r.setReversible(False), 'set reaction reversibility flag')
      check(r.setFast(False), 'set reaction "fast" attribute')

      # Create reactant: polymerase at locus
      reactant = r.createReactant()
      check(reactant, 'create reactant')
      check(reactant.setSpecies(state.id()), 'assign reactant species')
      check(reactant.setConstant(False), 'set "constant" on species ref 1')

      # Create product: active polymerase
      product = r.createProduct()
      check(product, 'create product')
      check(product.setSpecies(nxt.id()), 'assign product species')
      check(product.setConstant(False), 'set "constant" on species ref 2')

      math_ast = libsbml.parseL3Formula('1')
      check(math_ast, 'create AST for rate expression')

      kinetic_law = r.createKineticLaw()
      check(kinetic_law, 'create kinetic law')
      check(kinetic_law.setMath(math_ast), 'set math on kinetic law')

      counter += 1
      state = nxt
  except MissingItemError:
    # At end of elongation
    pass

sbmlstr = libsbml.writeSBMLToString(document)
with open('/tmp/tx.sbml', 'w') as f:
  f.write(sbmlstr)