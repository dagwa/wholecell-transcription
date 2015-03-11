### reads model components from tab limited text file and generates an SMBL file
### Pnar Pir 2015

import libsbml
import sys
#import numpy

def main():  
#def main(inputfile,outputfile):

### Generate the model using the lists derived from tablim text file
	outputfile='theModel.xml'

	sbmlDoc = libsbml.SBMLDocument(2, 4)
	model = sbmlDoc.createModel()
	model.setId('aModel')

	# unitdef = model.createUnitDefinition()
	# unitdef.setId("particles_per_cell_per_second")

	# #  Creates an Unit inside the UnitDefinition object ("litre_per_mole_per_second")

	# # unit = unitdef.createUnit()
	# # unit.setKind(unit_kind_dimensionless)
	# # unit.setExponent(1)

	# #  Creates an Unit inside the UnitDefinition object ("litre_per_mole_per_second")

	# unit = unitdef.createUnit()
	# unit.setKind(UNIT_KIND_LITRE)
	# unit.setExponent(-1)

	# #  Creates an Unit inside the UnitDefinition object ("litre_per_mole_per_second")

	# unit = unitdef.createUnit()
	# unit.setKind(UNIT_KIND_SECOND)
	# unit.setExponent(-1)


### Read the model components into Python lists from tablim text file

#	texts = open('modelComponents.txt', 'r') # lists of model components
	texts = open('modelPromoters.txt', 'r') # lists of model components
	for line in range (0,1000):
		fline=texts.readline() # read a line from the file
		cols=fline.split('\t') # split the line into its columns
		headerfound=cols[0].find('Header')
		if (headerfound<0) :
			compartmentfound=cols[0].find('Compartment')
			if (compartmentfound<0) :			
				speciesfound=cols[0].find('Species')
				if (speciesfound<0) :
					globalvariablefound=cols[0].find('GlobalVariable')
					if (globalvariablefound<0) :
						reactionfound=cols[0].find('Reaction')
						if (reactionfound<0) :
							print 'no reaction'
							eventfound=cols[0].find('Event')
							if (eventfound<0) :
								print 'either the last line is read or something went wrong'
							else:
	#							addEvent(model,cols)
								print 'will add event'
						else:
							#addReaction(model,cols)
							print 'will add reaction'
							rxn = model.createReaction()
							rxn.setId(cols[1])
							#rxn.setReversible(bool(cols[2])) 
							rxn.setReversible(bool('false')) 
							if (len(cols[3])>2):
								reactants=cols[3].replace('"','').split(',')
								print len(reactants)
								print reactants
								#daat
								for rs in range (0,len(reactants)):
									r=rxn.createReactant()
									r.setSpecies(reactants[rs])
									r.setStoichiometry(1) 
							if (len(cols[4])>2):
								p=rxn.createProduct()
								p.setSpecies(cols[4])
								p.setStoichiometry(1) 	

							par=model.createParameter()
							par.setId('k_'+reactants[1]+'_binds_'+reactants[0])
							par.setValue(0.0002)

							kinetic_string=str('k_'+reactants[1]+'_binds_'+reactants[0]+'*'+reactants[0]+'*'+reactants[1])	
							KL=rxn.createKineticLaw()
							aa=libsbml.parseFormula(kinetic_string)
							KL.setMath(aa)
					else:
	#					addGlobalVariable(model,cols)
						print 'adding the parameter', cols[2]
						par=model.createParameter()
						par.setId(cols[2])
						par.setName(cols[1])
						par.setValue(float(cols[3]))
						print 'lencol', len(cols[5])
						if (len(cols[5])>2):
							ra=model.createInitialAssignment()
							ra.setSymbol(cols[2])
							bb=libsbml.parseL3Formula(cols[5])
							ra.setMath(bb)	
						if (len(cols[6])>2):
							par.setConstant(bool(0))
							ra=model.createAssignmentRule()
							ra.setVariable(str(cols[2]))
							bb=libsbml.parseL3Formula(cols[6])
							ra.setMath(bb)  							

				else:
	#				addSpecies(model,cols)
					print 'adding the species', cols[2]
					sp=model.createSpecies()
					sp.setId(cols[2]) 
					sp.setName(cols[1]) 
					sp.setCompartment(cols[3])
					sp.setInitialConcentration(float(cols[4]))
					print 'lencol', len(cols[6])
					if (len(cols[6])>2):
						ra=model.createInitialAssignment()
						ra.setSymbol(str(cols[2]))
						bb=libsbml.parseL3Formula(cols[6])
						ra.setMath(bb)  
					if (len(cols[7])>2):
						ra=model.createAssignmentRule()
						ra.setVariable(str(cols[2]))
						bb=libsbml.parseL3Formula(cols[7])
						ra.setMath(bb)  	

			else:
	#			addCompartment(model,cols)
				print 'adding the comparment', cols[2]
				com = model.createCompartment()
				com.setId(cols[1])
				com.setSize(float(cols[3]))

	libsbml.writeSBMLToFile(sbmlDoc, str(outputfile))			

if(__name__ == '__main__'):
	main() 


