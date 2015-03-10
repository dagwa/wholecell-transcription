import pandas as pd
import libsbml
import re


def add_species_by_id(species_id, species_type, model, species_name=" ", compartment='c'):
    if species_id not in SBML_SPECIES:
        species = model.createSpecies()
        species.setId( species_id)
        species.setName( species_name)
        SBML_SPECIES.append(species_id)
        species.setCompartment(compartment)
        if species_type in ENTITY_TO_SBO_MAPPING.keys():
            species.setSBOTerm( ENTITY_TO_SBO_MAPPING[species_type])

def create_compartment(name_of_compartment, model):
    if name_of_compartment:
        if not model.getCompartment(name_of_compartment):
            NEW_COMPARTMENT = model.createCompartment()
            NEW_COMPARTMENT.setId(name_of_compartment)
            NEW_COMPARTMENT.setConstant(True)
            NEW_COMPARTMENT.setSize(1)
            NEW_COMPARTMENT.setUnits('volume')


reactionFile_data = pd.read_csv('Reactions.txt', delimiter='\t')
relevantReactions = reactionFile_data[reactionFile_data['Process'].isin(['Process_RNAModification','Process_RNADecay','Process_RNAProcessing','Process_Transcription' ])]

SBML_SPECIES = []


ENTITY_TO_SBO_MAPPING = {
    "Gene" : "SBO:0000354", # informational molecule segment
    "Complex" : "SBO:0000253", # non-covalent complex
    "Protein" : "SBO:0000252", # polypeptide chain
    "Dna" : "SBO:0000251", # deoxyribonucleic acid
    "DnaRegion" : "SBO:0000251", # deoxyribonucleic acid
    "Rna" : "SBO:0000250", # ribonucleic acid
    "RnaRegion" : "SBO:0000250", # ribonucleic acid
    "SmallMolecule" : "SBO:0000247", # simple chemical
    "Simple_molecule" : "SBO:0000247", # simple chemical
    "tRNA" : "SBO:0000313",
    "rRNA" : "SBO:0000314",
    "snRNA" : "SBO:0000318",
    "miRNA" : "SBO:0000316",
    "TSS": "SO:0000315",
    "mRNA": "SBO:0000278",
    "metabolite": "SBO:0000299"
}

DOCUMENT = libsbml.SBMLDocument( 2, 4)
MODEL = DOCUMENT.createModel()

Process_RNAModification_reactions = reactionFile_data[reactionFile_data['Process'].isin(['Process_RNAModification'])]

for index, row in Process_RNAModification_reactions.iterrows():
	if row['Molecule']:
		reaction_name = row['ID']
		print reaction_name
		m = re.search(r"\[([A-Za-z0]+)\]", row['Stoichiometry'])
		molecule = row['Molecule']
		modifier = row['Enzyme'] 
		modifier_name = row['Name']
		modifier_compartment = row['Enzyme Compartment'] 

		compartment = m.group(1)
		create_compartment(compartment, MODEL)
		create_compartment(modifier_compartment, MODEL)

		reaction_stoich = row['Stoichiometry'].rpartition(':')[-1].strip()
		if "<==>" in reaction_stoich:
			reversible = True
			LHS_reaction = reaction_stoich.split('<==>')[0].strip()
			RHS_reaction = reaction_stoich.split('<==>')[1].strip()
		elif "==>" in reaction_stoich:
			reversible = False
			LHS_reaction = reaction_stoich.split('==>')[0].strip()
			RHS_reaction = reaction_stoich.split('==>')[1].strip()
		reactants = [i.strip() for i in LHS_reaction.split('+')]
		products = [i.strip() for i in RHS_reaction.split('+')]

		reaction = MODEL.createReaction()
		reaction.setId(reaction_name)
		reaction.setReversible(reversible)

		for i in reactants:
			stoichiometry = 1
			add_reac = i
			if re.search(r"\((.*?)\)",i):
				m = re.search(r"\((.*?)\)",i)
				stoichiometry = int(m.group(1))
				add_reac = re.sub(r'\(.*?\)', '', i).strip()
			add_species_by_id(add_reac,'metabolite',MODEL,add_reac, compartment)
			reactant_ref = reaction.createReactant()
			reactant_ref.setSpecies(add_reac)
			reactant_ref.setStoichiometry(stoichiometry)

		for i in products:
			stoichiometry = 1
			add_prod = i
			if re.search(r"\((.*?)\)",i):
				m = re.search(r"\((.*?)\)",i)
				stoichiometry = int(m.group(1))
				add_prod = re.sub(r'\(.*?\)', '', i).strip()
			add_species_by_id(add_prod,'metabolite',MODEL,add_prod, compartment)
			product_ref = reaction.createProduct()
			product_ref.setSpecies(add_prod)
			product_ref.setStoichiometry(stoichiometry)

		RNA_molecule_reactant = 'rna'+ molecule
		RNA_molecule_product = 'rna'+ molecule + "_" + str(int(row['Position']))

		add_species_by_id(RNA_molecule_reactant,'Rna', MODEL, RNA_molecule_reactant, compartment)
		add_species_by_id(RNA_molecule_product,'Rna', MODEL, RNA_molecule_product, compartment)	

		reactant_ref = reaction.createReactant()
		reactant_ref.setSpecies(RNA_molecule_reactant)
		product_ref = reaction.createProduct()
		product_ref.setSpecies(RNA_molecule_product) 

		mod_ref = reaction.createModifier()
		add_species_by_id(modifier,'Protein',MODEL,modifier_name, modifier_compartment)
		mod_ref.setSpecies(modifier)

libsbml.writeSBMLToFile( DOCUMENT, "output.xml")




	
	

	
	

	
























