

## CH: 31.10.2018

class IHM_Chemical_descriptor():
	## TODO! THIS WILL BECOME PART OF THE IHM DICTIONARY. Thus, it should be removed from here if not needed anymore
	def __init__(self, auth_name=u"?",chem_comp_id=None,chemical_name=u"?",common_name=u"?",smiles=u"?",smiles_canonical=u"?",inchi=u"?",inchi_key=u"?"):
		self.auth_name = auth_name
		self.chem_comp_id = chem_comp_id
		self.chemical_name = chemical_name
		self.common_name = common_name
		self.smiles = smiles
		self.smiles_canonical = smiles_canonical
		self.inchi = inchi
		self.inchi_key = inchi_key


class Probe():
	""" 
	This class is not in the FLR dictionary, but it collects all the information connected by the probe_ids
	"""
	def __init__(self, probe_list_entry=None,probe_descriptor=None):
		self.probe_list_entry = probe_list_entry
		self.probe_descriptor = probe_descriptor
		
	def add_probe_list_entry(self,probe_list_entry):
		self.probe_list_entry = probe_list_entry
		
	def add_probe_descriptor(self,probe_descriptor):
		self.probe_descriptor = probe_descriptor
	
	def __str__(self):
		""" Example for testing"""
		return_string = u""
		return_string += u"Name: %s:"%(self.probe_list_entry.chromophore_name)
		return_string += u"\tSMILES (chromophore): %s"%(self.probe_descriptor.chromophore_chem_descriptor_id.smiles)
		return return_string
	
class Probe_descriptor():
	def __init__(self, reactive_probe_chem_descriptor_id=None, chromophore_chem_descriptor_id=None, chromophore_center_atom=None):
		self.reactive_probe_chem_descriptor_id = reactive_probe_chem_descriptor_id
		self.chromophore_chem_descriptor_id = chromophore_chem_descriptor_id
		self.chromophore_center_atom = chromophore_center_atom
		

class Probe_list():
	def __init__(self,chromophore_name=None,reactive_probe_flag=False, reactive_probe_name=None, probe_origin=None,probe_link_type=None):
		self.chromophore_name = chromophore_name
		self.reactive_probe_flag = reactive_probe_flag
		self.reactive_probe_name = reactive_probe_name
		self.probe_origin = probe_origin
		self.probe_link_type = probe_link_type
		

class Sample_probe_details_collection():
	""" 
	This class is not in the FLR dictionary.
	"""
	def __init__(self,sample_id=None,probe_id=None,fluorophore_type=None,description=None,poly_probe_position_id=None):
		self.sample_probe_details_list = []
		if sample_id != None and probe_id != None and fluorophore_type != None and poly_probe_position_id != None:
			self.add_sample_probe_combination(sample_id,probe_id,fluorophore_type,poly_probe_position_id,description)
		
	def add_sample_probe_combination(self,sample_id,probe_id,fluorophore_type,poly_probe_position_id,description=u""):
		self.sample_probe_details_list.append(Sample_probe_details(sample_id=sample_id,probe_id=probe_id,fluorophore_type=fluorophore_type,poly_probe_position_id=poly_probe_position_id,description=description))

	def __str__(self):
		return_string = u"Sample_probe_collection:\n"
		for entry in self.sample_probe_details_list:
			return_string += u"%s:\t"%(entry.sample_id.sample_description)
			return_string += u"%s\n"%(entry.description)
			
		return return_string
		
class Sample_probe_details():
	def __init__(self,sample_id=None,probe_id=None,fluorophore_type=None,description=None,poly_probe_position_id=None):
		self.sample_id = sample_id
		self.probe_id = probe_id
		self.fluorophore_type = fluorophore_type
		self.description = description
		self.poly_probe_position_id = poly_probe_position_id
		
	

class Poly_probe_conjugate():
	def __init__(self,sample_probe_id=None, chem_descriptor_id=None,ambiguous_stoichiometry=False,probe_stoichiometry=-1):
		self.sample_probe_id = sample_probe_id
		self.chem_descriptor_id = chem_descriptor_id
		self.ambiguous_stoichiometry = ambiguous_stoichiometry
		self.probe_stoichiometry = probe_stoichiometry

#class Poly_probe_position_collection():
#	"""
#	This class is not in the FLR dictionary.
#	"""
#	def __init__(self):
#		self.poly_probe_position_list = []
#
#	def add_poly_probe_position(self, entity_id=None, entity_description=u"",seq_id=-1,comp_id=u"",atom_id=u"",mutation_flag=False,modification_flag=False,auth_name=u"",chem_descriptor_mutated=None,chem_descriptor_modified=None):
#			if mutation_flag:
#				self.poly_probe_position_list.append(Poly_probe_position_mutated(entity_id=entity_id, \
#																			  entity_description=entity_description,\
#																			  seq_id=seq_id, \
#																			  comp_id=comp_id,\
#																			  atom_id=atom_id,\
#																			  mutation_flag=mutation_flag,\
#																			  modification_flag=modification_flag,\
#																			  auth_name=auth_name,\
#																			  chem_descriptor_id=chem_descriptor_mutated))
#			elif modification_flag:
#				self.poly_probe_position_list.append(Poly_probe_position_modified(entity_id=entity_id, \
#																			  entity_description=entity_description,\
#																			  seq_id=seq_id, \
#																			  comp_id=comp_id,\
#																			  atom_id=atom_id,\
#																			  mutation_flag=mutation_flag,\
#																			  modification_flag=modification_flag,\
#																			  auth_name=auth_name,\
#																			  chem_descriptor_id=chem_descriptor_modified))
#
#			else:
#				self.poly_probe_position_list.append(Poly_probe_position(entity_id=entity_id, \
#																			  entity_description=entity_description,\
#																			  seq_id=seq_id, \
#																			  comp_id=comp_id,\
#																			  atom_id=atom_id,\
#																			  mutation_flag=mutation_flag,\
#																			  modification_flag=modification_flag,\
#																			  auth_name=auth_name))
#


	
class Poly_probe_position():
	""" 
	This class combines Poly_probe_position, Poly_probe_position_modified, and Poly_probe_position_mutated from the FLR dictionary
	"""
	def __init__(self, entity_id=None, entity_description=u"",seq_id=-1,comp_id=u"",atom_id=u"",mutation_flag=False,modification_flag=False,auth_name=u"",mutated_chem_descriptor_id=None,modified_chem_descriptor_id=None):
		self.entity_id = entity_id
		self.entity_description = entity_description
		self.seq_id = seq_id
		self.comp_id = comp_id
		self.atom_id = atom_id
		self.mutation_flag = mutation_flag
		self.modification_flag = modification_flag
		self.auth_name = auth_name
		if self.mutation_flag:
			self.mutated_chem_descriptor_id = mutated_chem_descriptor_id
		if self.modification_flag:
			self.modified_chem_descriptor_id = modified_chem_descriptor_id
		

	def __str__(self):
		return_string = u"Poly_probe_position: "
		return_string += u"Name: %s"%(self.auth_name)
		return_string += u"\t%s%i\n"%(self.comp_id,self.seq_id)
		return return_string
		
#class Poly_probe_position_modified(Poly_probe_position):
#	def __init__(self,entity_id, entity_description,seq_id,comp_id,atom_id,mutation_flag,modification_flag,auth_name,chem_descriptor_id=None):
#		Poly_probe_position.__init__(entity_id=entity_id, entity_description=entity_description,seq_id=seq_id,comp_id=comp_id,atom_id=atom_id,mutation_flag=mutation_flag,modification_flag=modification_flag,auth_name=auth_name)
#		self.chem_descriptor_id = chem_descriptor_id
		
#class Poly_probe_position_mutated(Poly_probe_position):
#	def __init__(self,entity_id, entity_description,seq_id,comp_id,atom_id,mutation_flag,modification_flag,auth_name,chem_descriptor_id=None):
#		Poly_probe_position.__init__(entity_id=entity_id, entity_description=entity_description,seq_id=seq_id,comp_id=comp_id,atom_id=atom_id,mutation_flag=mutation_flag,modification_flag=modification_flag,auth_name=auth_name)
#		self.chem_descriptor_id = chem_descriptor_id



class Sample():
	def __init__(self, entity_assembly=None, num_of_probes=0,sample_condition=None, sample_description=u"",sample_details=u"",solvent_phase=u""):
		self.entity_assembly = entity_assembly
		self.num_of_probes = num_of_probes
		self.sample_condition = sample_condition
		self.sample_description= sample_description
		self.sample_details = sample_details
		self.solvent_phase = solvent_phase

	def __str__(self):
		return_string = u"Sample:\n"
		return_string += u"\t%s"%(self.entity_assembly)
		return_string += u"\t%s"%(self.sample_condition)
		return_string += u"\tNumber of probes:%i\n"%(self.num_of_probes)
		return_string += u"\tSample description:%s\n"%(self.sample_description)
		return_string += u"\tSample details:%s\n"%(self.sample_details)
		return_string += u"\tSolvent phase: %s\n"%(self.solvent_phase)
		return return_string

class Entity_assembly():
    def __init__(self,entity_id=None,num_copies=0,entity_description=u""):
        self.entity_list = []
        self.num_copies_list = []
        self.entity_description_list = []
        if entity_id != None and num_copies != 0:
            self.add_entity(entity_id,num_copies,entity_description)        
        
    def add_entity(self,entity_id,num_copies,entity_description=u""):
        if num_copies < 0:
            print(u"Error in Entity_assembly: Number of copies for Entity must be larger than zero.")
        self.entity_list.append(entity_id)
        self.num_copies_list.append(num_copies)
        self.entity_description_list.append(entity_description)

    def __str__(self):
        return_string = u"Entity_assembly:\n"
        for i in range(len(self.entity_list)):
            return_string += u"\t%ix %s:\t\"%s\"\n"%(self.num_copies_list[i],self.entity_list[i],self.entity_description_list[i])
        return return_string
            

class Sample_condition():
	def __init__(self,details=u""):
		self.details = details

	def __str__(self):
		return_string = u"Sample_condition:\n"
		return_string += u"\t%s"%(self.details)
		return return_string   

class Experiment():
	def __init__(self,instrument=None,exp_setting=None,sample=None,details=u"?"):
		self.instrument_list = []
		self.exp_setting_list = []
		self.sample_list = []
		self.details_list = []
		if instrument != None and exp_setting != None and sample != None:
			self.add_entry(instrument=instrument,exp_setting=exp_setting,sample=sample,details=details)
		
	def add_entry(self,instrument,exp_setting,sample,details=u"?"):
		self.instrument_list.append(instrument)
		self.exp_setting_list.append(exp_setting)
		self.sample_list.append(sample)
		self.details_list.append(details)
		
	def __str__(self):
		return_string = "Experiment:\n"
		for i in range(len(self.sample_list)):
			return_string += "\tSample: %s, Instrument: %s, Experimental setting: %s\n"%(self.sample_list[i],self.instrument_list[i],self.exp_setting_list[i])
		return return_string
		

		

class Instrument():
	def __init__(self, details = u""):
#		self.id = None
		self.details = details
	
	def __str__(self):
		return_string = u"%s"%(self.details)
		return return_string

class Exp_setting():
	def __init__(self, details = u""):
		self.details = details
	def __str__(self):			
		return u"%s"%(self.details)
		
####
class Fret_analysis():
	def __init__(self,experiment_id=None,sample_probe_id_1=None,sample_probe_id_2=None,forster_radius_id=None,calibration_parameters_id=None,method_name=None,chi_square_reduced=None,data_list_id=None,external_file_id=None,software_id=None):
		self.experiment_id = experiment_id
		self.sample_probe_id_1 = sample_probe_id_1
		self.sample_probe_id_2 = sample_probe_id_2
		self.forster_radius_id = forster_radius_id
		self.calibration_parameters_id = calibration_parameters_id
		self.method_name = method_name
		self.chi_square_reduced = chi_square_reduced
		self.dataset_list_id = data_list_id
		self.external_file_id = external_file_id
		self.software_id = software_id
        
class Fret_distance_restraint():
	def __init__(self,distance_restraint_id=None,group_id=None,sample_probe_id_1=None,sample_probe_id_2=None,state_id=None,analysis_id=None,distance=-1,distance_error_plus=0,distance_error_minus=0,distance_type=None,population_fraction=0,peak_assignment_id=None):
		self.distance_restraint_id = distance_restraint_id
		self.group_id = group_id
		self.sample_probe_id_1 = sample_probe_id_1
		self.sample_probe_id_2 = sample_probe_id_2
		self.state_id = state_id
		self.analysis_id = analysis_id
		self.distance = distance
		self.distance_error_plus = distance_error_plus
		self.distance_error_minus = distance_error_minus
		self.distance_type = distance_type
		self.population_fraction = population_fraction
		self.peak_assignment_id = peak_assignment_id
        
class Fret_forster_radius():
	def __init__(self,donor_probe_id=None,acceptor_probe_id=None,forster_radius=None,reduced_forster_radius=None):
		self.donor_probe_id = donor_probe_id
		self.acceptor_probe_id = acceptor_probe_id
		self.forster_radius = forster_radius
		self.reduced_forster_radius = reduced_forster_radius
        
class Fret_calibration_parameters():
	def __init__(self,phi_acceptor=None,alpha=None,alpha_sd=None,gG_gR_ratio=None,beta=None,gamma=None,delta=None,a_b=None):
		self.phi_acceptor = phi_acceptor
		self.alpha = alpha
		self.alpha_sd = alpha_sd
		self.gG_gR_ratio = gG_gR_ratio
		self.beta = beta
		self.gamma = gamma
		self.delta = delta
		self.a_b = a_b
        
class Peak_assignment():
	def __init__(self,method_name=None,details=None):
		self.method_name = method_name
		self.details = details

class Fret_model_quality():
	def __init__(self,chi_squared_reduced=None,dataset_group_id=None,method=None,details=None):
		self.chi_squared_reduced = chi_squared_reduced
		self.dataset_group_id = dataset_group_id
		self.method = method
		self.details = details
        

########## Example for testing
## Creating the Entity_assembly
Entity_Assembly1 = Entity_assembly()
Entity_Assembly1.add_entity(entity_id=u"1",num_copies=1,entity_description=u"RNA (37-MER) - alpha strand")
Entity_Assembly1.add_entity(entity_id=u"2",num_copies=1,entity_description=u"RNA (37-MER) - beta strand")
Entity_Assembly1.add_entity(entity_id=u"3",num_copies=1,entity_description=u"RNA (32-MER) - gamma strand")
Entity_Assembly1.add_entity(entity_id=u"4",num_copies=1,entity_description=u"RNA (32-MER) - delta strand")
## Creating the sample conditions
Sample_condition1 = Sample_condition(details=u"Description of sample conditions for (D)beta5c_(A)alpha12d")
Sample_condition2 = Sample_condition(details=u"Description of sample conditions for (D)beta5c_(A)gamma8b")
## Creating the samples
Sample1 = Sample(entity_assembly=Entity_Assembly1,num_of_probes=2,sample_condition=Sample_condition1,sample_description=u"(D)beta5c_(A)alpha12d",sample_details=u"?",solvent_phase=u"liquid")
Sample2 = Sample(entity_assembly=Entity_Assembly1,num_of_probes=2,sample_condition=Sample_condition2,sample_description=u"(D)beta5c_(A)gamma8b",sample_details=u"?",solvent_phase=u"liquid")

#print Sample1.entity_assembly
#print Sample1.sample_condition

print Sample1
print Sample2

## Experiment
Instrument1 = Instrument(u"Instrument")
Exp_settings1 = Exp_setting(u"Experimental settigns")
Experiment1 = Experiment()
Experiment1.add_entry(instrument=Instrument1,exp_setting=Exp_settings1,sample=Sample1)
Experiment1.add_entry(instrument=Instrument1,exp_setting=Exp_settings1,sample=Sample2)

print Experiment1.sample_list[0].sample_description
print Experiment1.sample_list[1].sample_description


## For Sample_probe_details, we require Sample, Probe, and Probe position
## Probe - Alexa488
Probe_Alexa488 = Probe()
Probe_list_Alexa488 = Probe_list(chromophore_name=u"Alexa488",reactive_probe_flag=True,reactive_probe_name=u"Alexa 488-tetrafluor-phenoxy ester triethylammonium",probe_origin=u"extrinsic",probe_link_type=u"covalent")
Probe_descriptor_reactive_Alexa488 = IHM_Chemical_descriptor(auth_name=u"Alexa488-tetrafluor-phenoxy ester triethylammonium",chem_comp_id=u"?",chemical_name=u"?",common_name=u"?",smiles=u"C1=CC(=C(C=C1C(=O)OC2=C(C(=CC(=C2F)F)F)F)C(=O)O[H])C4=C3C=CC(C(=C3OC5=C4C=CC(=C5[S](O[H])(=O)=O)N([H])[H])[S](O[H])(=O)=O)=N[H].CCN(CC)CC")
Probe_descriptor_chromophore_Alexa488 = IHM_Chemical_descriptor(auth_name=u"Alexa488 - free acid",common_name=u"Alexa Fluor 488",smiles=u"C1=CC(=C(C=C1C(=O)O[H])C(=O)O[H])C3=C2C=CC(=[N+]([H])[H])C(=C2OC4=C3C=CC(=C4[S](=O)(=O)O[H])N([H])[H])[S](=O)(=O)O[H]")
Probe_descriptor_Alexa488 = Probe_descriptor(reactive_probe_chem_descriptor_id=Probe_descriptor_reactive_Alexa488,chromophore_chem_descriptor_id=Probe_descriptor_chromophore_Alexa488,chromophore_center_atom=u"?")
Probe_Alexa488.add_probe_descriptor(Probe_descriptor_Alexa488)
Probe_Alexa488.add_probe_list_entry(Probe_list_Alexa488)
## Probe - Cy5 (different way of adding information)
Probe_list_Cy5 = Probe_list(chromophore_name=u"Cy5",reactive_probe_flag=True,reactive_probe_name=u"Sulfo-Cy5-N-hydroxy succinimidylester",probe_origin=u"extrinsic",probe_link_type=u"covalent")
Probe_descriptor_Cy5 = Probe_descriptor(reactive_probe_chem_descriptor_id=IHM_Chemical_descriptor(auth_name=u"Sulfo-Cy5-N-hydroxy succinimidylester",smiles=u"CCN2C1=C(C=C(C=C1)[S](=O)(=O)[O-])C(C2=CC=CC=CC3=[N+](C4=C(C3(C)C)C=C(C=C4)[S](=O)(=O)O)CCCCCC(=O)ON5C(CCC5=O)=O)(C)C")\
											  ,chromophore_chem_descriptor_id=IHM_Chemical_descriptor(auth_name=u"Cy5 - free acid",common_name=u"Cy5",smiles=u"CCN1C2=C(C=C(C=C2)S(=O)(=O)[O-])C(C1=CC=CC=CC3=[N+](C4=C(C3(C)C)C=C(C=C4)S(=O)(=O)O)CCCCCC(=O)O)(C)C")\
											  ,chromophore_center_atom=u"?")
Probe_Cy5 = Probe(probe_list_entry=Probe_list_Cy5,probe_descriptor=Probe_descriptor_Cy5)

print Probe_Alexa488
print Probe_Cy5

## Poly_probe_positions
Modified_chem_descriptor_Uridine=IHM_Chemical_descriptor(auth_name=u"Uridine-C6-Amino Linker",common_name=u"?",smiles=u"C1(=CN(C(=O)NC1=O)C2C(C(C(O2)CO[P](=O)(O)O)O)O)C=CC(=O)NCCCCCCN")
## As I do not have entities defined here, for testing, I simply refer to the entry of the entity assembly
## beta5c
Poly_probe_position_beta5c = Poly_probe_position(entity_id = Entity_Assembly1.entity_list[1],entity_description=u"RNA (37-MER) - beta strand",seq_id=5,comp_id=u"U",atom_id=u"C5",mutation_flag=False,modification_flag=True,auth_name="(D)beta5c",modified_chem_descriptor_id=Modified_chem_descriptor_Uridine)
## alpha12d
Poly_probe_position_alpha12d = Poly_probe_position(entity_id = Entity_Assembly1.entity_list[0],entity_description=u"RNA (37-MER) - alpha strand",seq_id=12,comp_id=u"U",atom_id=u"C5",mutation_flag=False,modification_flag=True,auth_name="(D)alpha12d",modified_chem_descriptor_id=Modified_chem_descriptor_Uridine)
## gamma8b
Poly_probe_position_gamma8b = Poly_probe_position(entity_id = Entity_Assembly1.entity_list[2],entity_description=u"RNA (32-MER) - gamma strand",seq_id=8,comp_id=u"U",atom_id=u"C5",mutation_flag=False,modification_flag=True,auth_name="(D)gamma8b",modified_chem_descriptor_id=Modified_chem_descriptor_Uridine)

Sample_probe_details_collection1 = Sample_probe_details_collection()
Sample_probe_details_collection1.add_sample_probe_combination(sample_id=Sample1,probe_id=Probe_Alexa488,fluorophore_type=u"donor",poly_probe_position_id=Poly_probe_position_beta5c,description=u"Donor in (D)beta5c_(A)alpha12d")
Sample_probe_details_collection1.add_sample_probe_combination(sample_id=Sample1,probe_id=Probe_Cy5,fluorophore_type=u"acceptor",poly_probe_position_id=Poly_probe_position_alpha12d,description=u"Acceptor in (D)beta5c_(A)alpha12d")
Sample_probe_details_collection1.add_sample_probe_combination(sample_id=Sample2,probe_id=Probe_Alexa488,fluorophore_type=u"donor",poly_probe_position_id=Poly_probe_position_beta5c,description=u"Donor in (D)beta5c_(A)gamma8b")
Sample_probe_details_collection1.add_sample_probe_combination(sample_id=Sample2,probe_id=Probe_Cy5,fluorophore_type=u"acceptor",poly_probe_position_id=Poly_probe_position_gamma8b,description=u"Acceptor in (D)beta5c_(A)gamma8b")
print Sample_probe_details_collection1

##  Poly_probe_conjugate

## Fret_analysis

## FRET_forster_radius

## FRET_calibration_parameters

## FRET_distance_restraint
## Peak_assignment

## FRET_model_quality

## FRET_model_distance
