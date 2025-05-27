'''
MaXLinker GUI - A GUI for the MaXLinker module

This script contains:
 (1) outputMaXLinkerCombinedResults - A wrapper function to combine multiple MaXLinker results
 (2) generateMaXLinker - A wrapper function to call MaXLinker function
 (3) runMaXLinker - The routine called when users press Run
 (4) __main__ The main routine to call the MaXLinker GUI
'''

__author__ = "Alden Leung"
__version__ = "1.0"

import os

import tkinter as tk
from tkinter import messagebox

from resources import *
from utilities import *

from MaXLinker1pt0_Submission import *

def outputMaXLinkerCombinedResults(filename, fraction_CSM_dicts, FDR_cutoff):
	'''
	A wrapper to process, combine and output all complete MaXLinker results
	'''
	# The suggested main function routine to combine all results according to Yugandhar
	tot_CSM_dict=defaultdict()
	for fraction_CSM_dict in fraction_CSM_dicts:
		tot_CSM_dict.update(fraction_CSM_dict)
	tot_uniq_XL=XL_recurrence_score(tot_CSM_dict)
	FDR_tot_XL=filter_XL_FDR(tot_uniq_XL, FDR_cutoff)
	write_XL(FDR_tot_XL, filename)

def generateMaXLinker(config_dict):	
	'''
	A process worker to obtain MaXLinker results
	'''
	# The parameter names are used here to partially match the MaxLinker main script functions, and hence, a bit messy parsing routine
	
	uniprot_dict = config_dict["fasta_dict"]
	FDR_cutoff = config_dict["FDR(%)"]
	MS2_rescue_mode = "NO" if config_dict["MS2_rescue_mode"] == 0 else "YES" # The main function accepts the string "YES" as True, instead of accepting a boolean
	file_prefix = config_dict["file_prefix"]
	
	# MS2_DeisoDeconv is input but not used?
	MS1_file_path = config_dict["MS1"]
	MS2_file_path = config_dict["MS2"]
	MS3_file_path = config_dict["MS3"]
	msn_1_MS3_path = config_dict["MS3_MSn-1"]
	PD_result_path = config_dict["MS3_PSM"]
	MS2_PD_result_path = config_dict["MS2_Rescue_PSM"]
	result_path = config_dict["Output_file_suffix"]	
	linker_mass = float(config_dict["Linker_mass"])	
	linker_Long = float(config_dict["Linker_long_mass"])
	linker_Short = float(config_dict["Linker_short_mass"])
	linker_dm = linker_Long - linker_Short
	linker_deficit = linker_mass - linker_Long - linker_Short
	linker_Long_name_SEQUEST = config_dict["Linker_long_name"]
	linker_Short_name_SEQUEST = config_dict["Linker_short_name"]
	
	

	prec_tol_ppm=20 # Hard-coded parameter
	
	# The suggested main function routine for each MS file according to Yugandhar
	maxlinkresult = find_XL_from_mgf(MS2_rescue_mode, MS1_file_path, MS2_file_path, MS3_file_path, msn_1_MS3_path, PD_result_path, MS2_PD_result_path, uniprot_dict, linker_mass, linker_Long, linker_Short, linker_dm, linker_deficit, linker_Long_name_SEQUEST, linker_Short_name_SEQUEST)
	fraction_CSM_dict=filter_pep_pair(maxlinkresult, mass_dict, mod_dict, linker_mass, prec_tol_ppm, file_prefix)
	fraction_uniq_XL=XL_recurrence_score(fraction_CSM_dict)	
	FDR_fraction_XL=filter_XL_FDR(fraction_uniq_XL, FDR_cutoff)		
	write_XL(FDR_fraction_XL, result_path) # output here
	
	return fraction_CSM_dict

def runMaXLinker(app, parameter_dict):
	'''
	Parse the parameters and deal with paths and file suffix issues. 
	Call batchRunGUI on the selected file prefices
	'''	
	config_dict, error_parameters = obtainParameters(parameter_dict, required_parameters = None)
	
	if len(error_parameters) == 0:
		
		inputfilepath = config_dict["Input_path"]		
		required_input_file_suffices_paramnames = [parameter.name for parameter in parameter_dict.values() if parameter.paramType == "inputsuffix"]
		
		outputfilepath = config_dict["Output_path"]
		required_output_file_suffices_paramnames = [parameter.name for parameter in parameter_dict.values() if parameter.paramType == "outputsuffix"]
		
		# Auto search for files first
		required_input_file_suffices = [config_dict[paramnames] for paramnames in required_input_file_suffices_paramnames]
		path_sets = walkForFilesGUI(app, inputfilepath, required_input_file_suffices)
		if path_sets is not None:
			# Handle the options first
			config_dict["fasta_dict"] = read_unip_fasta(config_dict.pop("Database"))			
			combinedresultfile = os.path.join(outputfilepath, config_dict.pop("Batch_output_file")) # Currently we have no such suffices here
			
			FDR_cutoff = config_dict["FDR(%)"]
			cpuno = int(config_dict.pop("No_of_threads"))
			
			# Customized config dicts for suffix
			customized_config_dicts = []			
			for file_prefix, path_set in path_sets.items():
				customized_config_dict = dict(config_dict)
				customized_config_dict["file_prefix"] = file_prefix
				for i in range(0, len(required_input_file_suffices_paramnames)):
					customized_config_dict[required_input_file_suffices_paramnames[i]] = path_set[i]
				for i in range(0, len(required_output_file_suffices_paramnames)):
					customized_config_dict[required_output_file_suffices_paramnames[i]] = os.path.join(outputfilepath, file_prefix + config_dict[required_output_file_suffices_paramnames[i]])
				customized_config_dicts.append(customized_config_dict)
			
			# generateMaXLinker(customized_config_dict)
			batchRunGUI(app, generateMaXLinker, customized_config_dicts, lambda maxlinkresults,filename=combinedresultfile: outputMaXLinkerCombinedResults(filename, maxlinkresults, FDR_cutoff), cpuno=cpuno)
		else:
			pass
		
	else:
		tk.messagebox.showinfo("MaXLinker", "The following parameters are improper: \n" + "\n".join(error_parameters))


#### MAIN ####
if __name__ == '__main__': 

	parameterSet = [(Category("1", "Mass spectrometry file input"), [Parameter("Input_path","","loadpath","."),
											Parameter("MS1", "", "inputsuffix", ".dta.zip"), 
											Parameter("MS2", "", "inputsuffix", "_[Node_07].mgf"), Parameter("MS3", "", "inputsuffix", "_[Node_09].mgf"), Parameter("MS3_MSn-1", "", "inputsuffix", "_MS3.mgf"), Parameter("MS2_DeisoDeconv", "", "inputsuffix", "_[Node_05].mgf"), 
											Parameter("MS3_PSM", "", "inputsuffix", "_MS3.tsv"), Parameter("MS2_Rescue_PSM", "", "inputsuffix", "_MS2_rescue.tsv")]), 											
					(Category("2", "Database file input"), [Parameter("Database", "", "loadfastafile", "database/TargetRandom_Human_uniprot_taxid_9606_20170623.fasta")]),
					(Category("3", "Output"), [Parameter("Output_path","","loadpath","."), Parameter("Output_file_suffix", "Result output prefix", "outputsuffix", "_out.tsv"), Parameter("Batch_output_file", "", "outputfile", "uniq.tsv")]),
					(Category("4", "linker parameters"), [Parameter("Linker_mass", "", "float", 158.003766), Parameter("Linker_long_mass", "", "float", 85.982637), Parameter("Linker_short_mass", "", "float", 54.010565), Parameter("Linker_long_name", "", "string", "DSSO_L"), Parameter("Linker_short_name", "", "string", "DSSO_S")]),
					(Category("5", "Scoring parameters"), [Parameter("FDR(%)", "", "float", "1.0")]), 
					(Category("6", "Module mode selection"), [Parameter("MS2_rescue_mode", "", "modulemode", "YES")]),
					(Category("7", "Multithreading"), [Parameter("No_of_threads", "Maximum no. of threads", "int", 1)])
					]	
	parameter_dict = {parameter.name:parameter for (category, parameters) in parameterSet for parameter in parameters}

	# GUI begins here
	PROGRAM_TITLE = "MaXLinker GUI"
	LABEL_WIDTH = 20
	ENTRY_WIDTH = 50	
	CHOOSE_FILE_PATH_BUTTON_WIDTH = 10	
	
	app = tk.Tk()
	app.title(PROGRAM_TITLE)
	iconimage = tk.PhotoImage(file=icon_path)
	app.tk.call("wm", "iconphoto", app._w, iconimage)
	
	app_width = 600
	app_height = 825
	app.geometry('%dx%d+%d+%d' % (app_width, app_height, (app.winfo_screenwidth() - app_width) / 2, (app.winfo_screenheight() - app_height) / 2))

	canvas=tk.Canvas(app)
	frame=tk.Frame(canvas)
	vertscrollbar=tk.Scrollbar(app,orient="vertical",command=canvas.yview)
	canvas.configure(yscrollcommand=vertscrollbar.set)
	vertscrollbar.pack(side="right",fill="y")
	canvas.pack(side="left", fill="both", expand=True)
	canvas.create_window((1,1), window=frame, anchor="nw")
	def onFrameConfigure(canvas):
		canvas.configure(scrollregion=canvas.bbox("all"))
	frame.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure(canvas))
	
	rowno = 0
	logoimage = tk.PhotoImage(file=logo_path)
	pimagelabel = tk.Label(frame, image=logoimage)
	pimagelabel.image = logoimage
	pimagelabel.grid(row=rowno)
	rowno += 1
	
	# Initialize main parameters frame	
	initializeParameterMenu(frame, parameterSet, rowno, LABEL_WIDTH, ENTRY_WIDTH, CHOOSE_FILE_PATH_BUTTON_WIDTH)	
	rowno += len(parameterSet)
	
	# Initialize command panel	
	runlabelframe = tk.LabelFrame(frame, text="Actions")
	runlabelframe.grid(row = rowno, padx=25, pady=5, ipadx=5, ipady=5)
	loadbutton = tk.Button(runlabelframe, text="Load configuration...", command = lambda:loadConfigFilePrompt(parameter_dict, [("MaXLinker Config File", "*.mcfg")]))
	loadbutton.grid(row = 0, column = 0, padx=10)
	savebutton = tk.Button(runlabelframe, text="Save configuration...", command = lambda:saveConfigFilePrompt(parameterSet, ".mcfg", [("MaXLinker Config File", "*.mcfg")]))
	savebutton.grid(row = 0, column = 1, padx=10)
	runbutton = tk.Button(runlabelframe, text="Run...", command = lambda pd=parameter_dict:runMaXLinker(app, pd))	
	runbutton.grid(row = 0, column = 2, padx=10)
	app.bind('<Return>', lambda e: runMaXLinker(app, parameter_dict))
	app.mainloop()
