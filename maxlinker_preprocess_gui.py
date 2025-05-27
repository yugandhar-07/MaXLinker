'''
MaXLinker Preprocessing GUI - A GUI for the MaXLinker Preprocessing module

This script contains:
 (1) generateRescueFile - A wrapper function to call MaXLinker preprocessing function
 (2) runPreProcessing - The routine called when users press Preprocess
 (3) __main__ The main routine to call the MaXLinker Preprocessing GUI
'''

__author__ = "Alden Leung"
__version__ = "1.0"


import tkinter as tk
from tkinter import messagebox

from resources import *
from utilities import *

from MaXLinker1pt0_Submission import *


#### MAIN ####
def generateRescueFile(config_dict):
	'''
	A process worker to perform pre-processing
	'''
	# The parameter names are used here to partially match the MaxLinker main script functions, and hence, a bit messy parsing routine
	MS2_DeisoDeconv_file_path = config_dict["MS2_DeisoDeconv"]
	MS2_file_path = config_dict["MS2"]
	MS3_file_path = config_dict["MS3"]
	msn_1_MS3_path = config_dict["MS3_MSn-1"]
	rescue_file_path = config_dict["Output_file_suffix"]
	linker_mass = float(config_dict["Linker_mass"])	
	linker_Long = float(config_dict["Linker_long_mass"])
	linker_Short = float(config_dict["Linker_short_mass"])
	linker_dm = linker_Long - linker_Short
	linker_deficit = linker_mass - linker_Long - linker_Short
	MaXLinker_Generate_Rescue_PrecMismatch_MS2mgf_Using_mgf(MS2_DeisoDeconv_file_path, MS2_file_path, MS3_file_path, msn_1_MS3_path, rescue_file_path, linker_dm, linker_deficit)


def runPreProcessing(app, param_dict):
	'''
	Parse the parameters and deal with paths and file suffix issues. 
	Call batchRunGUI on the selected file prefices
	'''
	config_dict, error_parameters = obtainParameters(param_dict, required_parameters = None)
	
	if len(error_parameters) == 0:
		inputfilepath = config_dict["Input_path"]
		required_input_file_suffices_paramnames = [parameter.name for parameter in parameter_dict.values() if parameter.paramType == "inputsuffix"]
		
		outputfilepath = config_dict["Output_path"]
		required_output_file_suffices_paramnames = [parameter.name for parameter in parameter_dict.values() if parameter.paramType == "outputsuffix"]
		
		required_input_file_suffices = [config_dict[paramnames] for paramnames in required_input_file_suffices_paramnames]
		path_sets = walkForFilesGUI(app, inputfilepath, required_input_file_suffices)
		if path_sets is not None:
			cpuno = int(config_dict.pop("No_of_threads"))

			customized_config_dicts = []
			for file_prefix, path_set in path_sets.items():
				customized_config_dict = dict(config_dict)
				for i in range(0, len(required_input_file_suffices_paramnames)):
					customized_config_dict[required_input_file_suffices_paramnames[i]] = path_set[i]
				for i in range(0, len(required_output_file_suffices_paramnames)):
					customized_config_dict[required_output_file_suffices_paramnames[i]] = os.path.join(outputfilepath, file_prefix + config_dict[required_output_file_suffices_paramnames[i]])
				customized_config_dicts.append(customized_config_dict)
			
			# generateRescueFile(customized_config_dict)
			batchRunGUI(app, generateRescueFile, customized_config_dicts, None, cpuno)
		else:
			# Cancel the run
			pass
		
	else:
		tk.messagebox.showinfo("Parameter Errors", "The following parameters are improper: \n" + "\n".join(error_parameters))

if __name__ == '__main__': 
	
	parameterSet = [(Category("1", "Mass spectrometry file input"), [Parameter("Input_path","","loadpath","."),  
											Parameter("MS2", "", "inputsuffix", "_[Node_07].mgf"), Parameter("MS3", "", "inputsuffix", "_[Node_09].mgf"), Parameter("MS3_MSn-1", "", "inputsuffix", "_MS3.mgf"), Parameter("MS2_DeisoDeconv", "", "inputsuffix", "_[Node_05].mgf")]),
					(Category("2", "Output"), [Parameter("Output_path","","loadpath","."), Parameter("Output_file_suffix", "", "outputsuffix", "_MS2_rescue.mgf")]),
					(Category("3", "Linker parameters"), [Parameter("Linker_mass", "", "float", 158.003766), Parameter("Linker_long_mass", "", "float", 85.982637), Parameter("Linker_short_mass", "", "float", 54.010565)]),
					(Category("4", "Multithreading"), [Parameter("No_of_threads", "Maximum no. of threads", "int", 1)])
					]
	parameter_dict = {parameter.name:parameter for (category, parameters) in parameterSet for parameter in parameters}
	
	# GUI begins here
	PROGRAM_TITLE = "MaXLinker Preprocessing GUI"
	LABEL_WIDTH = 20
	ENTRY_WIDTH = 50	
	CHOOSE_FILE_PATH_BUTTON_WIDTH = 10
	
	app = tk.Tk()
	app.title(PROGRAM_TITLE)
	iconimage = tk.PhotoImage(file=icon_path)
	app.tk.call("wm", "iconphoto", app._w, iconimage)

	app_width = 600
	app_height = 525
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
	loadbutton = tk.Button(runlabelframe, text="Load configuration...", command =  lambda:loadConfigFilePrompt(parameter_dict, [("MaXLinker Preprocessing Config File", "*.pcfg")]))
	loadbutton.grid(row = 0, column = 0, padx=10)
	savebutton = tk.Button(runlabelframe, text="Save configuration...", command = lambda:saveConfigFilePrompt(parameterSet, ".pcfg", [("MaXLinker Preprocessing Config File", "*.pcfg")]))	
	savebutton.grid(row = 0, column = 1, padx=10)
	preprocessbutton = tk.Button(runlabelframe, text="Preprocess...", command = lambda pd=parameter_dict:runPreProcessing(app, pd))	
	preprocessbutton.grid(row = 0, column = 2, padx=10)
	app.bind('<Return>', lambda e: runPreProcessing(app, parameter_dict))
	app.mainloop()
