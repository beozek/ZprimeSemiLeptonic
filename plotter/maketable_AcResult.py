import json
import os

def extract_ac_values(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    
    ac_values = {}
    for poi in data['POIs']:
        if poi['name'] == 'Ac':
            ac_values['measured_ac'] = poi['fit'][1] / 100
            ac_values['measured_ac_err_up'] = (poi['fit'][2] - poi['fit'][1]) / 100
            ac_values['measured_ac_err_down'] = (poi['fit'][1] - poi['fit'][0]) / 100
            break
    return ac_values

json_files = {
    "[0,500]": "/nfs/dust/cms/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/plotter/Combined_years/impacts_0_500.json",
    "[500,750]": "/nfs/dust/cms/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/plotter/Combined_years/impacts_500_750.json",
    "[750,1000]": "/nfs/dust/cms/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/plotter/Combined_years/impacts_750_1000.json",
    "[1000,1500]": "/nfs/dust/cms/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/plotter/Combined_years/impacts_1000_1500.json",
    "[1500,Inf]": "/nfs/dust/cms/user/beozek/uuh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/plotter/Combined_years/impacts_1500Inf.json"
}

sm_ac_values = {
    "[0,500]": {"sm_ac": 0.56, "sm_ac_err_up": 0.06, "sm_ac_err_down": 0.06},
    "[500,750]": {"sm_ac": 0.72, "sm_ac_err_up": 0.06, "sm_ac_err_down": 0.06},
    "[750,1000]": {"sm_ac": 0.79, "sm_ac_err_up": 0.04, "sm_ac_err_down": 0.06},
    "[1000,1500]": {"sm_ac": 0.96, "sm_ac_err_up": 0.09, "sm_ac_err_down": 0.09},
    "[1500,Inf]": {"sm_ac": 0.94, "sm_ac_err_up": 0.11, "sm_ac_err_down": 0.15}
}

table = ''
table += '\\begin{tabular}{|c|c|c|} \n'
table += '\\hline\n'
table += 'Mass Bin $(M_{tt})$ & Measured $A_C$ & SM Predicted $A_C$ \\\\ \n'
table += '\\hline\n'

for mass_bin, json_file in json_files.items():
    ac_values = extract_ac_values(json_file)
    
    measured_ac = ac_values['measured_ac']
    measured_ac_err_up = ac_values['measured_ac_err_up']
    measured_ac_err_down = ac_values['measured_ac_err_down']
    
    sm_ac = sm_ac_values[mass_bin]['sm_ac']
    sm_ac_err_up = sm_ac_values[mass_bin]['sm_ac_err_up']
    sm_ac_err_down = sm_ac_values[mass_bin]['sm_ac_err_down']
    
    table += '%s GeV & $%.4f^{+%.4f}_{-%.4f}$ & $%.2f^{+%.2f}_{-%.2f}$ \\\\ \n' % (
        mass_bin, measured_ac, measured_ac_err_up, measured_ac_err_down,
        sm_ac, sm_ac_err_up, sm_ac_err_down)

table += '\\hline\n'
table += '\\end{tabular}'

print(table)

with open('output_table.tex', 'w') as file:
    file.write(table)
