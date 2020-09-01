from model.SBML import SBML
from model.SBMLtoODES import sbml_to_odes

model_files = [
    "interesting/subset_0011010000/model.sbml",
    "interesting/subset_0011010010/model.sbml",
    "interesting/subset_0111010110/model.sbml",
    "interesting/subset_0111011110/model.sbml",
]

thetas = [
    [0.00033527383068416645, 2.0728186060235227, 1532.8500211526236,
        0.019863774716908453, 0.0447076409567],
    [0.0007059017327560171, 1.2829049910728634, 1057.6723127046432,
        0.04202162319506584, 0.00047947160579723944, 5318.083763606885,
        0.06541012419846036],
    [0.0005971463010398104, 2.2221083234333405, 1969.8728480793766,
        0.03559772947312278, 0.00016273705204659835, 2971.0153561336683,
        0.04861901211204291, 827.7313066555923, 2.327589966701676,
        0.047233928724508545], 
    [0.0006341602666315066, 2.108990957627937, 1881.8274751530353,
        0.03763254177550618, 0.000201408804476085, 5336.187245874958,
        0.055875410909535586, 675.2688527098721, 2.937579145757671,
        0.5703339139502793, 28822.260108945276, 0.07466928068149971],
]

time_points = [30, 60, 90, 120, 150, 180, 210, 240]
expression = "RasGTP"

for model_idx, model_file in enumerate(model_files):
    sbml = SBML()
    sbml.load_file(model_file)
    odes = sbml_to_odes(sbml)

    parameters = odes.get_all_parameters()
    parameter_names = parameters.keys()
    sample = thetas[model_idx]
    for parameter_idx, parameter in enumerate(parameters.keys()):
        odes.define_parameter(parameter, sample[parameter_idx])

    integration = odes.evaluate_exp_on(expression, time_points)
    print(model_file, "\t", integration)

print("\n\n")

for model_idx, model_file in enumerate(model_files):
    sbml = SBML()
    sbml.load_file(model_file)
    odes = sbml_to_odes(sbml)

    parameters = odes.get_all_parameters()
    parameter_names = list(parameters.keys())
    sample = thetas[model_idx]
    print(model_file, "\t", parameter_names, "\t", sample)
