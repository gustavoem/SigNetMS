from model.SBML import SBML
from model.SBMLtoODES import sbml_to_odes

model_file = "input/fail/model.sbml"
sbml = SBML()
sbml.load_file(model_file)
odes = sbml_to_odes(sbml)
time_points = [30, 60, 90, 120, 150, 180, 210, 240]
expression = "RasGTP"

theta = [1.59, 12000]
parameters = odes.get_all_parameters()
for idx, parameter in enumerate(parameters.keys()):
    odes.define_parameter(parameter, theta[idx])

integration = odes.evaluate_exp_on(expression, time_points)
print(integration)
