directory = ''
file_name = 'LifeIsEasy.m'

variables_to_update = {
    'harmonicOne': '9',
    'harmonicTwo': '11',
    'harmonicThree': '10',
    'gaussians': '3',
    'harmonics' : '3',
    'reconstruct_single_color': 'false',
    'pulse_single_color': 'false',
    'overload': 'false',
    'overload_data': '"No data"',
    'chirp': 'true',
    'remove_frequency': 'false',
    'gaussian_variation': 'true',
    'first_color': '2',
    'windowing': 'false',
    'checkpoint': 'false',
    'gaussian_blur': 'true',
    'ultimate': 'false'
}

with open(directory + file_name, 'r') as file:
    lines = file.readlines()

modified_lines = []
for line in lines:
    modified_line = line
    for var, new_value in variables_to_update.items():
        if line.strip().startswith(var + ' ='):
            modified_line = f"{var} = {new_value};\n"
            break
    modified_lines.append(modified_line)

with open(directory + file_name, 'w') as file:
    file.writelines(modified_lines)

print("The variables have been updated in " + directory + file_name)
