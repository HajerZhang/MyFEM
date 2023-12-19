import os

def read_inp_file(file_path):
    nodes = {}
    elements = []
    boundary_conditions = []
    forces = []

    with open(file_path, 'r') as file:
        reading_nodes = False
        reading_elements = False
        reading_boundary_conditions = False
        reading_forces = False

        for line in file:
            line = line.strip()

            if line.startswith('*NODE'):
                reading_nodes = True
                continue
            elif reading_nodes and line.startswith('*'):
                reading_nodes = False

            if reading_nodes:
                data = line.split(',')
                node_id = int(data[0].strip())
                x, y, z = float(data[1].strip()), float(data[2].strip()), float(data[3].strip())
                nodes[node_id] = (x, y, z)

            if line.startswith('*ELEMENT'):
                reading_elements = True
                continue
            elif reading_elements and line.startswith('*'):
                reading_elements = False

            if reading_elements:
                data = line.split(',')
                element_id = int(data[0].strip())
                node_ids = [int(node.strip()) for node in data[1:]]
                elements.append((element_id, node_ids))

            if line.startswith('*BOUNDARY'):
                reading_boundary_conditions = True
                continue
            elif reading_boundary_conditions and line.startswith('*'):
                reading_boundary_conditions = False

            if reading_boundary_conditions and line:
                boundary_conditions.append(line)

            if line.startswith('*CLOAD'):
                reading_forces = True
                continue
            elif reading_forces and line.startswith('*'):
                reading_forces = False

            if reading_forces and line:
                forces.append(line)

    return nodes, elements, boundary_conditions, forces

def write_dat_file(model_name, nodes, elements, boundary_conditions, forces, output_path):
    with open(output_path, 'w') as dat_file:
        # Write model name
        dat_file.write(f'model: {model_name}\n')

        # Write nodes and elements count
        dat_file.write(f'nodes: {len(nodes)}\n')
        dat_file.write(f'elements: {len(elements)}\n')

        # Write nodes list
        dat_file.write('nodes_list:\n')
        for node_id, (x, y, z) in nodes.items():
            dat_file.write(f'{node_id} {x} {y} {z}\n')

        # Write element list
        dat_file.write('element_list:\n')
        for element_id, node_ids in elements:
            dat_file.write(f'{element_id} {" ".join(map(str, node_ids))}\n')

        # Write boundary conditions
        dat_file.write('\nboundary:\n')
        dat_file.write('\n'.join(boundary_conditions))
        dat_file.write('\n')

        # Write forces
        dat_file.write('force:\n')
        dat_file.write('\n'.join(forces))
        dat_file.write('\n')

        # Write end
        dat_file.write('end\n')

if __name__ == "__main__":
    inp_file_path = "clip.inp"
    dat_output_path = "preprocess.dat"
    model_name = os.path.splitext(os.path.basename(inp_file_path))[0]

    nodes, elements, boundary_conditions, forces = read_inp_file(inp_file_path)
    write_dat_file(model_name, nodes, elements, boundary_conditions, forces, dat_output_path)
