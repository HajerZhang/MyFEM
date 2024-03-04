import matplotlib.pyplot as plt
from tqdm import tqdm

def read_dat_file(file_path):
    nodes = {}
    elements = []

    read_nodes = False
    read_elements = False

    with open(file_path, 'r') as file:
        lines = file.readlines()

        for line in lines:
            if line.startswith("nodes_list"):
                read_nodes = True
                read_elements = False
                continue
            elif line.startswith("element_list"):
                read_nodes = False
                read_elements = True
                continue
            elif line.startswith("boundary"):
                break

            if read_nodes:
                data = line.split()
                node_id = int(data[0].strip())
                x, y, z = float(data[1].strip()), float(data[2].strip()), float(data[3].strip())
                nodes[node_id] = (x, y, z)
            elif read_elements:
                data = line.split()
                element_id = int(data[0].strip())
                node_ids = [int(node.strip()) for node in data[1:]]
                elements.append((element_id, node_ids))

    return nodes, elements

def plot_elements(nodes, elements, model_name, save_path, dpi=600):
    fig, ax = plt.subplots()

    for element in tqdm(elements, desc='Plotting Elements', unit='element'):
        x_coords = [nodes[node_id][0] for node_id in element[1]]
        y_coords = [nodes[node_id][1] for node_id in element[1]]
        ax.plot(x_coords + [x_coords[0]], y_coords + [y_coords[0]], color='black', linestyle='-', linewidth=0.5)

    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')

    ax.text(1, 1.04, f'Model: {model_name}  |  Nodes: {len(nodes)}  |  Elements: {len(elements)}', transform=ax.transAxes, fontsize=6, verticalalignment='top', horizontalalignment='right')
    ax.set_title('Model Presentation', loc='center', pad=10)
    plt.savefig(save_path, bbox_inches='tight', dpi=dpi)

if __name__ == "__main__":
    dat_file_path = "./test.dat"
    save_path = "./mesh_model.png"
    model_name = "triangle"
    nodes, elements = read_dat_file(dat_file_path)
    plot_elements(nodes, elements, model_name, save_path, dpi=300)
