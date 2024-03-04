def generate_dat_file(model_name, E, nu, gauss, nelx, nely, domain_size):
    # Calculate node coordinates
    nodes = []
    node_id = 1
    for i in range(nelx + 1):
        for j in range(nely + 1):
            x = domain_size[0] + (domain_size[2] - domain_size[0]) * i / nelx
            y = domain_size[1] + (domain_size[3] - domain_size[1]) * j / nely
            nodes.append((node_id, x, y, 0.0))
            node_id += 1

    # Calculate element node lists
    elements = []
    element_id = 1
    for j in range(nelx):
        for i in range(nely):
            n1 = i + j * (nely + 1) + 1
            n2 = n1 + 1
            n4 = n1 + nely + 1
            n3 = n4 + 1
            elements.append((element_id, n1, n2, n3, n4))
            element_id += 1
        

    # Write to dat file
    with open(f"{model_name}.dat", "w") as file:
        file.write(f"model: {model_name}\n")
        file.write(f"E: {E}\n")
        file.write(f"nu: {nu}\n")
        file.write(f"gauss: {gauss}\n")
        file.write(f"nodes: {len(nodes)}\n")
        file.write(f"elements: {len(elements)}\n")
        file.write("nodes_list:\n")
        for node in nodes:
            file.write(f"{node[0]} {node[1]} {node[2]} {node[3]}\n")

        file.write("elements_list: 2D\n")
        for element in elements:
            file.write(f"{element[0]} {element[1]} {element[2]} {element[3]} {element[4]}\n")

        file.write("boundary:\n\n")
        file.write("force:\n\n")
        file.write("end")

# 示例调用
generate_dat_file("mesh_model", 2.1e5, 0.3, 2, 30, 10, (0, 0, 30, 10))