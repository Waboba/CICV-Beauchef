import os
import networkx as nx
import yaml

from generator.generator import GraphGenerator


if __name__ == '__main__':
    
    with open('config.yaml', 'r') as file:
        config = yaml.safe_load(file)
    
    # Check if folder for graphs exists
    
    savepath = config["generator"]["savepath"]
    
    if savepath is not None:
        if not os.path.exists(savepath):
            os.makedirs(savepath)
    
    
    graph_path = config["parameters"]["G"]
    if graph_path is not None:
        G = nx.read_gml(graph_path, label=None)
    else:
        G = None
    
    gen = GraphGenerator(G=G, graph_data=config["parameters"]["graph_data"], max_gl_size=config["parameters"]["n_g"])
    
    gen.get_properties(N=config["sampling"]["S"], savepath=config["sampling"]["savefile"], plot=config["sampling"]["plot"])
    
    gen.generate(model=config["generator"]["model"], graphs=config["generator"]["graphs"], savepath=config["generator"]["savepath"],
                 max_steps=config["GRAIP"]["max_s"], node_step=config["GRAIP"]["node_step"], w=config["GRAIP"]["w"], max_rej=config["GRAIP"]["max_rej"],
                 temperature=config["SwapCon"]["temperature"], cooling=config["SwapCon"]["cooling"], threshold=config["SwapCon"]["threshold"], max_reject=config["SwapCon"]["max_reject"])