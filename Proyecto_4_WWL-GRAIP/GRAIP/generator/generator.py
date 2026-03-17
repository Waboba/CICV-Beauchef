import numpy as np
import matplotlib.pyplot as plt

from generator.graphlet_counts import three_counts, four_counts, five_counts
from generator.sampling import sample
from generator.models import GRAIP, SwapCon, BA_graph



class GraphGenerator(object):
    
    def __init__(self, G=None, graph_data=None, max_gl_size=5):
        """
        Construct a graph generator object.

        Parameters
        ----------
        G : NetworkX graph, optional
            The probabilistic target graph. Every edge should have an attribute 'probability'. If a path for graph data is specified, this is ignored.
        graph_data : str, optional
            A path leading to an .npz file produced by the sample function. At least one of G and graph_data should be given.
        max_gl_size : int, optional
            The maximum order of graphlets taken into account. This can be 3, 4 or 5. The default is 5.
        """
        
        if graph_data is not None:
            self.data = graph_data
        elif G is not None:
            for edge in G.edges():
                try:
                    G.edges[edge]["probability"]
                except:
                    raise TypeError("Not all edges of the probabilistic graph G have attribute 'probability'. This is required for sampling.")
            self.probabilistic_graph = G
        if G is None and graph_data is None:
            raise ValueError("Provide either a probabilistic graph G or a path where the graph data from a previously sampled network is stored.")
        
        if max_gl_size == 3:
            self.graphlets = ['wedge', 'triangle']
            self.valid = {'wedge': {3,5,6}, 'triangle': {7}}
            self.count_func = three_counts
        elif max_gl_size == 4:
            self.graphlets = ['wedge', 'triangle', '4star', '4path', 'tailed_tri', '4cycle', 'diamond', '4clique']
            self.valid = {'wedge': {3,5,6}, 'triangle': {7}, '4star': {56,11,21,38}, '4path': {35,37,41,44,13,14,49,50,19,22,26,28},
                          'tailed_tri': {39,43,46,15,53,54,23,57,58,27,60,29}, '4cycle': {51,45,30}, 'diamond': {47,55,59,61,62,31},
                          '4clique': {63}}
            self.count_func = four_counts
        elif max_gl_size == 5:
            self.graphlets = ['wedge', 'triangle', '4star', '4path', 'tailed_tri', '4cycle', 'diamond', '4clique', '5star', 'prong',
                              '5path', 'fork_tailed_tri', 'long_tailed_tri', 'double_tailed_tri', 'tailed_cycle', '5cycle', 'hourglass',
                              'cobra', 'stingray', 'hatted_cycle', 'three_wedge', 'three_tri', 'tailed_clique', 'triangle_strip',
                              'diamond_wedge', 'wheel', 'hatted_clique', 'bipyramid', '5clique']
            self.valid = {'wedge': {3,5,6}, 'triangle': {7},
                          '4star': {56,11,21,38}, '4path': {35,37,41,44,13,14,49,50,19,22,26,28}, 'tailed_tri': {39,43,46,15,53,54,23,57,58,27,60,29},
                          '4cycle': {51,45,30}, 'diamond': {47,55,59,61,62,31}, '4clique': {63},
                          '5star': {960,294,75,149,568},
                          'prong': {897,898,645,774,904,267,523,139,141,270,401,147,533,277,278,792,150,538,540,156,418,291,165,550,166,680,553,293,556,300,177,306,562,561,312,184,833,706,579,708,836,456,329,202,77,78,464,848,83,212,85,90,480,736,99,356,102,105,624,120},
                          '5path': {771,643,773,646,777,393,650,652,780,269,142,394,785,402,786,531,526,534,275,408,154,282,284,417,673,674,547,549,163,424,297,169,172,305,178,581,582,204,332,337,210,594,340,596,86,525,344,92,353,609,226,228,612,101,232,108,240,368,113,114},
                          'fork_tailed_tri': {405,661,151,157,422,295,806,302,181,310,696,569,570,824,572,961,962,964,968,331,203,587,79,976,213,91,992,358,107,632},
                          'long_tailed_tri': {647,775,908,397,398,527,914,659,662,535,409,793,283,666,285,412,929,803,805,551,809,426,171,682,428,174,817,690,313,186,583,716,589,590,844,460,466,339,850,722,342,345,346,604,481,737,227,865,229,233,234,620,241,370,244,117,118,372,628,124},
                          'double_tailed_tri': {901,902,143,271,913,406,279,920,539,155,541,668,930,421,167,936,299,555,812,558,689,818,309,565,182,566,185,314,188,316,707,835,709,838,457,841,458,714,205,334,465,211,468,724,87,856,602,93,482,355,484,868,103,744,617,110,752,880,121,122},
                          'tailed_cycle': {899,905,906,395,651,653,782,779,403,789,790,794,796,542,286,158,419,677,678,681,684,173,557,301,433,434,179,307,563,440,837,710,333,206,849,595,852,597,341,214,472,218,220,94,738,611,376,357,230,614,488,361,740,364,109,496,625,626,115,248},
                          '5cycle': {675,613,425,236,781,654,369,242,787,598,410,348},
                          'hourglass': {993,807,235,972,430,591,978,245,374,663,825,698,347,636,413},
                          'cobra': {903,399,918,921,667,924,543,287,933,938,811,940,559,175,945,946,821,694,567,315,187,317,190,444,839,711,717,846,462,461,467,723,470,343,857,473,858,732,605,730,483,867,485,231,745,746,490,876,873,622,753,882,884,756,500,119,249,378,125,126},
                          'stingray': {917,407,437,438,669,159,934,423,814,303,693,822,183,311,697,826,571,828,189,318,573,574,952,700,963,965,966,969,970,459,715,843,335,207,977,980,469,725,215,984,603,219,221,95,994,996,486,870,1000,359,619,363,366,111,1008,888,760,633,634,123},
                          'hatted_cycle': {909,910,655,783,915,791,922,411,795,797,414,670,931,679,937,427,683,429,813,686,691,819,441,442,497,498,845,718,851,726,599,854,474,476,349,606,350,860,739,869,741,615,489,492,237,238,621,748,881,754,243,371,373,246,629,380,377,250,630,252},
                          'three_wedge': {742,907,365,685,798,627,435,853,504,222},
                          'three_tri': {998,971,367,981,439,1016,635,701,830,223},
                          'tailed_clique': {949,950,956,319,575,191,967,463,471,985,731,733,487,1002,875,878,1012,761,890,127},
                          'triangle_strip': {637,892,638,919,925,671,415,935,942,815,431,695,823,953,954,827,699,829,446,702,445,973,974,847,719,979,982,727,986,247,988,477,859,375,607,351,475,995,997,871,1001,491,1004,747,494,623,239,1009,1010,251,501,757,502,886,889,762,379,764,253,382},
                          'diamond_wedge': {911,923,926,799,939,941,687,947,443,855,861,862,478,734,743,749,750,493,877,499,755,883,758,631,885,505,506,508,381,254},
                          'wheel': {509,510,759,1005,766,943,751,1011,507,863,887,955,893,990,927},
                          'hatted_clique': {383,1020,639,951,957,958,447,703,831,975,983,987,989,735,479,999,1003,1006,495,879,1013,1014,503,1017,1018,891,763,765,894,255},
                          'bipyramid': {511,1007,767,1015,895,959,1019,1021,1022,991},
                          '5clique': {1023}}
            self.count_func = five_counts
        else:
            raise ValueError("max_gl_size can only be 3, 4 or 5. Larger graphlet sizes not implemented.")
    
    
    def get_properties(self, N=None, savepath=None, plot=False):
        """
        If a graph data path was specified upon generator object construction, read the graph properties from the .npz file.
        Otherwise, a number of samples N has to be provided and sampling is performed.

        Parameters
        ----------
        N : int, optional
            Number of samples. The default is None.
        savepath : str, optional
            If sampling is performed, the results are saved to the given path in .npz file format. The default is None.
        plot : bool, optional
            If True, plot the degree distribution and graphlet frequencies. Plots are saved under 'degree.PNG' and 'graphlets.PNG'. The default is False.
        """
        
        if hasattr(self, 'data'):
            print("Reading graph properties from %s" %self.data)
            with np.load(self.data) as npzfile:
                self.E_n = npzfile['E_nodes']
                self.std_n = npzfile['std_nodes']
                self.E_e = npzfile['E_edges']
                self.std_e = npzfile['std_edges']
                self.E_deg = npzfile['E_degrees']
                self.std_deg = npzfile['std_degrees']
                self.E_gl = npzfile['E_graphlets'][:len(self.graphlets)]
                self.std_gl = npzfile['std_graphlets'][:len(self.graphlets)]
                self.samples = npzfile['samples']
        
        else:
            if N is None:
                raise ValueError("Sampling required. Please provide a maximum number of samples N and (optionally) a path savepath to save the graph properties to.")
            self.E_n, self.std_n, self.E_e, self.std_e, self.E_gl, self.std_gl, self.E_deg, self.std_deg = sample(self.probabilistic_graph, N, len(self.graphlets), self.count_func, savepath=savepath)
            self.samples = N
        
        
        if plot:
            # Degree distribution
            
            P_target = np.trim_zeros(self.E_deg, 'b')/self.E_n
            P_bounds = 2*np.trim_zeros(self.std_deg, 'b')/self.E_n
            low = P_target - P_bounds
            up = P_target + P_bounds
            deg = np.arange(len(P_target))
            
            plt.figure(0)
            plt.scatter(deg, P_target, s=5, zorder=2, color='tab:blue')
            plt.gca().fill_between(deg, up, low, alpha=0.25, color='tab:blue', zorder=1)
            plt.title("Degree distribution")
            plt.xlabel("degree $k$")
            plt.ylabel("$P(k)$")
            plt.xscale('log')
            plt.yscale('log')
            plt.savefig("degree.PNG", bbox_inches='tight')
            
            # Graphlet counts
            
            plt.figure(1)
            plt.errorbar(self.graphlets, self.E_gl, yerr=2*self.std_gl, fmt='o', capsize=3, ms=2, elinewidth=1, zorder=1, color='tab:blue')
            plt.title("Graphlet counts")
            plt.yscale('log')
            plt.ylabel("counts")
            plt.xticks(rotation=90, fontsize=7)
            plt.savefig("graphlets.PNG", bbox_inches='tight')
    
    
    def generate(self, model='GRAIP', graphs=1, savepath=None, max_steps=None, node_step=5, w=2/3, max_rej=None, temperature=0.01, cooling=0.99, threshold=0.05, max_reject=None):
        """
        Generate graphs according to the properties of the target graph.

        Parameters
        ----------
        model : str, optional
            The generator model. Options are 'GRAIP', 'SwapCon' and 'BA'. The default is 'GRAIP'.
        graphs : int, optional
            The number of graphs to be generated. The default is 1.
        savepath : str, optional
            Generated graphs are saved in edgelist format under 'savepath/[model]_graph[n].txt', n = 0, .., (graphs-1). The default is None.
        
        
        Depending on the model, the following extra parameters can be specified:
        
        GRAIP:
            max_steps : int, optional
                Maximum number of iterations. If None, the maximum is set to 100 times the mean number of edges.
            node_step : int, optional
                A node is added or removed every node_step steps. The default is 5.
            w : float, optional
                Weight factor for cost function. The default is 2/3.
            max_rej : int, optional
                Number of rejected steps in a row before a step is automatically accepted. If None, this is set to 2% of the mean number of edges.
        
        SwapCon:
            temperature : float, optional
                Initial temperature for the simulated annealing algorithm. The default is 0.01.
            cooling : float, optional
                Cooling factor for the simulated annealing algorithm.. The default is 0.99.
            threshold : float, optional
                Energy threshold for convergence of the simulated annealing algorithm. The default is 0.05.
            max_reject : int, optional
                Alternative convergence criterion: the algorithm is stopped after max_reject swaps in a row have been rejected. If None, this is set to the mean number of edges.
        
        BA:
            No extra parameters needed.

        Returns
        -------
        generated_graphs : list
            List of generated graphs.
        """
        
        if not hasattr(self, 'E_n'):
            if not hasattr(self, 'data'):
                raise NameError("Graph properties not yet determined and no data path was specified upon generator object construction. Please perform sampling first using object.get_properties(N, savepath=path)")
            else:
                self.get_properties()
        
        properties = [self.E_n, self.std_n, self.E_e, self.std_e, self.E_gl, self.std_gl, self.E_deg, self.std_deg]
        graphlet_config = [self.graphlets, self.valid, self.count_func]
        
        generated_graphs = []
        
        for i in range(graphs):
            
            full_path = savepath + "/" + model + "_graph%s.txt" %i
        
            if model == "GRAIP":
                H = GRAIP(properties, self.samples, graphlet_config, max_steps=max_steps, node_step=node_step, w=w, max_rej=max_rej, savepath=full_path)
        
            elif model == "SwapCon":
                H = SwapCon(properties, self.count_func, temperature=temperature, cooling=cooling, threshold=threshold, max_reject=max_reject, savepath=full_path)
        
            elif model == "BA":
                H = BA_graph(self.E_n, self.std_n, self.E_e, savepath=full_path)
        
            else:
                raise ValueError("%s is not a valid model. Valid options are 'GRAIP', 'SwapCon' and 'BA'." %model)
            
            generated_graphs.append(H)
            
            print("Graph %s/%s generated." %(i+1, graphs))
        
        return generated_graphs
            
