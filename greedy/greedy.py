import sys
import time
import networkx as nx


# obtain furthest reachable from v with distance d

# obtain furthest reachable from v with distance d
def furthest_reachable_position(G, v, d, end_of_backbone):
    try:
        E_reachable = nx.bfs_edges(G, v, depth_limit=d)
        V_reachable_backbone = [v for e in E_reachable for v in e if int(v) <= end_of_backbone]
        return max(V_reachable_backbone)
    except nx.NetworkXError:  # Handle NetworkX error when node is not in the graph
        return v



# def furthest_reachable_position(G, v, d, end_of_backbone):
#     E_reachable = nx.bfs_edges(G, v, depth_limit=d)
#     V_reachable_backbone = [v for e in E_reachable for v in e if int(v) <= end_of_backbone]
#     return max(V_reachable_backbone)


def load_data(edges_file_name, location_substring_file_name, location_cost_file_name):

    # load graph edges
    edge_file = open(edges_file_name, 'r')
    E = []
    end_of_ref = -1
    for e in edge_file:
        e = e.strip().split()
        E.append((e[0], e[1], e[2], e[3]))

        if e[3] == '-' and int(e[1]) > end_of_ref:
            end_of_ref = int(e[1])

    edge_file.close()

    # load variant locations and costs
    location_cost_file = open(location_cost_file_name, 'r')
    locations = []
    cost_bounds = []
    for l in location_cost_file:
        arr = l.strip().split()
        locations.append(arr[0])
        cost_bounds.append(int(arr[1]))
    location_cost_file.close()

    # load number of variants at each location
    location_substring_file = open(location_substring_file_name, 'r')
    number_of_variants_at = [0]*len(locations)
    l_prev = ''

    i = 0
    for l in location_substring_file:
        l = l.strip()
        if len(l) > 0:
            if l == l_prev:
                number_of_variants_at[i] += 1
            if l != l_prev:
                number_of_variants_at[i] = 1
                i += 1
            l_prev = l

    return E, locations, cost_bounds, end_of_ref, number_of_variants_at


def construct_graph(E):
    G = nx.MultiDiGraph()
    for e in E:
        # start, end, symbol, variant #)
        G.add_edge(e[0], e[1], symbol=e[2], variant=e[3])
    return G


def greedy(G, locations, costs, alpha, delta, end_of_backbone):

    deleted = [0]*len(locations)

    for i in range(len(locations)-1, -1, -1):

        print('\rProcessing location :', i)

        furthest_reachable = furthest_reachable_position(G, locations[i], alpha, end_of_backbone)
        j = i+1
        deletion_sum = 0
        while j < len(locations) and int(locations[j]) <= int(furthest_reachable):
            if deleted[j] == 1:
                deletion_sum += costs[j]
            j += 1

        if deletion_sum + costs[i] <= delta:
            deleted[i] = 1
    return deleted


if __name__ == "__main__":

    edges_file_name = sys.argv[1]
    location_cost_file_name = sys.argv[2]
    location_substring_file_name = sys.argv[3]
    alpha = int(sys.argv[4])
    delta = int(sys.argv[5])

    print('Loading data...')
    E, locations, costs, end_of_backbone, number_of_variants_at = load_data(edges_file_name,
                                                                            location_substring_file_name,
                                                                            location_cost_file_name)

    print('Constructing graph...')
    G = construct_graph(E)

    print('Running Greedy...')
    start = time.time()
    deleted = greedy(G, locations, costs, alpha, delta, end_of_backbone)
    end = time.time()
    total_time = end - start
    print("Total time:", total_time)

    number_variants_deleted = sum([deleted[i] * number_of_variants_at[i] for i in range(len(deleted))])

    print('\nNumber of variants removed: ' + str(number_variants_deleted))
    print('Number variants retained: ' + str(sum(number_of_variants_at) - number_variants_deleted))

