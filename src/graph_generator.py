from networkx.generators.random_graphs import erdos_renyi_graph

def generator(node_num, p):
    g = erdos_renyi_graph(node_num, p)
    # [0, 1, 2, 3, 4, 5]
    print('nodes: ', g.nodes)
    print('# of nodes: ', g.number_of_nodes())
    # [(0, 1), (0, 2), (0, 4), (1, 2), (1, 5), (3, 4), (4, 5)]
    print('edges: ', g.edges)
    print('# of edges: ', g.number_of_edges())
    # Write .in file
    print('Input file is generating...')
    file_name = '../input/' + str(g.number_of_nodes()) + '_' + str(g.number_of_edges()) + '.in'
    with open(file_name, 'w') as f:
        # num of vertices
        f.write(str(g.number_of_nodes()) + '\n')
        # num of edges
        f.write(str(g.number_of_edges()) + '\n')
        # each edge
        for n1, n2, attr in g.edges(data=True):
            f.write("%d %d\n" % (n1, n2))
    print('Finish generating input file...')
    return

if __name__ == '__main__':
    while True:
        node_num = input('Please enter the number of nodes you want to generate (any positive integer): ')
        if node_num.isdigit() and eval(node_num) > 0:
            break
    while True:
        p = input('Please enter the random possibility (between 0 and 1): ')
        if (0 <= eval(p) <= 1):
            break
    generator(eval(node_num), eval(p))


# Example code
# n = 6
# p = 0.5
# g = erdos_renyi_graph(n, p)
# print(g.nodes)
# # [0, 1, 2, 3, 4, 5]
# print(g.edges)
# # [(0, 1), (0, 2), (0, 4), (1, 2), (1, 5), (3, 4), (4, 5)]