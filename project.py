import tkinter as tk
from tkinter import messagebox
import sympy as sp
import copy
from collections import deque


class GraphGUI:
    def __init__(self, master):
        self.master = master
        master.title("Graph GUI")

        self.canvas = tk.Canvas(master, width=400, height=300, bg="white")
        self.canvas.pack(pady=20)

        self.node_btn = tk.Button(master, text="Add Node", command=self.add_node_mode)
        self.node_btn.pack(side=tk.LEFT, padx=(20, 10))

        self.edge_btn = tk.Button(master, text="Add Edge", command=self.add_edge_mode)
        self.edge_btn.pack(side=tk.LEFT)

        self.clear_btn = tk.Button(master, text="Clear", command=self.clear_canvas)
        self.clear_btn.pack(side=tk.RIGHT, padx=(10, 20))
        
        self.calc_poly_btn = tk.Button(master, text="Calculate Chromatic Polynomial", command=self.display_chromatic_polynomial)
        self.calc_poly_btn.pack(side=tk.LEFT, padx=(10, 10))

        self.check_graph_type_btn = tk.Button(master, text="Check Graph Type", command=self.check_graph_type)
        self.check_graph_type_btn.pack(side=tk.LEFT, padx=(10, 10))

        self.mode = "node"  # Start in node-adding mode
        self.nodes = {}  # Store nodes with their IDs
        self.node_counter = 1
        self.edges = []
        self.temp_edge_start = None

        self.canvas.bind("<Button-1>", self.canvas_click)

    def add_node_mode(self):
        self.mode = "node"

    def add_edge_mode(self):
        self.mode = "edge"
        self.temp_edge_start = None

    def clear_canvas(self):
        self.canvas.delete("all")
        self.nodes.clear()
        self.edges.clear()
        self.node_counter = 1
        self.temp_edge_start = None
        
    def check_graph_type(self):
        adjacency_list = self.to_adjacency_list()
        
        # Check for a graph with only isolated vertices (no edges)
        if not self.edges:
            messagebox.showinfo("Graph Type", "The graph consists only of isolated vertices.")
            return
        
        if self.is_complete_graph(adjacency_list):
            messagebox.showinfo("Graph Type", "The graph is complete.")
        elif self.is_tree(adjacency_list):
            messagebox.showinfo("Graph Type", "The graph is a tree.")
        elif self.is_cycle_graph(adjacency_list):
            messagebox.showinfo("Graph Type", "The graph has a cycle.")
        elif self.is_bipartite(adjacency_list):
            messagebox.showinfo("Graph Type", "The graph is bipartite.")
        else:
            messagebox.showinfo("Graph Type", "The graph type is not specifically identified.")


    def canvas_click(self, event):
        if self.mode == "node":
            node_id = self.canvas.create_oval(event.x-10, event.y-10, event.x+10, event.y+10, fill="blue")
            self.canvas.create_text(event.x, event.y, text=str(self.node_counter), fill="white")
            self.nodes[self.node_counter] = (event.x, event.y, node_id)
            self.node_counter += 1
        elif self.mode == "edge":
            nearest_node = self.find_nearest_node(event.x, event.y)
            if nearest_node and self.temp_edge_start is None:
                self.temp_edge_start = nearest_node
            elif nearest_node and self.temp_edge_start:
                start_x, start_y, _ = self.nodes[self.temp_edge_start]
                end_x, end_y, _ = self.nodes[nearest_node]
            if self.temp_edge_start == nearest_node:  # If the start and end nodes are the same
               self.canvas.create_oval(start_x, start_y, start_x-20, start_y-40, outline="red")  # Draw a larger, red loop
               self.edges.append((self.temp_edge_start, self.temp_edge_start))  # Add a loop to the edges
               adjacency_list = self.to_adjacency_list()
               adjacency_list[nearest_node].append(nearest_node)  # Add the node to its own list of neighbors
            
            else:
                self.canvas.create_line(start_x, start_y, end_x, end_y, fill="black")
                self.edges.append((self.temp_edge_start, nearest_node))
                self.temp_edge_start = None

    def find_nearest_node(self, x, y):
        nearest_node = None
        min_distance = float("inf")
        for node_id, (node_x, node_y, _) in self.nodes.items():
            distance = (node_x - x)**2 + (node_y - y)**2
            if distance < min_distance:
                min_distance = distance
                nearest_node = node_id
        if min_distance <= 100:  # Arbitrary threshold for "closeness"
            return nearest_node
        return None

    def to_adjacency_list(self):
        adjacency_list = {node: [] for node in self.nodes}
        for node1, node2 in self.edges:
            adjacency_list[node1].append(node2)
            adjacency_list[node2].append(node1)
        return adjacency_list
    
    def chromatic_polynomial_generic(self, adjacency_list):
        x = sp.symbols('x')  # Symbolic variable for the polynomial
        
        # Check for an empty graph or a single node without edges
        if not any(adjacency_list.values()):
            return x**len(adjacency_list)
        
        # Check for a loop, which makes chromatic number 0
        for node, neighbors in adjacency_list.items():
            if node in neighbors:
                return 0
        
        # Choose an edge
        node1, node2 = next((node, neighbors[0]) for node, neighbors in adjacency_list.items() if neighbors)
        
        # Copy of the graph with the edge deleted
        adjacency_list_minus_e = copy.deepcopy(adjacency_list)
        if node2 in adjacency_list_minus_e[node1]:
            adjacency_list_minus_e[node1].remove(node2)
        if node1 in adjacency_list_minus_e[node2]:
            adjacency_list_minus_e[node2].remove(node1)
        
        # Copy of the graph with the edge contracted
        adjacency_list_over_e = copy.deepcopy(adjacency_list)
        
        # Merge node2's neighbors into node1 and remove node2
        neighbors2 = adjacency_list_over_e.pop(node2)
        adjacency_list_over_e[node1] = [n for n in (adjacency_list_over_e[node1] + neighbors2) if n != node1 and n != node2]
        
        # Replace node2 with node1 in other nodes' neighbor lists and remove duplicates
        for node, neighbors in adjacency_list_over_e.items():
            adjacency_list_over_e[node] = [n if n != node2 else node1 for n in neighbors]
            adjacency_list_over_e[node] = list(dict.fromkeys(adjacency_list_over_e[node]))  # Remove duplicates

        # Recursive call for the graph with the edge deleted
        P_minus_e = self.chromatic_polynomial_generic(adjacency_list_minus_e)

        # Recursive call for the graph with the edge contracted
        P_over_e = self.chromatic_polynomial_generic(adjacency_list_over_e)

        return P_minus_e - P_over_e

    
    def display_chromatic_polynomial(self):
        adjacency_list = self.to_adjacency_list()
        polynomial = self.compute_chromatic_polynomial(adjacency_list)
        messagebox.showinfo("Chromatic Polynomial", f"The chromatic polynomial is:\n{polynomial}")
    
    def compute_chromatic_polynomial(self, adjacency_list):
        x = sp.symbols('x')
        if not self.edges:  # Check if there are no edges
            # Chromatic polynomial for a graph with only isolated vertices is x^n
            return x**len(self.nodes)
        if self.is_complete_graph(adjacency_list):
            # Chromatic polynomial for complete graph K_n is x(x-1)...(x-(n-1))
            return sp.prod([x - i for i in range(len(adjacency_list))])
        elif self.is_tree(adjacency_list) or self.is_bipartite(adjacency_list):
            # Chromatic polynomial for a tree or a bipartite graph is x(x-1)^(n-1)
            return x * (x - 1)**(len(adjacency_list) - 1)
        elif self.is_cycle_graph(adjacency_list):
            # Chromatic polynomial for cycle C_n is (x-1)^n + (-1)^n * (x-1)
            n = len(adjacency_list)
            return (x - 1)**n + (-1)**n * (x - 1)
        else:
            # If the graph type is not recognized, use a generic method (deletion-contraction)
            return self.chromatic_polynomial_generic(adjacency_list)
        
    def is_complete_graph(self, adjacency_list):
        # A complete graph should have n*(n-1)/2 edges
        expected_num_edges = len(self.nodes) * (len(self.nodes) - 1) // 2
        if len(self.edges) != expected_num_edges:
            return False
        
        # Check that each node is connected to every other node
        for node, neighbors in adjacency_list.items():
            if len(neighbors) != len(self.nodes) - 1:
                return False
        
        return True


    def is_cycle_graph(self, adjacency_list):
        # Each node in a cycle graph has exactly 2 edges and the number of nodes equals the number of edges
        if len(self.edges) != len(adjacency_list):
            return False
        for edges in adjacency_list.values():
            if len(edges) != 2:
                return False
        return True
    

    def is_bipartite(self, adjacency_list):
        color = {}
        for node in adjacency_list:
            if node not in color:
                queue = deque([node])
                color[node] = 0  # Start coloring with color 0
                while queue:
                    current = queue.popleft()
                    for neighbour in adjacency_list[current]:
                        if neighbour not in color:
                            color[neighbour] = 1 - color[current]  # Assign an alternate color
                            queue.append(neighbour)
                        elif color[neighbour] == color[current]:
                            return False  # Same color on both ends of an edge means it's not bipartite
        return True



    def is_tree(self, adjacency_list):
        # A tree must have exactly n-1 edges where n is the number of nodes
        if len(self.edges) != len(adjacency_list) - 1:
            return False
        
        # Start DFS from the first node to check for connectivity and cycles
        visited = set()
        if self.has_cycle_dfs(next(iter(adjacency_list)), visited, adjacency_list, None):
            return False  # If a cycle is found, it's not a tree

        # If not all nodes were visited, the graph is not connected, hence not a tree
        return len(visited) == len(adjacency_list)

    def has_cycle_dfs(self, node, visited, adjacency_list, parent):
        visited.add(node)
        for neighbour in adjacency_list[node]:
            if neighbour not in visited:
                if self.has_cycle_dfs(neighbour, visited, adjacency_list, node):
                    return True
            elif parent is not None and neighbour != parent:
                # If an adjacent node is visited and not parent of the current node, there's a cycle
                return True
        return False

if __name__ == "__main__":
    root = tk.Tk()
    gui = GraphGUI(root)
    root.mainloop()