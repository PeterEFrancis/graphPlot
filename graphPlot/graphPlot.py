import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
import matplotlib.animation as animation
import random as r
from typing import *



class Node(object):
    def __init__(self, ID: int = -1, connections: list = None, pos: np.array = None):
        """
        Initialize Node object.

        Args:
            ID - node identifier
            connections - list of connected nodes
            pos - (x,y) position of node
        """
        self.ID = ID
        self.connections = [] if connections is None else connections
        self.degree = 0 if connections is None else len(connections)
        self.pos = np.array([0.0, 0.0]) if pos is None else pos
        self.pastPos = np.array([[0.0,0.0],[0.0,0.0]])


    def set_coord(self, pos: np.array):
        """
        Set node position

        Args:
            pos - (x,y) position of node
        """
        self.pos = pos

    def add_connection(self, connection: "Node"):
        """
        Add connection to node and update the node degree

        Args:
            connection - connected node
        """
        self.connections += [connection]
        self.degree += 1

    def add_connections(self, connections: list):
        """
        Add connection to node and update the node degree

        Args:
            connection - list of connected nodes
        """
        self.connections += connections
        self.degree += len(connections)

    def distance_to(self, node: "Node"):
        """
        Find the distance to node

        Args:
            node - Node object

        Return:
            Distance to node
        """
        return np.sqrt((self.pos[0] - node.pos[0]) ** 2 + (self.pos[1] - node.pos[1]) ** 2)

    # how would you check the connections being equal?
    def __eq__(self, node: "Node"):
        """
        Return if equal to node

        Args:
            node - Node object

        Return:
            boolean truth value of equality
        """
        boolID = self.ID == node.ID
        boolPos = all(self.pos == node.pos)
        boolDegree = self.degree == node.degree
        return boolID and boolPos and boolDegree

    def __repr__(self):
        """
        retrurn string representation of (self) node
        """
        return f"<Node {self.ID} @ ({round(self.pos[0],3)},{round(self.pos[1],3)})>"



class SpringBoard(object):
    def __init__(self, nodesDict: Dict[int, List[int]], k: float, Q: float, nodePosDict: Dict[int, Tuple[float, float]] = {}):
        """
        Construct a SpringBoard object.
        If the nodes are all centered at the orgin, spread them out.

        Args:
            nodesDict - an adjacency dictionary of first positive integers
            k - coefficient of spring
            Q - coefficient of electric field
            g - gravitational coefficient
            nodePosDict - (optional) dictionary specifying the position of nodes
        """

        # find list of supplied node IDs
        nodeIDs = list(nodesDict.keys())
        for nodeID in nodesDict:
            nodeIDs += nodesDict[nodeID]
        nodeIDs = sorted(list(set(nodeIDs)))


        # ensure that node labels are the first positive integers
        if list(range(1, len(nodeIDs) + 1)) != nodeIDs:
            raise ValueError("Node labels must be consecutive positive inetegers that include 1.")
        # ensure dictionary does not map a node to itself
        for nodeID in nodesDict:
            if nodeID in nodesDict[nodeID]:
                raise ValueError("Node in `nodesDict` maps to itself: not allowed.")
        # ensure that k > 0, Q < 0, and g < 0
        if (k <= 0):
            raise ValueError("k must be positive.")
        if (Q >= 0):
            raise ValueError("Q must be negative.")
        # if (g >= 0):
        #     raise ValueError("g must be negative.")


        # deal with nodePosDict if it is not empty
        if len(nodePosDict) != 0:
            # ensure that nodePosDict defines a position for every node specified in `nodesDict`
            if sorted(list(nodePosDict.keys())) != nodeIDs:
                raise ValueError("`nodePosDict` must define positions for every node specified in `nodesDict`")
            # ensure that no two positions are the same
            if len(set(nodePosDict.values())) != len(nodeIDs):
                raise ValueError("No two nodes may have the same position")

            # set list of node objects
            self.nodes = [Node(nodeID, pos = np.array(nodePosDict[nodeID], dtype=np.float)) for nodeID in nodeIDs]

        # otherwise, no positions are supplied
        else:
            self.nodes = [Node(nodeID) for nodeID in nodeIDs]
            self.encircle_nodes()


        # create a bidirectional dictionary of node numbers
        for nodeID in nodeIDs:
            if nodeID not in nodesDict:
                nodesDict[nodeID] = []
        graphNodesDict = {}
        for nodeIDa in nodeIDs:
            connected_to_a = lambda nodeIDb: (nodeIDa in nodesDict[nodeIDb]) or (nodeIDb in nodesDict[nodeIDa])
            graphNodesDict[nodeIDa] = [nodeIDb for nodeIDb in nodesDict if (connected_to_a(nodeIDb))]

        # add connections to nodes using dictionary and edges to springboard object
        self.edges = []
        for nodeA in self.nodes:
            nodeA.add_connections([self.nodes[nodeIDb - 1] for nodeIDb in graphNodesDict[nodeA.ID]])
            self.edges += [(nodeA, self.nodes[nodeIDb - 1]) for nodeIDb in graphNodesDict[nodeA.ID] if nodeA.ID  < nodeIDb]

        # set spring and field constants
        self.k = k
        self.Q = Q
        # self.g = g

    def _increment(self, deltaT: float):
        """
        Increment timestep simulation by one step

        Args:
            deltaT - simulation time step
        """

        for node in self.nodes:
            change = np.array([0.0, 0.0])

            # add the spring forces
            for connection in node.connections:
                deltaD = deltaT ** 2 * self.k * (1 - node.distance_to(connection)) / node.degree
                vec = node.pos - connection.pos
                vec *= deltaD / np.sqrt(vec[0] ** 2 + vec[1] ** 2)
                change += vec

            # add the repellant forces
            for other in self.nodes:
                if node != other:
                    deltaD = self.Q * (deltaT / node.distance_to(other)) ** 2 * other.degree
                    vec = other.pos - node.pos
                    vec *= deltaD / np.sqrt(vec[0] ** 2 + vec[1] ** 2)
                    change += vec

            # add gravitational force
            # change += node.pos * node.degree * self.g / np.sqrt(0.001 + sum(node.pos ** 2))

            # set displacemnts
            node.pos += change

            # the "second order backwards" appromimation step
            # node.pos += node.pastPos[0] - node.pastPos[1]
            # node.pastPos[1] = node.pastPos[0].copy()
            # node.pastPos[0] = node.pos.copy()

    def move(self, deltaT: float, n: int):
        """
        Iterate _increment()

        Args:
            deltaT - simulation time step
            n - number of time steps
        """

        for _ in range(n):
            self._increment(deltaT)

    def plot(self, saveAs: str = ""):
        """
        Plot the graph

        Args:
            saveAs - (optional) file path to save
        """
        fig, ax = plt.subplots(figsize=(7, 7))
        ax.set_aspect("equal")
        ax.autoscale()
        for (nodeA, nodeB) in self.edges:
            ax.annotate("", xytext=nodeA.pos, xy=nodeB.pos, arrowprops={"arrowstyle": "-"}, va="center")
        for node in self.nodes:
            x = [node.pos[0] for node in self.nodes]
            y = [node.pos[1] for node in self.nodes]
        ax.plot(x, y, "o")
        plt.show()
        if saveAs != "":
            plt.savefig(saveAs)

    def move_plot_save(self, deltaT: float, n: int, saveAs: str):
        """
        Move forward in simulation and save image after each timestep

        Args:
            deltaT - simulation time step
            n - number of time steps
            saveAs - folder path to save images
        """
        for i in range(n):
            self.plot(f"{saveAs}/{i}")
            self._increment(deltaT)
        self.plot(f"{saveAs}/{n}")

    def settle(self, deltaT: float):
        """
        Increment timestep simulation until objects have settled

        Args:
            deltaT - simulation time step
        """
        sumDiff = 1  # > 0.05
        while sumDiff > 0.05:
            last = {node.ID: node.pos for node in self.nodes}
            self.move(deltaT, 500)
            sumDiff = sum([abs(node.pos[0] - last[node.ID][0]) + abs(node.pos[1] - last[node.ID][1]) for node in self.nodes])

    def random_reset(self):
        """
        Randomly reset node positions
        """
        for node in self.nodes:
            node.pos = np.array([r.uniform(0, 1), r.uniform(0, 1)])

    def encircle_nodes(self):
        """
        Arrange node positions into a circle.
        """
        for (node, i) in zip(self.nodes, range(len(self.nodes))):
            arg = 2 * np.pi * i / len(self.nodes)
            node.set_coord(np.array([np.cos(arg), np.sin(arg)]))

    def animate(self, deltaT: float, numFrames: int, movesPerFrame: int, xlim: float, ylim: float, size: Tuple[int, int]):
        """
        Increment timestep simulation until objects have settled

        Args:
            deltaT - simulation time step
            numFrames - number of frames in the animation
            movesPerFrame - value of `n` when `move(deltaT, n)` is called between frames
            xlim - [X lower bound, X upper bound]
            ylim - [Y lower bound, Y upper bound]
            size - figsize (x,y)
        Return:
            animation object
        """

        fig, ax = plt.subplots(figsize=size)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        edgeLines = [None] * len(self.edges)
        for i in range(len(self.edges)):
            A, B = self.edges[i]
            edgeLines[i], = ax.plot([A.pos[0], B.pos[0]],[A.pos[1], B.pos[1]])
        nodePoints, = ax.plot([node.pos[0] for node in self.nodes],[node.pos[1] for node in self.nodes], "o")

        def _next_frame(start):
            nodePoints.set_data([node.pos[0] for node in self.nodes],[node.pos[1] for node in self.nodes])
            for i in range(len(self.edges)):
                A, B = self.edges[i]
                edgeLines[i].set_data([A.pos[0], B.pos[0]],[A.pos[1], B.pos[1]])
            if start:
                self.move(deltaT, movesPerFrame)
            start = True

        start = False

        return animation.FuncAnimation(fig, _next_frame, fargs = (start), frames=numFrames, interval=30)


class Graph(object):
    def __init__(self, nodesDict: dict, isDigraph: bool = False):
        """
        Construct Graph object

        Args:
            nodesDict - adjacency of first positive integers
            isDigraph (bool) - boolean value to declare Graph type
        """

        self.nodesDict = nodesDict
        self.isDigraph = isDigraph

        # check to make sure that the keys and values are not skipping any
        # positive integers and that they are the first positive integers
        testIDs = []
        for key in nodesDict:
            testIDs += [key] + [value for value in nodesDict[key]]
        if set(testIDs) != set(range(1, len(list(set(testIDs))) + 1)):
            raise ValueError("Error in node keys and values: " + "missing number or disallowed character.")

        # make sure that the nodes dictionary is bidirectional
        graphNodesDict = {}
        connectedToA = lambda nodeB: (nodeA in nodesDict[nodeB]) or (nodeB in nodesDict[nodeA])
        for nodeA in nodesDict:
            graphNodesDict[nodeA] = [nodeB for nodeB in nodesDict if (connectedToA(nodeB) and nodeA != nodeB)]

        # make adjacency matrix
        if isDigraph:
            self.adjacencyMatrix = np.vstack([np.array([1 if nodeB in self.nodesDict[nodeA] else 0 for nodeB in self.nodesDict]) for nodeA in self.nodesDict]).T
        else:
            self.adjacencyMatrix = np.vstack([np.array([1 if nodeB in graphNodesDict[nodeA] else 0 for nodeB in graphNodesDict]) for nodeA in graphNodesDict]).T


        # use SpringBoard to find good coordinates
        self.springBoard = SpringBoard(graphNodesDict, 1, -1, -0.001)
        self.springBoard.move(0.1, 8000)
        self._normalize_pos()

    def _normalize_pos(self):
        """
        Normalize the positions springboard nodes for plotting
        """
        # collect all X and Y coordinates
        X = [node.pos[0] for node in self.springBoard.nodes]
        Y = [node.pos[1] for node in self.springBoard.nodes]
        # sutract out minmum of each
        for node in self.springBoard.nodes:
            node.pos -= np.array([min(X), min(Y)])
        # recollect all X and Y coordinates
        X = [node.pos[0] for node in self.springBoard.nodes]
        Y = [node.pos[1] for node in self.springBoard.nodes]
        # Scale by a little more than the max of each collection, X and Y
        for node in self.springBoard.nodes:
            node.pos = np.array([node.pos[0] / (max(X) + 1), node.pos[1] / (max(Y) + 1)])

    def plot(self, saveAs: str = "_"):
        """
        Plot the Graph

        Args:
            saveAs - (optional) a file path to save the plot
        """

        fig, ax = plt.subplots(figsize=(7, 7))
        plt.axis("off")
        ax.set_aspect("equal")
        r = 0.04

        for node in self.springBoard.nodes:
            X1, Y1 = node.pos[0], node.pos[1]
            # TODO: structure allows for other names, but circles won't adjust

            # add circle
            ax.add_artist(plt.Circle((X1, Y1), r, color="b", fill=False, clip_on=False))
            ax.text(X1, Y1, str(node.ID), fontsize=15, horizontalalignment="center", verticalalignment="center")

            # add lines per circle
            if self.isDigraph: # arrows
                for connectionIDNumber in self.nodesDict[node.ID]:
                    connection = self.springBoard.nodes[connectionIDNumber - 1]
                    X2, Y2 = connection.pos[0], connection.pos[1]
                    d = np.sqrt((X2 - X1) ** 2 + (Y2 - Y1) ** 2)
                    ax.annotate("", xytext=(X1, Y1), xy=(X2, Y2), arrowprops={"width": 0.01, "shrink": 1.2 * r / d})
            else: # lines
                for connection in node.connections:
                    if node.ID < connection.ID: # this makes each connection only graph once
                        X2, Y2 = connection.pos[0], connection.pos[1]
                        d = np.sqrt((X2 - X1) ** 2 + (Y2 - Y1) ** 2)
                        x = r * ((X2 - X1) / d)
                        y = r * ((Y2 - Y1) / d)
                        ax.annotate("", xytext=(X1 + x, Y1 + y), xy=(X2 - x, Y2 - y), arrowprops={"arrowstyle": "-"})
        if saveAs != "_":
            plt.savefig(saveAs)

    def force(self, n: int):
        """
        Move forward time step simulation
        n - number of time steps
        """
        self.springBoard.move(0.1, n)
        self._normalize_pos()

    def random_reset(self):
        """
        Randomly reset node positions and let time step simulation resettle
        """
        self.springBoard.random_reset()
        self.springBoard.move(0.1,8000)
        self._normalize_pos()
