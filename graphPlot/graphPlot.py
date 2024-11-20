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
        self.disp = np.array([0.0, 0.0])
        self.prevDisp = np.array([0.0, 0.0])
        self.fixed = False

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
        # np.sqrt((self.pos[0] - node.pos[0]) ** 2 + (self.pos[1] - node.pos[1]) ** 2)
        return np.sqrt(np.sum((self.pos - node.pos) ** 2))

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
        return f"<{self.ID} @ ({round(self.pos[0],3)},{round(self.pos[1],3)})>"



class SpringBoard(object):
    def __init__(self, nodesDict: Dict[int, List[int]], k: float = 1, Q: float = -1, nodePosDict: Dict[int, Tuple[float, float]] = {}):
        """
        Construct a SpringBoard object.
        If the nodes are all centered at the orgin, spread them out.

        Args:
            nodesDict - an adjacency dictionary of first positive integers
            k - (optional, 1 by default) coefficient of spring
            Q - (optional, -1 by default) coefficient of electric field
            nodePosDict - (optional) dictionary specifying the position of nodes
        """

        # find list of supplied node IDs
        self.nodeIDs = list(nodesDict.keys())
        for nodeID in nodesDict:
            self.nodeIDs += nodesDict[nodeID]
        self.nodeIDs = sorted(list(set(self.nodeIDs)))


        # ensure that node labels are the first positive integers
        if list(range(1, len(self.nodeIDs) + 1)) != self.nodeIDs:
            raise ValueError("Node labels must be consecutive positive inetegers that include 1. Currently, node labels are " + str(self.nodeIDs))
        # ensure dictionary does not map a node to itself
        for nodeID in nodesDict:
            if nodeID in nodesDict[nodeID]:
                raise ValueError(f"Node in `nodesDict` maps {nodeID} to itself: not allowed.")
        # ensure that k > 0, Q < 0
        if (k <= 0):
            raise ValueError("k must be positive.")
        if (Q >= 0):
            raise ValueError("Q must be negative.")


        # deal with nodePosDict if it is not empty
        if len(nodePosDict) != 0:
            # ensure that no two positions are the same
            if len(set(nodePosDict.values())) != len(nodePosDict.values()):
                raise ValueError("No two nodes may have the same position")
            # ensure that only nodeIDs are given positions
            for nodeID in nodePosDict:
                if nodeID not in self.nodeIDs:
                    raise ValueError(f"{nodeID} is not a node in this SpringBoard")

            # make sure that nodePosDict defines a position for every node specified in `nodesDict`
            if sorted(list(nodePosDict.keys())) != self.nodeIDs:
                for nodeID in self.nodeIDs:
                    if nodeID not in nodePosDict:
                        testPos = (r.uniform(0,1), r.uniform(0,1))
                        while testPos in nodePosDict.keys():
                            testPos = (r.uniform(0,1), r.uniform(0,1))
                        nodePosDict[nodeID] = testPos

            # set list of node objects
            self.nodes = [Node(nodeID, pos = np.array(nodePosDict[nodeID], dtype=float)) for nodeID in self.nodeIDs]

        # otherwise, no positions are supplied
        else:
            self.nodes = [Node(nodeID) for nodeID in self.nodeIDs]
            self.encircle_nodes()


        # fill in mapped to but not mapped from ids
        self.nodesDict = nodesDict
        for nodeID in self.nodeIDs:
            if nodeID not in self.nodesDict:
                self.nodesDict[nodeID] = []
        # create a bidirectional dictionary of node numbers
        self.graphNodesDict = {}
        for nodeIDa in self.nodeIDs:
            connected_to_a = lambda nodeIDb: (nodeIDa in nodesDict[nodeIDb]) or (nodeIDb in nodesDict[nodeIDa])
            self.graphNodesDict[nodeIDa] = [nodeIDb for nodeIDb in nodesDict if (connected_to_a(nodeIDb))]

        # add connections to nodes using dictionary and edges to springboard object
        self.edges = []
        for nodeA in self.nodes:
            nodeA.add_connections([self.nodes[nodeIDb - 1] for nodeIDb in self.graphNodesDict[nodeA.ID]])
            self.edges += [(nodeA, self.nodes[nodeIDb - 1]) for nodeIDb in self.graphNodesDict[nodeA.ID] if nodeA.ID  < nodeIDb]

        # set spring and field constants
        self.k = k
        self.Q = Q

    def _increment(self, deltaT: float):
        """
        Increment timestep simulation by one step

        Args:
            deltaT - simulation time step
        """

        for node in filter(lambda node: not node.fixed, self.nodes):

            # add the repellant forces
            rep_change = np.array([0.0, 0.0])
            for other in self.nodes:
                if other != node:
                    rep_change += other.degree * (other.pos - node.pos) / (node.distance_to(other)**3)

            # add the spring forces
            spring_change = np.array([0.0, 0.0])
            for connection in node.connections:
                dist = node.distance_to(connection)
                spring_change += (connection.pos - node.pos) * (dist - 1) / (node.degree * dist)

            node.disp = deltaT**2 * (self. Q * rep_change + self.k * spring_change)

        # set displacemnts
        for node in filter(lambda node: not node.fixed, self.nodes):
            node.pos = node.pos + node.disp + node.prevDisp
            node.prevDisp = node.disp
            node.disp = np.array([0.0, 0.0])
            

    def move(self, deltaT: float, n: int):
        """
        Iterate _increment()

        Args:
            deltaT - simulation time step
            n - number of time steps
        """

        for _ in range(n):
            self._increment(deltaT)

    def plot(self, size: Tuple[int, int] = (7,7), saveAs: str = ""):
        """
        Plot the graph

        Args:
            size - (optional, (7,7) by default) size tuple for plot image
            saveAs - (optional) file path to save
        """
        fig, ax = plt.subplots(figsize=size)
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

    def get_fixed_nodes(self):
        """
        find list of nodes that are fixed

        Return:
            list of node objects that are fixes
        """
        return list(filter(lambda node: node.fixed, self.nodes))

    def set_node_pos(self, nodeID: int, pos: Tuple[float, float]):
        """
        Set a node position

        Args:
            nodeID - ID integer of the node to be set
            pos - tuple of floats defining the position to set the node
        """
        self.nodes[nodeID - 1].pos = np.array(pos)

    def fix_nodes(self, nodeIDList: List[int]):
        """
        fix node positions

        Args:
            nodeIDList - list of nodeID ints to be fixed
        """
        for nodeID in nodeIDList:
            if nodeID not in self.nodeIDs:
                raise ValueError(f"Node ID {nodeID} does not match a node in this SpringBoard. As a result, no nodes were fixed.")
        for nodeID in nodeIDList:
            self.nodes[nodeID - 1].fixed = True


class Graph(object):
    def __init__(self, nodesDict: dict, isDigraph: bool = False):
        """
        Construct Graph object

        Args:
            nodesDict - adjacency of first positive integers
            isDigraph (bool) - boolean value to declare Graph type
        """

        self.isDigraph = isDigraph

        # use SpringBoard to find good coordinates
        self.springBoard = SpringBoard(nodesDict, 1, -1)
        self.springBoard.move(0.1, 8000)
        self._normalize_pos()

        self.nodesDict = dict(self.springBoard.nodesDict)

        self.graphNodesDict = dict(self.springBoard.graphNodesDict)

        # make adjacency matrix
        if self.isDigraph:
            self.adjacencyMatrix = np.vstack([np.array([1 if nodeB in self.nodesDict[nodeA] else 0 for nodeB in self.nodesDict]) for nodeA in self.nodesDict]).T
        else:
            self.adjacencyMatrix = np.vstack([np.array([1 if nodeB in self.graphNodesDict[nodeA] else 0 for nodeB in self.graphNodesDict]) for nodeA in self.graphNodesDict]).T

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













