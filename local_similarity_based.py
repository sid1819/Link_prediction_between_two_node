import random
import sys
import math

from collections import deque

random_links_to_be_predicted = {}
NumberOfNOdes = [-1]
NumberOfEdges = [-1]
graph = {}
num = [-1]
delta = 10 ** -6

def re():
    file = open("karate.net", "r+")
    ar = file.readline().split()
    NumberOfNOdes[0], NumberOfEdges[0] = int(ar[0]), int(ar[1])

    graph.clear()
    for i in range(1, NumberOfNOdes[0] + 1):
        graph[i] = set([])

    for _ in range(NumberOfEdges[0]):
        ar = file.readline().split()
        graph[int(ar[0])].add(int(ar[1]))
        graph[int(ar[1])].add(int(ar[0]))

    file.close()

    random_links_to_be_predicted.clear()
    num[0] = 0
    for _ in range(max(NumberOfEdges[0] // 10, 2)):
        x = random.randint(1, NumberOfNOdes[0])
        y = random.randint(1, NumberOfNOdes[0])
        if x == y:
            continue
        x, y = min(x, y), max(x, y)
        if (x, y) in random_links_to_be_predicted:
            continue
        num[0] += 1
        if y in graph[x]:
            graph[x].remove(y)
            graph[y].remove(x)
            random_links_to_be_predicted[(x, y)] = 0
        else:
            random_links_to_be_predicted[(x, y)] = 1


re()
NetworkDensity = NumberOfEdges[0] / (NumberOfNOdes[0] * (NumberOfNOdes[0] - 1) / 2)

def CommonNeighbour(a, b):
    s1 = set(graph[a])
    s2 = set(graph[b])
    common = s1.intersection(s2)
    return len(common)


def JaccardCoefficient(a, b):
    s1 = set(graph[a])
    s2 = set(graph[b])
    common = s1.union(s2)
    return CommonNeighbour(a, b) / len(common)


 
def Preferentialattachment(a, b):
    return len(graph[a]) * len(graph[b])



def AdamicAdarCoefficient(a, b):
    val = 0
    s1 = set(graph[a])
    s2 = set(graph[b])
    common = s1.intersection(s2)
    for node in common:
        val += 1 / math.log(len(graph[node]))
    return val

 
def ResourceAllocationIndex(a, b):
    val = 0
    s1 = set(graph[a])
    s2 = set(graph[b])
    common = s1.intersection(s2)
    for node in common:
        val += 1 / len(graph[node])
    return val


 
def SaltonIndex(a, b):
    return CommonNeighbour(a, b) / (math.sqrt(len(graph[a]) * len(graph[b])) + delta)

 
def SorensenIndex(a, b):
    return 2 * (CommonNeighbour(a, b) / ((len(graph[a]) * len(graph[b])) + delta))



def Gammma(common, z):
    val = 0
    for node in graph[z]:
        if node in common:
            val += 1
    return val

 
def CARBasedCommonNeighborIndex(a, b):
    s1 = set(graph[a])
    s2 = set(graph[b])
    common = s1.intersection(s2)
    val = 0
    for node in common:
        val += (Gammma(common, node) / 2)
    return val * CommonNeighbour(a, b)


 

def CARBasedAdamicAdarIndex(a, b):
    s1 = set(graph[a])
    s2 = set(graph[b])
    common = s1.intersection(s2)
    val = 0
    for node in common:
        val += ((Gammma(common, node) / 2) / math.log2(len(graph[node])))
    return val

 
def CARBasedResourceAllocationIndex(a, b):
    s1 = set(graph[a])
    s2 = set(graph[b])
    common = s1.intersection(s2)
    val = 0
    for node in common:
        val += ((Gammma(common, node) / 2) / len(graph[node]))
    return val


def CAR(a, b):
    val = 0
    s1 = set(graph[a])
    s2 = set(graph[b])
    common = s1.intersection(s2)
    for node in common:
        val += Gammma(common, node)
    val //= 2
    return CommonNeighbour(a, b) * val


 
def CARBasedPreferentialAttachmentIndex(a, b):
    s1 = set(graph[a])
    s2 = set(graph[b])
    common = s1.intersection(s2)
    ea = len(s1 - common)
    eb = len(s2 - common)
    ca = CAR(a, b)
    return ea * eb + ea * ca + eb * ca + ca * ca



def HubPromotedIndex(a, b):
    s1 = set(graph[a])
    s2 = set(graph[b])
    common = s1.intersection(s2)
    return len(common) / min(len(graph[a]), len(graph[b]))


 
def HubDepressedIndex(a, b):
    s1 = set(graph[a])
    s2 = set(graph[b])
    common = s1.intersection(s2)
    return len(common) / max(len(graph[a]), len(graph[b]))



def ClusteringCoefficient(z):
    s = set(graph[z])
    val = 0
    for node in s:
        for i in graph[node]:
            if i in s:
                val += 1
    val //= 2
    return val / (len(graph[z]) * (len(graph[z]) - 1))


 
def LocalNaiveBayesBasedCommonNeighbors(a, b):
    s1 = set(graph[a])
    s2 = set(graph[b])
    common = s1.intersection(s2)
    val = 0
    for node in common:
        val += math.log(ClusteringCoefficient(node) / (1 - ClusteringCoefficient(node))) + math.log((1 - NetworkDensity) / NetworkDensity)
    return val

 

def LeichtHolmeNewmanLocalIndex(a, b):
    s1 = set(graph[a])
    s2 = set(graph[b])
    common = s1.intersection(s2)
    return len(common) / (len(s1) * len(s2))



def NodeClusteringCoefficient(a, b):
    s1 = set(graph[a])
    s2 = set(graph[b])
    common = s1.intersection(s2)
    val = 0
    for node in common:
        val += ClusteringCoefficient(node)
    return val


 
def NodeAndLinkClusteringCoefficient(a, b):
    s1 = set(graph[a])
    s2 = set(graph[b])
    common = s1.intersection(s2)
    val = 0
    for node in common:
        c = ClusteringCoefficient(node)
        s3 = set(graph[node])
        k = len(graph[node]) - 1
        val += (len(s1.intersection(s3)) / k) * c + (len(s2.intersection(s3)) / k) * c
    return val



v = [CommonNeighbour, JaccardCoefficient, Preferentialattachment, AdamicAdarCoefficient, ResourceAllocationIndex, SaltonIndex, SorensenIndex, CARBasedCommonNeighborIndex, CARBasedAdamicAdarIndex, CARBasedResourceAllocationIndex, CARBasedPreferentialAttachmentIndex, HubPromotedIndex, HubDepressedIndex, LocalNaiveBayesBasedCommonNeighbors, LeichtHolmeNewmanLocalIndex, NodeClusteringCoefficient, NodeAndLinkClusteringCoefficient]
for i in v:
    re()
    print(i.__name__)
    metric = []
    cur = 0
    for node1, node2 in random_links_to_be_predicted:
        metric.append([i(node1, node2), random_links_to_be_predicted[(node1, node2)]])
        cur += random_links_to_be_predicted[(node1, node2)]
    metric.sort()
    ac = 0
    while len(metric) > 0:
        a, b = metric.pop()
        if cur > 0:
            cur -= 1
            if b == 1:
                ac += 1
        else:
            if b == 0:
                ac += 1
    print("Accuracy : ", ac / num[0])
