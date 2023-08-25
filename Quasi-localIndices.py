import math
import random
import sys

from collections import deque, Counter

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




def LocalPathIndex(a, b):
    epsilon = 0.25
    s1 = set(graph[a])
    s2 = set(graph[b])
    common = s1.intersection(s2)
    c = Counter()
    for node in s1:
        for i in graph[node]:
            c[i] += 1
    val = 0
    for node in s2:
        val += c[node]
    return len(common) + epsilon * val


def CH2_L2(a, b):
    val = 0
    s1 = set(graph[a])
    s2 = set(graph[b])
    common = s1.intersection(s2)
    for node in common:
        s3 = set(graph[node])
        p = len(s3.intersection(common))
        val += (1 + p) / (1 + len(s3) - p - 2)
    return val


def CH2_L3(a, b):
    val = 0
    s1 = list(graph[a])
    s2 = list(graph[b])
    s = set([])
    for node in s1:
        for i in graph[node]:
            if i in s2:
                s.add(node)
                s.add(i)

    if a in s:
        s.remove(a)
    if b in s:
        s.remove(b)

    cz1, cz2, oz1, oz2 = [], [], [], []
    for node in s1:
        s3 = s.intersection(graph[node])
        cz1.append(len(s3))
        flag = 0
        if a in graph[node]:
            flag += 1
        if b in graph[node]:
            flag += 1
        oz1.append(len(graph[node]) - len(s3) - flag)

    for node in s2:
        s3 = s.intersection(graph[node])
        cz2.append(len(s3))
        flag = 0
        if a in graph[node]:
            flag += 1
        if b in graph[node]:
            flag += 1
        oz2.append(len(graph[node]) - len(s3) - flag)

    for i in range(len(s1)):
        for j in range(len(s2)):
            if s2[j] not in graph[s1[i]]:
                continue
            val += math.sqrt((1 + cz1[i]) * (1 + cz2[j])) / math.sqrt((1 + oz1[i]) * (1 + oz2[j]))

    return val


def MatrixTranspose(Ma):
    Tra = []
    for _ in range(len(Ma[0])):
        Tra.append([])

    for i in range(len(Ma[0])):
        for j in range(len(Ma)):
            Tra[i].append(Ma[j][i])

    return Tra


def SuperposedRandomWalk(a, b):
    P = []
    for node in range(1, NumberOfNOdes[0] + 1):
        q = []
        for i in range(1, NumberOfNOdes[0] + 1):
            if i in graph[node]:
                q.append(1 / len(graph[node]))
            else:
                q.append(0)
        P.append(q)

    PT = MatrixTranspose(P)
    Pia = [[0 for _ in range(NumberOfNOdes[0])]]
    Pia[0][a - 1] = 1
    for _ in range(5):
        q = []
        for i in range(NumberOfNOdes[0]):
            val = 0
            for j in range(NumberOfNOdes[0]):
                val += PT[i][j] * Pia[-1][j]
            q.append(val)
        Pia.append(q)

    Pib = [[0 for _ in range(NumberOfNOdes[0])]]
    Pib[0][b - 1] = 1
    for _ in range(5):
        q = []
        for i in range(NumberOfNOdes[0]):
            val = 0
            for j in range(NumberOfNOdes[0]):
                val += PT[i][j] * Pib[-1][j]
            q.append(val)
        Pib.append(q)

    SRW = 0
    for t in range(1, 5):
        LRWt = ((len(graph[a]) * Pia[t][b - 1]) + (len(graph[b]) * Pib[t][a - 1])) / (2 * NumberOfEdges[0])
        SRW +=  LRWt
    return SRW


v = [LocalPathIndex, CH2_L2, CH2_L3, SuperposedRandomWalk]
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

