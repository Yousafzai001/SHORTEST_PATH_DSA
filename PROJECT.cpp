#include <iostream>
#include <limits>
#include <climits>
#include <iomanip>
#include <string>

using namespace std;

/*
 * Forward declarations
 */
class Edge;
class Vertex;
class DynamicEdgeArray;
class DynamicVertexArray;
class Graph;
class MinHeapNode;
class MinHeap;
class ShortestPathFinder;

// ===========================
// EDGE CLASS
// ===========================
class Edge {
private:
    int dest;      // index of destination vertex in Graph's vertex array
    int weight;    // weight of this edge

public:
    Edge() : dest(-1), weight(0) {}
    Edge(int d, int w) : dest(d), weight(w) {}

    int getDest() const { return dest; }
    int getWeight() const { return weight; }

    void setDest(int d) { dest = d; }
    void setWeight(int w) { weight = w; }
};

// ===========================
// DYNAMIC ARRAY OF EDGES
// ===========================
class DynamicEdgeArray {
private:
    Edge* arr;         // pointer to dynamic array holding Edge objects
    int capacity;      // current capacity of the array
    int size;          // number of valid elements in the array

    void resize() {
        int newCapacity = capacity * 2;
        Edge* newArr = new Edge[newCapacity];
        for (int i = 0; i < size; ++i) {
            newArr[i] = arr[i];
        }
        delete[] arr;
        arr = newArr;
        capacity = newCapacity;
    }

public:
    DynamicEdgeArray() {
        capacity = 4;
        size = 0;
        arr = new Edge[capacity];
    }

    ~DynamicEdgeArray() {
        delete[] arr;
    }

    int getSize() const {
        return size;
    }

    bool isEmpty() const {
        return size == 0;
    }

    Edge get(int index) const {
        if (index < 0 || index >= size) {
            return Edge(-1, 0);
        }
        return arr[index];
    }

    void pushBack(const Edge& e) {
        if (size == capacity) {
            resize();
        }
        arr[size++] = e;
    }

    // Remove element at index (shifts subsequent elements left)
    void removeAt(int index) {
        if (index < 0 || index >= size) {
            return;
        }
        for (int i = index; i < size - 1; ++i) {
            arr[i] = arr[i + 1];
        }
        --size;
    }

    // Find index of edge whose destination equals destID; return -1 if not found
    int findEdgeIndex(int destID) const {
        for (int i = 0; i < size; ++i) {
            if (arr[i].getDest() == destID) {
                return i;
            }
        }
        return -1;
    }

    // Update weight if edge with destID exists; return true if updated, false otherwise
    bool updateEdgeWeight(int destID, int newWeight) {
        int idx = findEdgeIndex(destID);
        if (idx == -1) return false;
        arr[idx].setWeight(newWeight);
        return true;
    }

    // Display all edges in this adjacency list
    void display() const {
        for (int i = 0; i < size; ++i) {
            cout << " -> [Dest IDX: " << arr[i].getDest()
                << ", W: " << arr[i].getWeight() << "]";
        }
    }
};

// ===========================
// VERTEX CLASS
// ===========================
class Vertex {
private:
    int id;                     // Unique ID to represent this vertex
    DynamicEdgeArray adjList;   // adjacency list of edges

public:
    Vertex() : id(-1) {}
    Vertex(int idVal) : id(idVal) {}

    int getID() const { return id; }
    void setID(int idVal) { id = idVal; }

    DynamicEdgeArray& getAdjList() { return adjList; }
    const DynamicEdgeArray& getAdjList() const { return adjList; }

    // Add an edge from this vertex to destination index with weight
    void addEdge(int destIdx, int weight) {
        int existing = adjList.findEdgeIndex(destIdx);
        if (existing != -1) {
            // Edge already exists; update weight
            adjList.updateEdgeWeight(destIdx, weight);
        }
        else {
            Edge e(destIdx, weight);
            adjList.pushBack(e);
        }
    }

    // Remove an edge to destIdx if exists
    bool removeEdge(int destIdx) {
        int idx = adjList.findEdgeIndex(destIdx);
        if (idx == -1) return false;
        adjList.removeAt(idx);
        return true;
    }

    // Update weight of edge to destIdx
    bool updateEdge(int destIdx, int newWeight) {
        return adjList.updateEdgeWeight(destIdx, newWeight);
    }
};

// ===========================
// DYNAMIC ARRAY OF VERTICES
// ===========================
class DynamicVertexArray {
private:
    Vertex** arr;      // pointer to dynamic array holding pointers to Vertex
    int capacity;      // capacity of the array
    int size;          // current number of vertices

    void resize() {
        int newCapacity = capacity * 2;
        Vertex** newArr = new Vertex * [newCapacity];
        for (int i = 0; i < size; ++i) {
            newArr[i] = arr[i];
        }
        delete[] arr;
        arr = newArr;
        capacity = newCapacity;
    }

public:
    DynamicVertexArray() {
        capacity = 4;
        size = 0;
        arr = new Vertex * [capacity];
    }

    ~DynamicVertexArray() {
        for (int i = 0; i < size; ++i) {
            delete arr[i];
        }
        delete[] arr;
    }

    int getSize() const { return size; }

    bool isEmpty() const { return size == 0; }

    Vertex* get(int index) const {
        if (index < 0 || index >= size) return nullptr;
        return arr[index];
    }

    // Add a new Vertex pointer (ownership is transferred)
    void pushBack(Vertex* v) {
        if (size == capacity) {
            resize();
        }
        arr[size++] = v;
    }

    // Remove vertex at index; also delete object
    void removeAt(int index) {
        if (index < 0 || index >= size) return;
        delete arr[index];
        for (int i = index; i < size - 1; ++i) {
            arr[i] = arr[i + 1];
        }
        --size;
    }

    // Find index of vertex with given ID, using linear search
    int findVertexIndexLinear(int vertexID) const {
        for (int i = 0; i < size; ++i) {
            if (arr[i]->getID() == vertexID) return i;
        }
        return -1;
    }

    // Sort vertices by ID (ascending) using QuickSort
    void quickSortVertices(int left, int right) {
        if (left >= right) return;
        int pivot = arr[(left + right) / 2]->getID();
        int i = left, j = right;
        while (i <= j) {
            while (arr[i]->getID() < pivot) i++;
            while (arr[j]->getID() > pivot) j--;
            if (i <= j) {
                Vertex* temp = arr[i];
                arr[i] = arr[j];
                arr[j] = temp;
                i++;
                j--;
            }
        }
        if (left < j) quickSortVertices(left, j);
        if (i < right) quickSortVertices(i, right);
    }

    // MergeSort helper: merge two sorted halves
    void merge(int left, int mid, int right) {
        int n1 = mid - left + 1;
        int n2 = right - mid;
        Vertex** L = new Vertex * [n1];
        Vertex** R = new Vertex * [n2];
        for (int i = 0; i < n1; ++i) L[i] = arr[left + i];
        for (int j = 0; j < n2; ++j) R[j] = arr[mid + 1 + j];
        int i = 0, j = 0, k = left;
        while (i < n1 && j < n2) {
            if (L[i]->getID() <= R[j]->getID()) {
                arr[k++] = L[i++];
            }
            else {
                arr[k++] = R[j++];
            }
        }
        while (i < n1) arr[k++] = L[i++];
        while (j < n2) arr[k++] = R[j++];
        delete[] L;
        delete[] R;
    }

    // MergeSort vertices by ID
    void mergeSortVertices(int left, int right) {
        if (left >= right) return;
        int mid = left + (right - left) / 2;
        mergeSortVertices(left, mid);
        mergeSortVertices(mid + 1, right);
        merge(left, mid, right);
    }

    // Public wrappers for sorting full array
    void sortVerticesQuick() {
        if (size > 0) quickSortVertices(0, size - 1);
    }
    void sortVerticesMerge() {
        if (size > 0) mergeSortVertices(0, size - 1);
    }

    // Binary Search for vertex ID (array must be sorted by ID)
    int findVertexIndexBinary(int vertexID) const {
        int left = 0, right = size - 1;
        while (left <= right) {
            int mid = left + (right - left) / 2;
            if (arr[mid]->getID() == vertexID) return mid;
            if (arr[mid]->getID() < vertexID) left = mid + 1;
            else right = mid - 1;
        }
        return -1;
    }

    // Display the IDs of all vertices
    void displayVertexIDs() const {
        cout << "[ ";
        for (int i = 0; i < size; ++i) {
            cout << arr[i]->getID();
            if (i < size - 1) cout << ", ";
        }
        cout << " ]" << endl;
    }
};

// ===========================
// GRAPH CLASS
// ===========================
class Graph {
private:
    DynamicVertexArray vertices;

    // Helper: Check if a vertex ID already exists
    bool vertexExists(int vertexID) const {
        return (vertices.findVertexIndexLinear(vertexID) != -1);
    }

public:
    Graph() {}
    ~Graph() {}

    // Create (add) a new vertex with given ID
    bool addVertex(int vertexID) {
        if (vertexExists(vertexID)) {
            return false;   // already exists
        }
        Vertex* v = new Vertex(vertexID);
        vertices.pushBack(v);
        return true;
    }

    // Remove a vertex and all edges pointing to or from it
    bool removeVertex(int vertexID) {
        int idx = vertices.findVertexIndexLinear(vertexID);
        if (idx == -1) return false;
        // Before removing vertex, remove edges in other adjacency lists that point to it
        for (int i = 0; i < vertices.getSize(); ++i) {
            if (i == idx) continue;
            vertices.get(i)->removeEdge(idx);
        }
        // Remove the vertex itself
        vertices.removeAt(idx);
        // Because indices changed, we must rebuild all adjacency lists: 
        // any dest index > idx must be decremented by 1.
        for (int i = 0; i < vertices.getSize(); ++i) {
            DynamicEdgeArray& arr = vertices.get(i)->getAdjList();
            for (int j = 0; j < arr.getSize(); ++j) {
                Edge e = arr.get(j);
                int d = e.getDest();
                if (d > idx) {
                    arr.get(j).setDest(d - 1);
                }
            }
        }
        return true;
    }

    // Read (retrieve) index of vertex with ID (linear search)
    int getVertexIndex(int vertexID) const {
        return vertices.findVertexIndexLinear(vertexID);
    }

    // Add or update an edge between two vertices
    bool addEdge(int srcID, int destID, int weight) {
        int srcIdx = getVertexIndex(srcID);
        int destIdx = getVertexIndex(destID);
        if (srcIdx == -1 || destIdx == -1) return false;
        vertices.get(srcIdx)->addEdge(destIdx, weight);
        return true;
    }

    // Remove edge between two vertices
    bool removeEdge(int srcID, int destID) {
        int srcIdx = getVertexIndex(srcID);
        int destIdx = getVertexIndex(destID);
        if (srcIdx == -1 || destIdx == -1) return false;
        return vertices.get(srcIdx)->removeEdge(destIdx);
    }

    // Update weight of an existing edge
    bool updateEdge(int srcID, int destID, int newWeight) {
        int srcIdx = getVertexIndex(srcID);
        int destIdx = getVertexIndex(destID);
        if (srcIdx == -1 || destIdx == -1) return false;
        return vertices.get(srcIdx)->updateEdge(destIdx, newWeight);
    }

    // Display the entire graph (each vertex with its adjacency list)
    void displayGraph() const {
        cout << "\nGraph Representation (Adjacency Lists):\n";
        for (int i = 0; i < vertices.getSize(); ++i) {
            cout << "Vertex ID " << vertices.get(i)->getID() << " (Idx " << i << "):";
            vertices.get(i)->getAdjList().display();
            cout << endl;
        }
        cout << endl;
    }

    // Return reference to DynamicVertexArray (for sorting/searching utilities)
    DynamicVertexArray& getVerticesArray() {
        return vertices;
    }

    // Return const reference
    const DynamicVertexArray& getVerticesArray() const {
        return vertices;
    }
};

// ===========================
// MIN-HEAP NODE CLASS
// ===========================
class MinHeapNode {
public:
    int v;      // vertex index
    int dist;   // distance (key)

    MinHeapNode() : v(-1), dist(INT_MAX) {}
    MinHeapNode(int vertex, int distance) : v(vertex), dist(distance) {}
};

// ===========================
// MIN-HEAP CLASS FOR DIJKSTRA
// ===========================
class MinHeap {
private:
    MinHeapNode* harr;   // pointer to array of heap nodes
    int capacity;        // maximum possible size of min heap
    int heapSize;        // current number of elements in min heap
    int* pos;            // pos[i] gives position of vertex i in harr array, needed for decreaseKey

    void swapMinHeapNode(int i, int j) {
        MinHeapNode temp = harr[i];
        harr[i] = harr[j];
        harr[j] = temp;

        // Update positions
        pos[harr[i].v] = i;
        pos[harr[j].v] = j;
    }

    // Heapify at index i (smallest distance at top)
    void minHeapify(int i) {
        int smallest = i;
        int left = 2 * i + 1;
        int right = 2 * i + 2;

        if (left < heapSize && harr[left].dist < harr[smallest].dist)
            smallest = left;
        if (right < heapSize && harr[right].dist < harr[smallest].dist)
            smallest = right;

        if (smallest != i) {
            swapMinHeapNode(i, smallest);
            minHeapify(smallest);
        }
    }

public:
    MinHeap(int cap) : capacity(cap), heapSize(0) {
        harr = new MinHeapNode[cap];
        pos = new int[cap];
    }

    ~MinHeap() {
        delete[] harr;
        delete[] pos;
    }

    bool isEmpty() const {
        return heapSize == 0;
    }

    // Insert a new node (vertex v with distance dist)
    void insertKey(int v, int dist) {
        if (heapSize == capacity) return; // cannot insert
        int i = heapSize++;
        harr[i] = MinHeapNode(v, dist);
        pos[v] = i;

        // Fix the min heap property if violated
        while (i != 0 && harr[(i - 1) / 2].dist > harr[i].dist) {
            int parentIdx = (i - 1) / 2;
            swapMinHeapNode(i, parentIdx);
            i = parentIdx;
        }
    }

    // Decrease key value for vertex v to newDist
    void decreaseKey(int v, int newDist) {
        int i = pos[v];
        harr[i].dist = newDist;

        while (i != 0 && harr[(i - 1) / 2].dist > harr[i].dist) {
            int parentIdx = (i - 1) / 2;
            swapMinHeapNode(i, parentIdx);
            i = parentIdx;
        }
    }

    // Extract node with minimum dist (root of the heap)
    MinHeapNode extractMin() {
        if (isEmpty()) return MinHeapNode(-1, INT_MAX);

        // Store the root node
        MinHeapNode root = harr[0];

        // Replace root with last node
        MinHeapNode lastNode = harr[heapSize - 1];
        harr[0] = lastNode;
        pos[lastNode.v] = 0;

        --heapSize;
        minHeapify(0);

        return root;
    }

    // Check if a given vertex 'v' is in this heap
    bool isInMinHeap(int v) const {
        if (pos[v] < heapSize) return true;
        return false;
    }
};

// ===========================
// SHORTEST PATH FINDER (DIJKSTRA)
// ===========================
class ShortestPathFinder {
private:
    const Graph* graph;

    // Utility to reconstruct path from parent array
    void printPath(int parent[], int j, const DynamicVertexArray& vertices) const {
        if (parent[j] == -1) {
            cout << vertices.get(j)->getID();
            return;
        }
        printPath(parent, parent[j], vertices);
        cout << " -> " << vertices.get(j)->getID();
    }

public:
    ShortestPathFinder(const Graph* g) : graph(g) {}

    // Dijkstra's algorithm from sourceID to targetID
    void dijkstra(int sourceID, int targetID) const {
        const DynamicVertexArray& vertices = graph->getVerticesArray();
        int V = vertices.getSize();
        if (V == 0) {
            cout << "Graph is empty.\n";
            return;
        }
        int srcIdx = vertices.findVertexIndexLinear(sourceID);
        int tgtIdx = vertices.findVertexIndexLinear(targetID);
        if (srcIdx == -1 || tgtIdx == -1) {
            cout << "Source or target vertex not found.\n";
            return;
        }

        // dist[i] holds shortest distance from src to i
        int* dist = new int[V];
        int* parent = new int[V];   // to store shortest path tree

        // Initialize all distances as infinite and parents as -1
        for (int i = 0; i < V; ++i) {
            dist[i] = INT_MAX;
            parent[i] = -1;
        }

        // MinHeap for vertices not yet processed
        MinHeap minHeap(V);

        // Insert all vertices into min heap with infinite distance initially
        for (int v = 0; v < V; ++v) {
            minHeap.insertKey(v, dist[v]);
        }

        // Make distance of source vertex as 0 so it is extracted first
        dist[srcIdx] = 0;
        minHeap.decreaseKey(srcIdx, 0);

        while (!minHeap.isEmpty()) {
            // Extract vertex with minimum distance value
            MinHeapNode minNode = minHeap.extractMin();
            int u = minNode.v;

            // If we reached target, break
            if (u == tgtIdx) break;

            // Iterate through all adjacent vertices of u
            const DynamicEdgeArray& adj = vertices.get(u)->getAdjList();
            for (int i = 0; i < adj.getSize(); ++i) {
                Edge e = adj.get(i);
                int v = e.getDest();
                int weight = e.getWeight();

                // If shorter path found
                if (minHeap.isInMinHeap(v) && dist[u] != INT_MAX && dist[u] + weight < dist[v]) {
                    dist[v] = dist[u] + weight;
                    parent[v] = u;
                    minHeap.decreaseKey(v, dist[v]);
                }
            }
        }

        // Print shortest distance and path
        if (dist[tgtIdx] == INT_MAX) {
            cout << "No path exists from " << sourceID << " to " << targetID << ".\n";
        }
        else {
            cout << "Shortest distance from " << sourceID << " to " << targetID << " is " << dist[tgtIdx] << ".\n";
            cout << "Path: ";
            printPath(parent, tgtIdx, vertices);
            cout << "\n";
        }

        delete[] dist;
        delete[] parent;
    }
};

// ===========================
// SORTING ALGORITHMS (UTILITIES)
// ===========================
namespace SortSearchUtils {

    // QuickSort for int array
    void quickSort(int arr[], int left, int right) {
        if (left >= right) return;
        int pivot = arr[(left + right) / 2];
        int i = left, j = right;
        while (i <= j) {
            while (arr[i] < pivot) i++;
            while (arr[j] > pivot) j--;
            if (i <= j) {
                int temp = arr[i];
                arr[i] = arr[j];
                arr[j] = temp;
                i++; j--;
            }
        }
        if (left < j) quickSort(arr, left, j);
        if (i < right) quickSort(arr, i, right);
    }

    // MergeSort helpers
    void merge(int arr[], int l, int m, int r) {
        int n1 = m - l + 1;
        int n2 = r - m;
        int* L = new int[n1];
        int* R = new int[n2];
        for (int i = 0; i < n1; ++i) L[i] = arr[l + i];
        for (int j = 0; j < n2; ++j) R[j] = arr[m + 1 + j];
        int i = 0, j = 0, k = l;
        while (i < n1 && j < n2) {
            if (L[i] <= R[j]) {
                arr[k++] = L[i++];
            }
            else {
                arr[k++] = R[j++];
            }
        }
        while (i < n1) arr[k++] = L[i++];
        while (j < n2) arr[k++] = R[j++];
        delete[] L;
        delete[] R;
    }

    void mergeSort(int arr[], int l, int r) {
        if (l >= r) return;
        int m = l + (r - l) / 2;
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);
        merge(arr, l, m, r);
    }

    // Linear Search: return index or -1
    int linearSearch(int arr[], int size, int target) {
        for (int i = 0; i < size; ++i) {
            if (arr[i] == target) return i;
        }
        return -1;
    }

    // Binary Search: array must be sorted; return index or -1
    int binarySearch(int arr[], int left, int right, int target) {
        while (left <= right) {
            int mid = left + (right - left) / 2;
            if (arr[mid] == target) return mid;
            if (arr[mid] < target) left = mid + 1;
            else right = mid - 1;
        }
        return -1;
    }
}

// ===========================
// MAIN MENU AND CLI
// ===========================
void displayMenu() {
    cout << "\n========= Shortest Path Finder Menu =========\n";
    cout << "1.  Add Vertex\n";
    cout << "2.  Remove Vertex\n";
    cout << "3.  Add Edge\n";
    cout << "4.  Remove Edge\n";
    cout << "5.  Update Edge Weight\n";
    cout << "6.  Display Graph\n";
    cout << "7.  Find Shortest Path (Dijkstra)\n";
    cout << "8.  Sort Vertex IDs (QuickSort)\n";
    cout << "9.  Sort Vertex IDs (MergeSort)\n";
    cout << "10. Search for Vertex (Linear Search)\n";
    cout << "11. Search for Vertex (Binary Search)\n";
    cout << "12. Exit\n";
    cout << "==============================================\n";
    cout << "Enter your choice: ";
}

int main() {
   // ios::sync_with_stdio(false);
    // cin.tie(nullptr);

    Graph graph;
    ShortestPathFinder spFinder(&graph);

    int choice = 0;

    while (true) {
        displayMenu();
        cin >> choice;
        if (cin.fail()) {
            cin.clear();
            cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            cout << "Invalid input. Please enter a number between 1 and 12.\n";
            continue;
        }

        if (choice == 1) {
            // Add Vertex
            cout << "Enter new vertex ID (integer): ";
            int vid;
            cin >> vid;
            if (graph.addVertex(vid)) {
                cout << "Vertex " << vid << " added.\n";
            }
            else {
                cout << "Vertex " << vid << " already exists.\n";
            }
        }
        else if (choice == 2) {
            // Remove Vertex
            cout << "Enter vertex ID to remove: ";
            int vid;
            cin >> vid;
            if (graph.removeVertex(vid)) {
                cout << "Vertex " << vid << " removed (and associated edges).\n";
            }
            else {
                cout << "Vertex " << vid << " does not exist.\n";
            }
        }
        else if (choice == 3) {
            // Add Edge
            cout << "Enter source vertex ID: ";
            int src;
            cin >> src;
            cout << "Enter destination vertex ID: ";
            int dst;
            cin >> dst;
            cout << "Enter weight of edge: ";
            int w;
            cin >> w;
            if (graph.addEdge(src, dst, w)) {
                cout << "Edge added/updated: " << src << " -> " << dst << " with weight " << w << ".\n";
            }
            else {
                cout << "Failed to add edge. One or both vertices do not exist.\n";
            }
        }
        else if (choice == 4) {
            // Remove Edge
            cout << "Enter source vertex ID: ";
            int src;
            cin >> src;
            cout << "Enter destination vertex ID: ";
            int dst;
            cin >> dst;
            if (graph.removeEdge(src, dst)) {
                cout << "Edge removed: " << src << " -> " << dst << ".\n";
            }
            else {
                cout << "Failed to remove edge. Edge may not exist or vertices invalid.\n";
            }
        }
        else if (choice == 5) {
            // Update Edge Weight
            cout << "Enter source vertex ID: ";
            int src;
            cin >> src;
            cout << "Enter destination vertex ID: ";
            int dst;
            cin >> dst;
            cout << "Enter new weight for edge: ";
            int nw;
            cin >> nw;
            if (graph.updateEdge(src, dst, nw)) {
                cout << "Edge weight updated: " << src << " -> " << dst << " new weight " << nw << ".\n";
            }
            else {
                cout << "Failed to update. Edge may not exist or vertices invalid.\n";
            }
        }
        else if (choice == 6) {
            // Display Graph
            graph.displayGraph();
        }
        else if (choice == 7) {
            // Find Shortest Path using Dijkstra
            cout << "Enter source vertex ID: ";
            int src;
            cin >> src;
            cout << "Enter target vertex ID: ";
            int tgt;
            cin >> tgt;
            spFinder.dijkstra(src, tgt);
        }
        else if (choice == 8) {
            // Sort Vertex IDs using QuickSort
            DynamicVertexArray& vArrQ = graph.getVerticesArray();
            int nQ = vArrQ.getSize();
            if (nQ == 0) {
                cout << "No vertices to sort.\n";
            }
            else {
                // Extract IDs into simple array
                int* ids = new int[nQ];
                for (int i = 0; i < nQ; ++i) ids[i] = vArrQ.get(i)->getID();
                SortSearchUtils::quickSort(ids, 0, nQ - 1);
                cout << "Vertex IDs sorted (QuickSort): [ ";
                for (int i = 0; i < nQ; ++i) {
                    cout << ids[i] << (i < nQ - 1 ? ", " : " ");
                }
                cout << "]\n";
                delete[] ids;
            }
        }
        else if (choice == 9) {
            // Sort Vertex IDs using MergeSort
            DynamicVertexArray& vArrM = graph.getVerticesArray();
            int nM = vArrM.getSize();
            if (nM == 0) {
                cout << "No vertices to sort.\n";
            }
            else {
                int* ids = new int[nM];
                for (int i = 0; i < nM; ++i) ids[i] = vArrM.get(i)->getID();
                SortSearchUtils::mergeSort(ids, 0, nM - 1);
                cout << "Vertex IDs sorted (MergeSort): [ ";
                for (int i = 0; i < nM; ++i) {
                    cout << ids[i] << (i < nM - 1 ? ", " : " ");
                }
                cout << "]\n";
                delete[] ids;
            }
        }
        else if (choice == 10) {
            // Search for Vertex (Linear)
            DynamicVertexArray& vArrL = graph.getVerticesArray();
            int nL = vArrL.getSize();
            if (nL == 0) {
                cout << "No vertices present.\n";
            }
            else {
                int* ids = new int[nL];
                for (int i = 0; i < nL; ++i) ids[i] = vArrL.get(i)->getID();
                cout << "Enter vertex ID to search (Linear Search): ";
                int target;
                cin >> target;
                int idx = SortSearchUtils::linearSearch(ids, nL, target);
                if (idx == -1) {
                    cout << "Vertex ID " << target << " not found (Linear Search).\n";
                }
                else {
                    cout << "Vertex ID " << target << " found at array index " << idx << ".\n";
                }
                delete[] ids;
            }
        }
        else if (choice == 11) {
            // Search for Vertex (Binary) - first sort, then search
            DynamicVertexArray& vArrB = graph.getVerticesArray();
            int nB = vArrB.getSize();
            if (nB == 0) {
                cout << "No vertices present.\n";
            }
            else {
                // Extract and sort IDs
                int* ids = new int[nB];
                for (int i = 0; i < nB; ++i) ids[i] = vArrB.get(i)->getID();
                SortSearchUtils::quickSort(ids, 0, nB - 1);
                cout << "Enter vertex ID to search (Binary Search): ";
                int target;
                cin >> target;
                int idxB = SortSearchUtils::binarySearch(ids, 0, nB - 1, target);
                if (idxB == -1) {
                    cout << "Vertex ID " << target << " not found (Binary Search).\n";
                }
                else {
                    cout << "Vertex ID " << target << " found at sorted-array index " << idxB << ".\n";
                }
                delete[] ids;
            }
        }
        else if (choice == 12) {
            cout << "Exiting application. Goodbye!\n";
            break;
        }
        else {
            cout << "Invalid choice. Please select between 1 and 12.\n";
        }
    }

    return 0;
}