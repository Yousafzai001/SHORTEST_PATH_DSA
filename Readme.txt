 Dynamic Graph & Shortest Path Finder (C++ Console App)
📌 Overview
This is a C++ console-based application that allows users to create and manipulate a dynamic graph with full support for:

Vertex and edge operations

Dijkstra’s shortest path algorithm

QuickSort and MergeSort on vertex IDs

Linear and Binary search

Custom dynamic array implementations

It features an interactive menu-driven interface for real-time graph construction and analysis.

🚀 Features
✅ Add / Remove vertices

✅ Add / Remove / Update edges

✅ Display graph (adjacency list)

✅ Compute shortest path using Dijkstra's algorithm

✅ Sort vertex IDs using QuickSort and MergeSort

✅ Search vertex using Linear and Binary Search

⚙️ How It Works
The project includes:

Custom dynamic arrays for edges and vertices

A graph class with adjacency list representation

A min-heap-based Dijkstra algorithm

Built-in sorting & searching utilities

🧪 Menu Options
Option	Action
1	Add Vertex
2	Remove Vertex
3	Add Edge
4	Remove Edge
5	Update Edge Weight
6	Display Graph
7	Find Shortest Path (Dijkstra)
8	Sort Vertex IDs (QuickSort)
9	Sort Vertex IDs (MergeSort)
10	Search Vertex (Linear Search)
11	Search Vertex (Binary Search)
12	Exit Application

🛠️ How to Compile
Make sure you have a C++ compiler (like g++) installed.

bash
Copy
Edit
g++ PROJECT.cpp -o graph_app
./graph_app
👥 Team Members
Mueez Waqar – 231672

Faizan Mehdi – 231257

Abu Huraira Yousafzai – 231257

📚 Algorithms Used
Dijkstra's Algorithm (Min-Heap based)

QuickSort / MergeSort for Vertex IDs

Linear and Binary Search for Vertex Lookup

📈 Future Improvements
Add file-based input/output support

Implement more graph algorithms (e.g., Bellman-Ford, Floyd-Warshall)

GUI-based frontend using Qt or SFML