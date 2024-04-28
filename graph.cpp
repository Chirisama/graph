#include "graph.h"
#include <algorithm>
#include <climits>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <unordered_set>
#include <utility>
#include <vector>

// NOLINTNEXTLINE
using namespace std;

/**
 * Constructor, empty graph where directionalEdges defaults to true
 */
Graph::Graph(bool directionalEdges) {
  this->directionalEdges = directionalEdges;
}

/**
 * Destructor for graph
 */
Graph::~Graph() = default;

/**
 * This method returns the number of vertices.
 * No input parameters.
 * @return total number of vertices
 */
int Graph::verticesSize() const {
  int count = 0;
  for (auto vertice = graph.begin(); vertice != graph.end(); vertice++) {
    count++;
  }
  return count;
}

/**
 * This method returns the number of edges.
 * No input parameters.
 * @return total number of edges
 */
int Graph::edgesSize() const {
  int count = 0;
  for (auto const &i : graph) {
    auto const &v = i.second;
    for (int j = 0; j < v.size(); j++) {
      count++;
    }
  }
  return count;
}

/**
 * This method returns the number of edges from given vertex.
 * input parameter string label which is the string representing the vertex of a
 * graph.
 * @return number of edges from given vertex, -1 if vertex not found
 */
int Graph::vertexDegree(const string &label) const {
  if (graph.find(label) != graph.end()) {
    vector<pair<string, int>> v = graph.find(label)->second;
    int count = 0;
    for (int i = 0; i < v.size(); i++) {
      count++;
    }
    return count;
  }
  return -1;
}

/**
* This method returns true if the vertex is added.
* This method has one input parameter which is string label which represents the
vertex of the graph.
* @return true if vertex added, false if it already is in the graph
*/
bool Graph::add(const string &label) {
  if (graph.find(label) == graph.end()) {
    vector<pair<string, int>> v;
    graph.insert(make_pair(label, v));
    return true;
  }
  return false;
}

/**
* This method returns true if the vertex is added.
* This method has one input parameter which is string label which represents the
vertex of the graph.
* @return true if vertex already in graph
*/
bool Graph::contains(const string &label) const {
  return (graph.find(label) != graph.end());
}

/**
* This method returns the edges as a string.
* This method has one input parameter which is string label which represents the
vertex of the graph.
* @return string representing edges and weights, "" if vertex not found
*/
// A-3->B, A-5->C should return B(3),C(5)
string Graph::getEdgesAsString(const string &label) const {
  string res;
  if (contains(label)) {
    vector<pair<string, int>> v = graph.find(label)->second;
    priority_queue<pair<string, int>, vector<pair<string, int>>,
                   greater<pair<string, int>>>
        q;
    // NOLINTNEXTLINE, the solution clang-tidy suggests crashes.
    for (int i = 0; i < v.size(); i++) {
      pair<string, int> k = make_pair(v[i].first, v[i].second);

      q.push(k);
    }

    while (!q.empty()) {
      string first = q.top().first;
      string second = to_string(q.top().second);
      res += first += "(" + second + "),";
      q.pop();
    }
    if (!res.empty()) {
      res.pop_back();
    }
  }
  return res;
}

/**
* This method returns true if successfully connected.
* This method has three input parameters: strings from and to where the
connection is happenning and weight of the vertices.
* @return true if successfully connected
*/
// NOLINTNEXTLINE, Intend to use recursive function.
bool Graph::connect(const string &from, const string &to, int weight) {
  if (from == to)
    return false;

  if (graph.find(from) == graph.end()) {
    add(from);
    return connect(from, to, weight);
  }

  if (graph.find(to) == graph.end()) {
    add(to);
    return connect(from, to, weight);
  }
  if (!directionalEdges) {
    vector<pair<string, int>> j = graph.find(to)->second;
    // NOLINTNEXTLINE, the solution clang-tidy suggests crashes.
    for (int i = 0; i < j.size(); i++) {
      if (j[i].first == from) {
        return false;
      }
    }
    graph.find(to)->second.emplace_back(from, weight);
  }
  vector<pair<string, int>> v = graph.find(from)->second;
  // NOLINTNEXTLINE, the solution clang-tidy suggests crashes.
  for (int i = 0; i < v.size(); i++) {
    if (v[i].first == to) {
      return false;
    }
  }
  graph.find(from)->second.emplace_back(to, weight);
  return true;
}

/**
* This method returns true if successfully disconnected.
* This method has three input parameters: strings from and to where the
connection is disconnected.
* @return true if successfully disconnected
*/
bool Graph::disconnect(const string &from, const string &to) {
  // Check if both vertices exist in the graph
  if (graph.find(from) == graph.end() || graph.find(to) == graph.end()) {
    // One or both vertices not found
    return false;
  }

  // Find and remove the edge from 'from' to 'to'
  if (!directionalEdges) {
    auto &edges = graph[to];
    for (auto node = edges.begin(); node != edges.end(); node++) {
      if (node->first == from) {
        edges.erase(node);
      }
    }
  }

  auto &edges = graph[from];
  for (auto node = edges.begin(); node != edges.end(); node++) {
    if (node->first == to) {
      edges.erase(node);
      return true;
    }
  }

  // Edge not found
  return false;
}

/**
* This method performs depth-first traversal starting from given startLabel
* This method has two input parameters: string startLabel which denotes the
starting point of the graph or the first vertex and void visit which is a
method.
*/
void Graph::dfs(const string &startLabel, void visit(const string &label)) {
  if (contains(startLabel)) {
    vector<string> res;
    map<string, bool> is_passed;
    res.push_back(startLabel);
    is_passed.insert(make_pair(startLabel, true));
    stack<string> stk;
    vector<pair<string, int>> v = graph[startLabel];
    priority_queue<string> q;
    // NOLINTNEXTLINE, the solution clang-tidy suggests crashes.
    for (int i = 0; i < v.size(); i++) {
      q.push(v[i].first);
    }
    while (!q.empty()) {
      stk.push(q.top());
      is_passed.insert(make_pair(q.top(), true));
      q.pop();
    }
    while (!stk.empty()) {
      string &top = stk.top();
      v = graph[top];
      stk.pop();
      res.push_back(top);
      is_passed.insert(make_pair(top, true));
      // NOLINTNEXTLINE, the solution clang-tidy suggests crashes.
      for (int i = 0; i < v.size(); i++) {
        if (is_passed.find(v[i].first) == is_passed.end()) {
          q.push(v[i].first);
        }
      }
      while (!q.empty()) {
        stk.push(q.top());
        q.pop();
      }
    }
    // NOLINTNEXTLINE, the solution clang-tidy suggests crashes.
    for (int i = 0; i < res.size(); i++) {
      visit(res[i]);
    }
  }
}

/**
* This method performs breadth-first traversal starting from given startLabel
* This method has two input parameters: string startLabel which denotes the
starting point of the graph or the first vertex and void visit which is a
method.
*/
void Graph::bfs(const string &startLabel, void visit(const string &label)) {
  queue<string> q;
  priority_queue<string, vector<string>, greater<string>> pq;
  vector<pair<string, int>> neighbor = graph[startLabel];
  unordered_set<string> visited;
  q.push(startLabel);
  visited.insert(startLabel);
  while (!q.empty()) {
    string current_label = q.front();
    q.pop();
    visit(current_label);
    vector<pair<string, int>> neighbor = graph[current_label]; // here
    for (int i = neighbor.size() - 1; i >= 0; i--) {
      string &neighbor_label = neighbor[i].first;
      // If the neighbor has not been visited yet, enqueue and mark as visited
      if (visited.find(neighbor_label) == visited.end()) {
        pq.push(neighbor_label);
      }
    }
    while (!pq.empty()) {
      q.push(pq.top());
      visited.insert(pq.top());
      pq.pop();
    }
  }
}

/**
* This method performs dijkstra traversal starting from given startLabel, store
the weights in a map and store the previous label in a map.
* This method has one input parameter: string startLabel which denotes the
starting point of the graph or the first vertex.
*/
pair<map<string, int>, map<string, string>>
Graph::dijkstra(const string &startLabel) const {
  map<string, int> weights;
  map<string, string> previous;
  // TODO(student) Your code here
  string str = startLabel;
  if (!contains(startLabel))
    return make_pair(weights, previous);
  map<string, int> stack;
  vector<pair<string, int>> neighbor = graph.find(str)->second;
  int prev_cost = 0;
  map<string, int>::iterator it;
  while (!stack.empty() || prev_cost == 0) {

    neighbor = graph.find(str)->second;
    // NOLINTNEXTLINE, the solution clang-tidy suggests crashes.
    for (int i = 0; i < neighbor.size(); i++) {
      if (neighbor[i].first == startLabel ||
          weights.find(neighbor[i].first) != weights.end()) {
        continue;
      }
      if (stack.find(neighbor[i].first) == stack.end() &&
          neighbor[i].first != startLabel &&
          previous.find(neighbor[i].first) == previous.end()) {
        stack.insert(
            make_pair(neighbor[i].first, neighbor[i].second + prev_cost));
        previous.insert(make_pair(neighbor[i].first, str));
      } else if (directionalEdges) {
        if (stack.find(neighbor[i].first)->second >
            neighbor[i].second + prev_cost) {
          stack.find(neighbor[i].first)->second =
              neighbor[i].second + prev_cost;
          previous.find(neighbor[i].first)->second = str;
        }
      } else {
        if (previous.find(neighbor[i].first)->first != str) {
          if (stack.find(neighbor[i].first)->second >
              (neighbor[i].second + prev_cost)) {
            stack.find(neighbor[i].first)->second =
                (neighbor[i].second + prev_cost);
            previous.find(neighbor[i].first)->second = str;
          }
        }
      }
    }
    for (it = stack.begin(); it != stack.end(); it++) {
      if (it == stack.begin()) {
        str = it->first;
      }
      if (it->second < stack.find(str)->second) {
        str = it->first;
      }
    }

    weights.insert(make_pair(str, stack.find(str)->second));
    prev_cost = weights.find(str)->second;
    stack.erase(stack.find(str));
  }
  return make_pair(weights, previous);
}

/**
* This method performs minimum spanning tree using Prim's algorithm.
* This method has two input parameters: string startLabel which denotes the
starting point of the graph or the first vertex and the visit method.
* @return returns an integer which is the total weight.
*/
int Graph::mstPrim(const string &startLabel,
                   void visit(const string &from, const string &to,
                              int weight)) const {
  if (!contains(startLabel))
    return -1;
  priority_queue<pair<int, string>, vector<pair<int, string>>,
                 greater<pair<int, string>>>
      prim;
  set<string> mst_set;
  map<string, int> key;
  map<string, string> parent;

  for (const auto &vertex : graph) {
    key[vertex.first] = INT_MAX;
    parent[vertex.first] = "";
  }
  key[startLabel] = 0;
  prim.emplace(0, startLabel);

  int total_weight = 0;

  while (!prim.empty()) {
    string u = prim.top().second;
    prim.pop();

    if (mst_set.find(u) == mst_set.end()) {
      mst_set.insert(u);
      if (!parent[u].empty()) {
        visit(parent[u], u, key[u]);
        total_weight += key[u];
      }

      for (const auto &edge : graph.at(u)) {
        string v = edge.first;
        int weight = edge.second;
        if (mst_set.find(v) == mst_set.end() && weight < key[v]) {
          key[v] = weight;
          parent[v] = u;
          prim.emplace(key[v], v);
        }
      }
    }
  }

  return total_weight;
}

/**
 * This method reads a text file and create the graph
 * This method has one input parameter: string filename which denotes the
 * @return true if graph is created
 */
bool Graph::readFile(const string &filename) {
  ifstream myfile(filename);
  if (!myfile.is_open()) {
    cerr << "Failed to open " << filename << endl;
    return false;
  }
  int edges = 0;
  int weight = 0;
  string from_vertex;
  string to_vertex;
  myfile >> edges;
  for (int i = 0; i < edges; ++i) {
    myfile >> from_vertex >> to_vertex >> weight;
    connect(from_vertex, to_vertex, weight);
  }
  myfile.close();
  return true;
}