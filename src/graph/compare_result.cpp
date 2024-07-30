#include "StringTokenizer.h"
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

using namespace std;

template <typename V_T, typename E_T> struct edge_counter {
  typedef pair<pair<V_T, V_T>, E_T> EDGE_T;
  typedef typename map<EDGE_T, unsigned int>::iterator IT;
  typedef typename map<EDGE_T, unsigned int>::const_iterator CIT;

  map<EDGE_T, unsigned int> _counter;

  unsigned int get_count(V_T src_l, V_T dest_l, E_T edge_l) {
    EDGE_T e;
    if (src_l < dest_l)
      e = make_pair(make_pair(src_l, dest_l), edge_l);
    else
      e = make_pair(make_pair(dest_l, src_l), edge_l);
    CIT cit = _counter.find(e);
    if (cit != _counter.end())
      return cit->second;
    else
      return 0;
  }
  void insert(V_T src_l, V_T dest_l, E_T edge_l) {
    EDGE_T e;
    if (src_l < dest_l)
      e = make_pair(make_pair(src_l, dest_l), edge_l);
    else
      e = make_pair(make_pair(dest_l, src_l), edge_l);
    IT it = _counter.find(e);
    if (it == _counter.end())
      _counter.insert(make_pair(e, 1));
    else
      it->second++;
  }
  void print() const {
    CIT cit = _counter.begin();
    for (; cit != _counter.end(); cit++)

      cout << "(" << cit->first.first.first << " " << cit->first.second << " "
           << cit->first.first.second << "):" << cit->second << endl;
  }

  double operator-(const edge_counter &other) const {
    const map<EDGE_T, unsigned int> &other_cnt = other._counter;
    CIT cit;
  }

  bool operator<(const edge_counter &other) const {
    const map<EDGE_T, unsigned int> &other_cnt = other._counter;
    return (_counter < other_cnt);
  }
  bool operator==(const edge_counter &other) const {
    const map<EDGE_T, unsigned int> &other_cnt = other._counter;
    return !(_counter < other_cnt || other_cnt < _counter);
  }
};

int load_spin_result(const char *node_file, const char *edge_file,
                     vector<edge_counter<int, int> *> &all_graphs) {

  ifstream nodefile(node_file, ios::in);
  ifstream edgefile(edge_file, ios::in);

  vector<vector<int> *> vertex_label;
  string one_line;
  int linesize, last_graph_id = -1;
  vector<int> *a_graph = 0;

  if (!nodefile) {
    cout << "Input file:" << node_file << " could not open to read from!"
         << endl;
    exit(1);
  }

  if (!edgefile) {
    cout << "Input file:" << edge_file << " could not open to read from!"
         << endl;
    exit(1);
  }

  //////////// ///////////// Reading Node file
  ///////////////////////////////////////
  while (!nodefile.eof()) {
    std::getline(nodefile, one_line);
    linesize = one_line.length();
    if (linesize < 5)
      break;
    StringTokenizer strtok = StringTokenizer(one_line, " ");
    int cnt = strtok.countTokens();
    assert(cnt == 4);
    strtok.nextToken(); // skiping the word "node"
    int graph_id = strtok.nextIntToken();
    if (graph_id != last_graph_id) {
      if (a_graph)
        vertex_label.push_back(a_graph);
      a_graph = new vector<int>();
      last_graph_id = graph_id;
    }
    strtok.nextToken(); // skipping the vertex id, as it is always ascending
    int v_label = strtok.nextIntToken();
    a_graph->push_back(v_label);
  }
  vertex_label.push_back(a_graph);
  nodefile.close();

  //////////// ///////////// Reading Edge file
  ///////////////////////////////////////
  last_graph_id = -1;
  edge_counter<int, int> *one_edc = 0;

  while (!edgefile.eof()) {
    std::getline(edgefile, one_line);
    linesize = one_line.length();
    if (linesize < 5)
      break;
    StringTokenizer strtok = StringTokenizer(one_line, " ");
    int cnt = strtok.countTokens();
    assert(cnt == 5);
    strtok.nextToken(); // skiping the word "edge"
    int graph_id = strtok.nextIntToken();
    if (graph_id != last_graph_id) {
      if (one_edc)
        all_graphs.push_back(one_edc);
      one_edc = new edge_counter<int, int>();
      last_graph_id = graph_id;
    }
    int v1 = strtok.nextIntToken();
    int v2 = strtok.nextIntToken();
    int e_label = strtok.nextIntToken();
    vector<int> *&v_vector = vertex_label[graph_id];
    int v1_label = (*v_vector)[v1];
    int v2_label = (*v_vector)[v2];
    one_edc->insert(v1_label, v2_label, e_label);
  }
  all_graphs.push_back(one_edc);
  edgefile.close();

  ////////////////////// Freeing memory
  /////////////////////////////////////////////////
  for (int i = 0; i < vertex_label.size(); i++)
    delete vertex_label[i];
}

int load_our_result(const char *our_file,
                    vector<edge_counter<int, int> *> &all_graphs) {

  ifstream ourfile(our_file, ios::in);

  string one_line;
  int linesize;

  edge_counter<int, int> *one_edc = 0;

  if (!ourfile) {
    cout << "Input file:" << our_file << " could not open to read from!"
         << endl;
    exit(1);
  }

  while (!ourfile.eof()) {
    std::getline(ourfile, one_line);
    linesize = one_line.length();
    if (linesize < 5)
      break;
    StringTokenizer strtok = StringTokenizer(one_line, " ");
    int cnt = strtok.countTokens();
    assert(cnt == 5);
    int vid1 = strtok.nextIntToken(); // finding the vid of edge-1
    int vid2 = strtok.nextIntToken(); // finding the other vid of edge-1
    if (vid1 == 0 && vid2 == 1) {     // new graph starts
      if (one_edc)
        all_graphs.push_back(one_edc);
      one_edc = new edge_counter<int, int>();
    }
    int v1_label = strtok.nextIntToken();
    int e_label = strtok.nextIntToken();
    int v2_label = strtok.nextIntToken();
    one_edc->insert(v1_label, v2_label, e_label);
  }
  all_graphs.push_back(one_edc);
  ourfile.close();
}

int main(int argc, char *argv[]) {
  if (argc < 3) {
    cout << "usage:";
    cout << "\t " << argv[0] << " spin_node_file spin_edge_file our-output-file"
         << endl;
    exit(0);
  } else {
    vector<edge_counter<int, int> *> all_graphs;
    vector<edge_counter<int, int> *> all_graphs2;
    load_spin_result(argv[1], argv[2], all_graphs);
    load_our_result(argv[3], all_graphs2);
    for (int i = 0; i < all_graphs.size(); i++) {
      all_graphs[i]->print();
      cout << "\n--------------------------------\n";
    }
    cout << "\n***********************************\n";
    for (int i = 0; i < all_graphs2.size(); i++) {
      all_graphs2[i]->print();
      cout << "\n--------------------------------\n";
    }
    ////////// Now comparing the two //////////////////////
    int common = 0;
    for (int i = 0; i < all_graphs.size(); i++) {
      for (int j = 0; j < all_graphs2.size(); j++) {
        edge_counter<int, int> &p = *all_graphs[i];
        edge_counter<int, int> &q = *all_graphs2[j];
        if (p < q) {
          common++;
          break;
        }
      }
    }
    double coverage = (double)common / all_graphs.size();
    cout << setprecision(4) << "Coverage:" << coverage << endl;
  }
} // main
