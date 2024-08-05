# ORIGAMI: Representative Orthogonal Graph Mining

Origami uses random walks over the graph partial order to mine a representative sample of the maximal frequent subgraph patterns. It next selects only the orthogonal representative patterns via a local clique-finding method.

**Relevant Publications** 

* [2008-origami:sadm] Vineet Chaoji, Mohammad Al Hasan, Saeed Salem, Jeremy Besson, and Mohammed J. Zaki. ORIGAMI: A Novel and Effective Approach for Mining Representative Orthogonal Graph Patterns. Statistical Analysis and Data Mining, 1(2):67–84, June 2008. doi:10.1002/sam.10004. 

* [2007-origami] Mohammad Hasan, Vineet Chaoji, Saeed Salem, Jeremy
    Besson, and Mohammed J. Zaki. ORIGAMI: mining representative
    orthogonal graph patterns. In 7th IEEE International Conference on
    Data Mining. October 2007.

## How to
```sh
cd test && mkdir build && cd build
cmake ..
make
```
