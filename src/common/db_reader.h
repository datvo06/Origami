/*
 *  Copyright (C) 2005 M.J. Zaki <zaki@cs.rpi.edu> Rensselaer Polytechnic
 * Institute Written by parimi@cs.rpi.edu Updated by chaojv@cs.rpi.edu,
 * alhasan@cs.rpi.edu, salems@cs.rpi.edu
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
 */
/** \file db_reader.h - user defined class to parse input database and retrieve
 * level-1 patterns */
#ifndef _DB_READER_H
#define _DB_READER_H

#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// using namespace std;

// adding all files manually for now, TODO: improve this include system

#include "generic_classes.h"
#include "helper_funs.h"
#include "pat_fam.h" // added later to make .cpp files
/**
 * \brief Database Reader class, to read the input file.
 *
 * This class reads the database, using the tokenizer class and populate the
 * level-1 VAT.
 */
template <typename PATTERN, typename TOKENIZER> class db_reader {
public:
  typedef vat<typename PATTERN::PAT_PROPS, typename PATTERN::MINE_PROPS,
              std::vector>
      VAT; /// vat typedef
  typedef tokenizer<PATTERN, TOKENIZER>
      TKNZ; /// tokenizer class for this pattern-type, its method
            /// parse_next_line() is invoked by db_reader
  typedef map<pair<pair<typename PATTERN::VERTEX_T, typename PATTERN::VERTEX_T>,
                   typename PATTERN::EDGE_T>,
              int>
      FREQ_MAP;

  /** \fn db_reader(const char* infile_name)
   * \brief Constructor
   * \param infile_name Name of the input database (flat) file
   */
  db_reader(const char *infile_name)
      : _in_db(infile_name), filename(std::string(infile_name)) {}

  /** \fn db_reader(const char* infile_name, int mem_size)
   * \brief Constructor_for_gigabase
   * \param infile_name Name of the input database (flat) file.
   * \param mem_size Maximum size of memory vat for gigabase backend.
   */
  db_reader(const char *infile_name, int mem_size) : _in_db(infile_name) {
    filename = std::string(infile_name);
    std::cout << "Filename: " << filename << std::endl;
    _max_mem = mem_size;
  }

  /** \fn ~db_reader()
   * \brief destructor
   */
  ~db_reader() { close(); }

  /** \fn void open(const char* infile_name)
   * \brief Opens the specified input file. This is an alternative to
   * the parameterized constructor
   * \param infile_name Name of the input file
   */
  void open(const char *infile_name) {
    if (is_open())
      _in_db.close();
    _in_db.open(infile_name);
  }

  /** void close()
   * \brief Closes the file associated with this class
   */
  void close() { _in_db.close(); }

  /** bool is_open
   * \brief returns whether file associated with this object is open
   */
  bool is_open() { return _in_db.is_open(); }

  /** void get_length_one(pat_fam<PATTERN>& freq_pats, vat_db<PATTERN, VAT,0>&
   * vat_hmap, int minsup) \brief obtain length one frequent patterns in sorted
   * order, and populate vat_db with their vats \param freq_pats Pattern Family
   * which is populated with the frequent patterns \param vat_hmap The hashmap
   * used to store pattern-to-VAT mappings \param minsup Minimum support
   * threshold
   */
  template <class SM_T>
  void get_length_one(pat_fam<PATTERN> &freq_pats,
                      storage_manager<PATTERN, VAT, SM_T> &vat_hmap,
                      const int &minsup, FREQ_MAP &fm) {

    int tid;
    VAT *ivat;
    typename pat_fam<PATTERN>::IT pf_it;

    if (!is_open()) {
      // stream not open
      std::cerr << "db_reader: file stream not open in get_length_one()"
                << std::endl;
      return;
    }

    tid = tknz.parse_next_trans(_in_db, freq_pats, vat_hmap, fm);
    int i = 0; // i keep track of total transaction read
    while (tid != -1) {
      // cout << "Read " << i << " transaction\n";
      i++;
      tid = tknz.parse_next_trans(_in_db, freq_pats, vat_hmap, fm);
    }
    _trans_cnt = i;

    // fill in support of level-1, discarding infrequent ones
    for (pf_it = freq_pats.begin(); pf_it != freq_pats.end(); ++pf_it) {
      if (!(ivat = vat_hmap.get_vat(*pf_it))) {
        std::cerr << "db_reader.get_length_one: VAT not found for " << *pf_it
                  << std::endl;
        return;
      }

      // cout << "LEVEL 1 " << *ivat << endl;

      if ((ivat->size()) >= minsup)
        (*pf_it)->set_sup(make_pair(ivat->size(), 0));
      else {
        // Delete the pattern and the vat.
        vat_hmap.delete_vat(*pf_it);
        delete (*pf_it);

        freq_pats.erase(pf_it);
        pf_it--;
      }
    } // end for

    // sort level-1 patterns
    // typename pat_fam<PATTERN>::IT b=freq_pats.begin(), e=freq_pats.end();
    // sort(b, e, less_than<PATTERN>());

  } // end get_length_one()
  static void print_edge_freq_map(const FREQ_MAP &fm) {

    cout << "Inside print_edge_freq_map." << endl;
    typename FREQ_MAP::const_iterator cit = fm.begin();
    for (; cit != fm.end(); cit++)
      cout << "(" << cit->first.first.first << ", " << cit->first.second << ", "
           << cit->first.first.second << ") -->" << cit->second << endl;
  }

  unsigned int get_transaction_count() const { return _trans_cnt; }

private:
  std::ifstream _in_db;
  std::string filename; // Holds the file name of the dataset
  unsigned long _max_mem;
  TKNZ tknz; // An object of Tokenizer class
  unsigned int _trans_cnt;
}; // end class db_reader<itemset>

#endif
