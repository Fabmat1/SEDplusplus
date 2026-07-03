// Minimal VOTable parser for the query layer. Handles the multi-TABLE VOTable
// returned by VizieR's CGI `viz-bin/votable` endpoint (one TABLE per catalog,
// tagged with ID="II_335_galex_ais") as well as single-TABLE TAP responses.
// String-based and tolerant, mirroring parse_vizier_votable_multi in
// query_vizier_catalogue.sl but parsing every table in one pass.
#pragma once

#include <string>
#include <vector>

namespace sed::votable {

struct Table {
  std::string id;    // TABLE ID attribute ("" if absent)
  std::string name;  // TABLE name attribute ("" if absent)
  std::vector<std::string> fields;               // FIELD names, in order
  std::vector<std::vector<std::string>> rows;    // cells, aligned to fields

  // Column index by exact name, then case-insensitive; -1 if not found.
  int field_index(const std::string& name) const;
  int field_index_ci(const std::string& name) const;
};

struct Document {
  std::vector<Table> tables;

  // Locate a catalog's table: exact ID == cat-with-'/'-replaced-by-'_', else
  // exact name == cat, else ID/name suffix match. nullptr if none.
  const Table* find(const std::string& catalogue) const;
};

// Parse a full VOTable document. Never throws; returns empty on garbage.
Document parse(const std::string& content);

}  // namespace sed::votable
