#include "votable.hpp"

#include <algorithm>
#include <cctype>

namespace sed::votable {

namespace {

std::string to_lower(std::string s) {
  for (char& c : s) c = (char)std::tolower((unsigned char)c);
  return s;
}

// Decode the XML entities that appear in VOTable text/cell content.
std::string xml_unescape(const std::string& s) {
  std::string out;
  out.reserve(s.size());
  for (size_t i = 0; i < s.size();) {
    if (s[i] == '&') {
      if (s.compare(i, 5, "&amp;") == 0) { out += '&'; i += 5; continue; }
      if (s.compare(i, 4, "&lt;") == 0) { out += '<'; i += 4; continue; }
      if (s.compare(i, 4, "&gt;") == 0) { out += '>'; i += 4; continue; }
      if (s.compare(i, 6, "&quot;") == 0) { out += '"'; i += 6; continue; }
      if (s.compare(i, 6, "&apos;") == 0) { out += '\''; i += 6; continue; }
    }
    out += s[i++];
  }
  return out;
}

std::string trim(const std::string& s) {
  size_t a = 0, b = s.size();
  while (a < b && std::isspace((unsigned char)s[a])) ++a;
  while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
  return s.substr(a, b - a);
}

// Extract attribute value: attr="..." within a tag substring. "" if absent.
std::string attr(const std::string& tag, const std::string& key) {
  std::string needle = key + "=\"";
  size_t p = tag.find(needle);
  if (p == std::string::npos) return "";
  p += needle.size();
  size_t e = tag.find('"', p);
  if (e == std::string::npos) return "";
  return tag.substr(p, e - p);
}

// Parse FIELD names within a table block. Matches "<FIELD" followed by a
// space/newline (excludes <FIELDref>), reads its name attribute.
std::vector<std::string> parse_fields(const std::string& block) {
  std::vector<std::string> fields;
  size_t pos = 0;
  while (true) {
    size_t p = block.find("<FIELD", pos);
    if (p == std::string::npos) break;
    char after = block[p + 6];
    if (after != ' ' && after != '\t' && after != '\n' && after != '\r') {
      pos = p + 6;  // e.g. <FIELDref>
      continue;
    }
    size_t e = block.find('>', p);
    if (e == std::string::npos) break;
    std::string tag = block.substr(p, e - p);
    fields.push_back(xml_unescape(attr(tag, "name")));
    pos = e + 1;
  }
  return fields;
}

// Parse TABLEDATA rows into aligned string cells. Handles <TD>v</TD> and the
// self-closing <TD/> (empty). Cells beyond the row are left empty.
std::vector<std::vector<std::string>> parse_data(const std::string& block,
                                                 size_t n_fields) {
  std::vector<std::vector<std::string>> rows;
  size_t td_start = block.find("<TABLEDATA>");
  size_t td_end = block.find("</TABLEDATA>");
  if (td_start == std::string::npos || td_end == std::string::npos) return rows;
  td_start += 11;
  std::string data = block.substr(td_start, td_end - td_start);

  size_t pos = 0;
  while (true) {
    size_t tr = data.find("<TR", pos);
    if (tr == std::string::npos) break;
    size_t tr_open_end = data.find('>', tr);
    if (tr_open_end == std::string::npos) break;
    size_t tr_end = data.find("</TR>", tr_open_end);
    if (tr_end == std::string::npos) break;
    std::string row_block = data.substr(tr_open_end + 1, tr_end - tr_open_end - 1);
    pos = tr_end + 5;

    std::vector<std::string> cells(n_fields, "");
    size_t cpos = 0, j = 0;
    while (j < n_fields) {
      size_t td = row_block.find("<TD", cpos);
      if (td == std::string::npos) break;
      size_t tag_end = row_block.find('>', td);
      if (tag_end == std::string::npos) break;
      bool self_closing = (tag_end > td && row_block[tag_end - 1] == '/');
      if (self_closing) {
        cells[j] = "";
        cpos = tag_end + 1;
        ++j;
        continue;
      }
      size_t close = row_block.find("</TD>", tag_end + 1);
      if (close == std::string::npos) break;
      cells[j] = xml_unescape(trim(row_block.substr(tag_end + 1, close - tag_end - 1)));
      cpos = close + 5;
      ++j;
    }
    rows.push_back(std::move(cells));
  }
  return rows;
}

}  // namespace

int Table::field_index(const std::string& name) const {
  for (size_t i = 0; i < fields.size(); ++i)
    if (fields[i] == name) return (int)i;
  return field_index_ci(name);
}

int Table::field_index_ci(const std::string& name) const {
  std::string low = to_lower(name);
  for (size_t i = 0; i < fields.size(); ++i)
    if (to_lower(fields[i]) == low) return (int)i;
  return -1;
}

const Table* Document::find(const std::string& catalogue) const {
  std::string id_form = catalogue;
  std::replace(id_form.begin(), id_form.end(), '/', '_');
  for (const auto& t : tables)
    if (t.id == id_form || t.name == catalogue) return &t;
  // suffix fallback (some servers prefix table ids)
  for (const auto& t : tables) {
    if (!t.id.empty() && t.id.size() >= id_form.size() &&
        t.id.compare(t.id.size() - id_form.size(), id_form.size(), id_form) == 0)
      return &t;
    if (!t.name.empty() && t.name == catalogue) return &t;
  }
  return nullptr;
}

Document parse(const std::string& content) {
  Document doc;
  size_t pos = 0;
  while (true) {
    size_t p = content.find("<TABLE", pos);
    if (p == std::string::npos) break;
    char after = content[p + 6];
    if (after != ' ' && after != '\t' && after != '\n' && after != '\r' &&
        after != '>') {
      pos = p + 6;  // e.g. <TABLEDATA>
      continue;
    }
    size_t open_end = content.find('>', p);
    if (open_end == std::string::npos) break;
    size_t close = content.find("</TABLE>", open_end);
    if (close == std::string::npos) break;

    std::string open_tag = content.substr(p, open_end - p);
    std::string block = content.substr(p, close + 8 - p);
    pos = close + 8;

    Table t;
    t.id = attr(open_tag, "ID");
    t.name = xml_unescape(attr(open_tag, "name"));
    t.fields = parse_fields(block);
    if (!t.fields.empty()) t.rows = parse_data(block, t.fields.size());
    doc.tables.push_back(std::move(t));
  }
  return doc;
}

}  // namespace sed::votable
