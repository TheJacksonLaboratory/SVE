#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

struct Coordinate {
  string chr;
  unsigned int begin;
  unsigned int end;
};

template <class T>
bool OpenFile(const string& filename, T& file) {

  file.open(filename.c_str());
  if (!file.good()) {
    cerr << "ERROR: " << filename << " cannot open." << endl;
    return false;
  }
  return true;
}

void LoadCoordinateVector(ifstream& file, vector<Coordinate>& cor_map) {
  Coordinate temp;
  string temp_s;
  while (!file.eof()) {
    getline(file, temp_s);
    if (file.eof()) break;
    if (temp_s[0] == '#') continue; // skip the comment line
    stringstream ss(temp_s);
    ss >> temp.chr >> temp.begin >> temp.end;
    cor_map.push_back(temp);
  }
}

int SearchCoordinate(const Coordinate& pos, const vector<Coordinate>& cor_map) {
  for (unsigned int i = 0; i < cor_map.size(); ++i) {
    Coordinate pos1 = cor_map[i];
    if (pos.chr == pos1.chr && pos.begin == pos1.begin && pos.end == pos1.end) {
      return i;
    }
  }
  return -1;
}

void RemoveUnmapped(const vector<Coordinate>& unmapped_map, vector<Coordinate>& cor_map){
  for (const auto pos : unmapped_map) {
    int unmapped_id = SearchCoordinate(pos, cor_map);
    if (unmapped_id != -1) cor_map.erase(cor_map.begin() + unmapped_id);
  }
}

bool LiftOverFna(const vector<Coordinate>& cor1_map, const vector<Coordinate>& cor2_map, 
                 const vector<Coordinate>& unmapped_map, 
                 const string& bplib_prefix, const string& output_prefix, const string& fnaOrIns) {
  ifstream file;
  ofstream output;
  if (fnaOrIns == "fna") {
    if (!OpenFile(bplib_prefix+".fna", file)) return false;
    if (!OpenFile(output_prefix+".fna", output)) return false;
  } else if (fnaOrIns == "ins") {
    if (!OpenFile(bplib_prefix+".ins", file)) return false;
    if (!OpenFile(output_prefix+".ins", output)) return false;
  } else {
    return false;
  }

  while (!file.eof()) {
    string temp_s;
    getline(file, temp_s);
    if (file.eof()) break;
    if (temp_s[0] == '>') {// Name
      temp_s.erase(0, 1); // remove '>'
      size_t pos = temp_s.find(':');
      Coordinate temp_cor;
      temp_cor.chr = "chr" + temp_s.substr(0, pos);
      size_t pos1 = temp_s.find('-');
      temp_cor.begin = atoi(temp_s.substr(pos + 1, pos1 - pos).c_str());
      pos = temp_s.find(':', pos1 + 1);
      temp_cor.end = atoi(temp_s.substr(pos1 + 1, pos - pos1).c_str());
      int id = SearchCoordinate(temp_cor, cor1_map);
      if (id == -1) {
        int unmapped_id = SearchCoordinate(temp_cor, unmapped_map);
        getline(file, temp_s); // get rid of the seq as well
        if (unmapped_id == -1)
          cerr << "ERROR: " << temp_cor.chr << ":" << temp_cor.begin << "-" << temp_cor.end 
               << " cannot be transferred and is either not in unmapped_map." << endl;
      }
      else {
        output << ">" << cor2_map[id].chr << ":" << cor2_map[id].begin 
               << "-" << cor2_map[id].end << temp_s.substr(pos, temp_s.size() - pos + 1) << endl;;
      }
    }
    else {
      output << temp_s << endl;
    }
  }
}

bool LiftOverGff(const vector<Coordinate>& cor1_map, const vector<Coordinate>& cor2_map, 
                 const vector<Coordinate>& unmapped_map, 
                 const string& bplib_prefix, const string& output_prefix) {
  ifstream file;
  ofstream output;
  if (!OpenFile(bplib_prefix+".gff", file)) return false;
  if (!OpenFile(output_prefix+".gff", output)) return false;
  while (!file.eof()) {
    Coordinate temp_cor;
    string temp1, temp2, temp3, temp4, temp5;
    file >> temp_cor.chr >> temp1 >> temp2 >> temp_cor.begin >> temp_cor.end >> temp3 >> temp4 >> temp5;
    if (file.eof()) break;
    temp_cor.chr = "chr" + temp_cor.chr;
    int id = SearchCoordinate(temp_cor, cor1_map);
    if (id == -1) {
      int unmapped_id = SearchCoordinate(temp_cor, unmapped_map);
      if (unmapped_id == -1)
        cerr << "ERROR: " << temp_cor.chr << ":" << temp_cor.begin << "-" << temp_cor.end 
             << " cannot be transferred and is either not in unmapped_map." << endl;
    }
    else {
      output << cor2_map[id].chr << "\t" << temp1 << "\t" << temp2 << "\t" << cor2_map[id].begin 
             << "\t" << cor2_map[id].end << "\t" << temp3 << "\t" << temp4 << "\t" << temp5 << endl;
    }
  }
}

int main (int argc, char** argv) {
  if (argc != 6) {
    cerr << "Usage: " << argv[0] << " coordinate1" << " coordinate2" 
         << " unmapped" << " bplib_prefix" << " output_prefix" << endl;
    return 1;
  }

  ifstream cor1, cor2, unmapped;
  if (!OpenFile(argv[1], cor1)) return 1;
  if (!OpenFile(argv[2], cor2)) return 1;
  if (!OpenFile(argv[3], unmapped)) return 1;
  string in_filename = argv[4], out_filename = argv[5];

  // corrdinate and get_mapped
  vector<Coordinate> cor1_map, cor2_map, unmapped_map;
  LoadCoordinateVector(cor1, cor1_map);
  LoadCoordinateVector(cor2, cor2_map);
  LoadCoordinateVector(unmapped, unmapped_map);
  // Remove unmapped positions from cor1
  RemoveUnmapped(unmapped_map, cor1_map);
  // Lifeover
  LiftOverGff(cor1_map, cor2_map, unmapped_map, in_filename, out_filename);
  LiftOverFna(cor1_map, cor2_map, unmapped_map, in_filename, out_filename,"fna");
  LiftOverFna(cor1_map, cor2_map, unmapped_map, in_filename, out_filename,"ins");
  //cout << cor_map.size() << endl;
  
  return 0;
}
