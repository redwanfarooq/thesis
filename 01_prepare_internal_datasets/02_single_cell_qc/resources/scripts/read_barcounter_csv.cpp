#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
List read_barcounter_csv(string file) {
  // read header line and extract feature names
  ifstream infile(file);
  string line;
  getline(infile, line);

  line.erase(remove(line.begin(), line.end(), '\r'), line.end()); // remove carriage return characters

  istringstream iss(line);
  string feature;
  vector<string> features;
  while(getline(iss, feature, ',')) {
    features.emplace_back(feature);
  }
  features.erase(features.begin()); // remove barcode column header

  // initialise vectors to store barcodes, row indices, values, and column pointers
  vector<string> barcodes;
  vector<long> i, x;
  vector<long> p = {0};

  // read remaining lines and extract barcodes, row indices, values, and column pointers
  while (getline(infile, line)) {
    // extract barcode
    istringstream iss(line);
    string barcode;
    getline(iss, barcode, ',');
    barcodes.emplace_back(barcode);

    // read values
    vector<long> y;
    string val;
    while(getline(iss, val, ',')) {
      y.emplace_back(stoi(val));
    }

    // extract row indices and values for non-zero entries
    for (size_t j = 0; j < y.size(); ++j) {
      if (y[j] > 0) {
        i.emplace_back(j);
        x.emplace_back(y[j]);
      }
    }

    // set next column pointer
    p.emplace_back(i.size());
  }

  return List::create(
    Named("i") = i,
    Named("p") = p,
    Named("x") = x,
    Named("features") = features,
    Named("barcodes") = barcodes
  );
}
